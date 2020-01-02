library(dplyr)
library(RPostgreSQL)
library(coloc)
library(ggplot2)
library(gridExtra)


loadRegions <- function(filename) {
    read.table(filename, header=TRUE, sep=",", stringsAsFactors = FALSE)
}

fetchVariantPosition  <- function(refsnp, user, pwd, dbname, host, port) {
    sth <- dbConnect(PostgreSQL(), user=user, password=pwd, dbname=dbname, host=host, port=port)
    SQL  <- "SELECT replace(chromosome, 'chr', '') || ':' || position::text AS position FROM NIAGADS.Variant
WHERE variant_id IN (SELECT find_variant_by_refsnp($1::text)) LIMIT 1"

    qh <- dbSendQuery(sth, SQL, c(refsnp))
    r  <- fetch(qh, n = -1)
    dbDisconnect(sth)
    r$position
}


fetchGWASFromDB <- function(accession,user,pwd,dbname,host,port) {
    sth <- dbConnect(PostgreSQL(), user=user, password=pwd, dbname=dbname, host=host, port=port)

    ##on.exit(dbDisconnect(sth))

    start <- Sys.time()
    print(paste("Fetching", accession))

    SQL2  <- "WITH datasets AS (SELECT $1::text AS track)
SELECT DISTINCT
CASE WHEN split_part(metaseq_id, ':', 1) = 'X' THEN 23
WHEN split_part(metaseq_id, ':', 1) = 'Y' THEN 24
WHEN split_part(metaseq_id, ':', 1) = 'M' THEN 25
ELSE split_part(metaseq_id, ':', 1)::integer END AS chr,

split_part(metaseq_id, ':', 2)::integer AS pos,

allele AS testAllele,

metaseq_id AS variant,
CASE WHEN r.source_id LIKE 'rs%' THEN r.source_id  ELSE NULL END AS marker,

max(r.neg_log10_pvalue) OVER w AS neg_log10_pvalue,
first_value(r.pvalue_display::text) OVER w AS display_pvalue,
first_value(r.frequency) OVER w AS frequency,

first_value(r.restricted_stats->>'beta') OVER w AS beta,
first_value(power((r.restricted_stats->>'freq_se')::numeric, 2)) OVER w AS variance

FROM
Results.VariantGWAS r,
Study.ProtocolAppNode pan,
Datasets
WHERE pan.source_id = datasets.track
AND r.protocol_app_node_id = pan.protocol_app_node_id

WINDOW w as (PARTITION BY metaseq_id)"

    qh <- dbSendQuery(sth, SQL2, c(accession))
    data <- fetch(qh, n = -1) # extract all rows


    dbDisconnect(sth)

    print(Sys.time() - start)

    print("Unlogging p-value")
    data$pvalue <- 10 ^ (-1 * data$neg_log10_pvalue)

    data
}


assembleGeneRegions <- function(genes, phenotype) {
    print("Fetching gene regions")
    regions  <- NULL
    for (g in genes) {
        print(g)
        loc  <- fetchGeneLocation(g)
        result  <- list(g, paste(g, phenotype, sep="_"), loc$chr, loc$start, loc$end)
        regions  <- rbind.data.frame(regions, result, stringsAsFactors=FALSE, make.row.names=FALSE)
    }

    colnames(regions) <- c("gene", "label", "chr", "start", "end")

    regions
}


fetchGeneLocation  <- function(gene, user, pwd,dbname, host, port) {
    sth <- dbConnect(PostgreSQL(), user=user, password=pwd, dbname=dbname, host=host, port=port)
    SQL  <- "SELECT replace(chromosome, 'chr', '')::integer AS chr, location_start AS start, location_end AS end
FROM CBIL.GeneAttributes WHERE gene_symbol = trim(replace($1::text,'*','')) ORDER BY transcript_count DESC LIMIT 1"

    qh <- dbSendQuery(sth, SQL, c(gene))
    geneLoc  <- fetch(qh, n = -1)
    dbDisconnect(sth)
    geneLoc
}

loadStroke <- function(datasetIndex) {
    fileName <- paste(STROKE$path, STROKE$file[datasetIndex], sep="/")

    print(paste("Loading megastroke data:", fileName))
    start <- Sys.time()
    data <- read.table(fileName,  header=TRUE, stringsAsFactors= FALSE) # nrows=50000,
    print(Sys.time() - start)

    print(dim(data))

    colnames(data) <- c("marker", "testallele", "allele2", "frequency", "beta", "variance", "pvalue")
    data$variance  <- data$variance^2
    data$testallele  <- toupper(data$testallele)
    data$allele2  <- toupper(data$allele2)
    data$neg_log10_pvalue  <- -log10(data$pvalue)

    data
}

generateAnnotation <- function(data, regions, flanking){
    print("generating annotation")
    annotation <- rep(1, length(data$pvalue))
    for (row in 1:nrow(regions)) {
        ## print(as.character(regions$variant)[row])
        annotation[with(data, chr==regions$chr[row] & pos >= (regions$start[row] - flanking) & pos < (regions$end[row] + flanking))] <- row
    }
    print(paste("Done. Generated annotation for", nrow(regions), "regions."))
    annotation <- factor(annotation, levels=1:nrow(regions), labels=regions[["label"]])
    annotation
}


doColoc <- function(regions, data1, details1, data2, details2, flanking=0, logFile) {
    ## subset
    result <- NULL


    if (!is.null(logFile)) {
        sink(file = logFile, append = FALSE, type = c("output", "message"), split = FALSE)
    }

    MAF <- TRUE
    if (is.na(data1[1,"frequency"]) || is.na(data2[1,"frequency"])) {
        MAF <- FALSE
    }


    BETA  <- TRUE
    if (is.na(data1[1,"beta"]) || is.na(data2[1,"beta"])) {
        BETA <- FALSE
    }

    for (row in 1:nrow(regions)) {

        label  <- regions[["label"]][row]

        print(paste(row, label, sep=": "))
        print(paste(regions$chr[row], regions$start[row] - flanking, regions$end[row] + flanking))

        if (regions$start[row] - flanking  == regions$end[row] + flanking) {
            print("1bp region; skipping")
            result[[label]][["result"]] <- NULL
            result[[label]][["span"]]  <- list(chr=regions$chr[row], start=regions$start[row] - flanking, end=regions$end[row] + flanking)
            result[[label]][["nvariants"]] <- list(trait1=1)
            result[[label]][["summary"]]  <- NULL
            next
        }

        subset1 <- data1[with(data1, chr==regions$chr[row] & pos >= (regions$start[row] - flanking) & pos < (regions$end[row] + flanking)), ,drop=FALSE]
        overlap <- intersect(subset1$marker, data2$marker)
        subset2 <- data2[data2$marker %in% overlap, ,drop=FALSE]
        row.names(subset2)  <- subset2$marker

        result[[label]][["nvariants"]] <- list(trait1=nrow(subset1),trait2=nrow(subset2),overlap=length(overlap))
        print(paste("NVARIANTS (1,2,overlap):", nrow(subset1), nrow(subset2), length(overlap)))

        subset1 <- subset1[subset1$marker %in% overlap, ,drop=FALSE]
        row.names(subset1)  <- subset1$marker

        adj  <- adjustDirectionality(subset1, subset2, !(MAF && BETA))


        trait1 <- list(snp=adj$subset1$marker, pvalues=adj$subset1$pvalue, type=details1$type, N=details1$nsamples, s=details1$ncases/details1$nsamples)
        trait2 <- list(snp=adj$subset2$marker,  pvalues=adj$subset2$pvalue, type=details2$type, N=details2$nsamples, s=details2$ncases/details2$nsamples)

        if (MAF) {
            trait1[["MAF"]]  <- adj$subset1$maf
            trait2[["MAF"]]  <- adj$subset2$maf
        }

        if (BETA) {
            trait1[["beta"]]  <-  adj$subset1$beta
            trait1[["betavar"]]  <- adj$subset1$variance
            trait2[["beta"]] <-  adj$subset2$beta
            trait2[["betavar"]] <-  adj$subset2$variance
        }

        r <- coloc.abf(trait1, trait2)

        print(summary(r))

        result[[label]][["result"]] <- r
        result[[label]][["span"]]  <- list(chr=regions$chr[row], start=regions$start[row] - flanking, end=regions$end[row] + flanking)
        result[[label]][["summary"]] <- summary(r)
    }
    # print(str(result))
    print("Done.")

    if (!is.null(logFile)) {
        sink()
    }
    result
}



adjustDirectionality  <- function(data1, data2, skipAdjustment) {
    if (!skipAdjustment) {
        for (marker in row.names(data1)) {
            testAllele.data1  <- data1[marker, "testallele"]
            testAllele.data2  <- data2[marker, "testallele"]
            freq.data1  <- data1[marker, "frequency"]
            freq.data2  <- data2[marker, "frequency"]

            if ((freq.data1 <= 0.40 && freq.data2 >= 0.5) || (freq.data1 >= 0.5 && freq.data2 <= 0.4)) {
                print(paste("WARNING: frequency mismatch found for marker: ", marker, " - ", freq.data1 , "/", freq.data2, " - ", testAllele.data1 , "/", testAllele.data2, sep=""))
            }

            if (testAllele.data1 != testAllele.data2) {
                print(paste("WARNING: Allele mismatch found for marker: ", marker, " - ", testAllele.data1 , "/", testAllele.data2, sep=""))
                beta.data2  <- data2[marker, "beta"]

                data2[marker, "testallele"]  <- testAllele.data1
                data2[marker, "beta"]  <- 1 * beta.data2
                data2[marker, "frequency"]  <- 1.0 - freq.data2
                data2[marker, "maf"]  <- freq.data2
            }
        }
    }
    list(subset1=data1, subset2=data2)
}

writeColocResult  <- function(colocR, trait1, trait2, fileName) {

    cnames  <- c("tag_region_label", "chr", "start", "end",
                 paste("nvariants", trait1, sep="_"), paste("nvariants", trait2, sep="_"), paste("ppH0", "no_causal", sep="_"),
                 paste("ppH1", trait1, "only", sep="_"), paste("ppH2", trait2, "only", sep="_"),
                 paste("ppH3", "more_than_one_shared_causal", sep="_"), paste("ppH4", "one_shared_causal", sep="_"))
    df <- NULL
    for (label in names(colocR)) {

        if (colocR[[label]]$nvariants$trait1 == 1) {
            print(paste("WARNING: no result for", label, "as region span was 1bp."))
            next
        }
        row  <-  list(label, colocR[[label]]$span$chr, colocR[[label]]$span$start, colocR[[label]]$span$end,
                  colocR[[label]]$nvariants$trait1, colocR[[label]]$nvariants$trait2,
                  colocR[[label]]$result$summary[["PP.H0.abf"]], colocR[[label]]$result$summary[["PP.H1.abf"]],
                  colocR[[label]]$result$summary[["PP.H2.abf"]], colocR[[label]]$result$summary[["PP.H3.abf"]],
                  colocR[[label]]$result$summary[["PP.H4.abf"]])

        df <- rbind.data.frame(df, row, stringsAsFactors=FALSE, make.row.names=FALSE)
    }

    colnames(df)  <- cnames
    df  <- df[with(df, order(-ppH4_one_shared_causal)), ]
    write.table(df, file=fileName, sep="\t", quote=FALSE, row.names=FALSE)
    df
}

extractPP <- function(colocR, hindex) {
    pp <- vector()
    ppH  <- paste("PP.H", hindex, ".abf", sep="")
    for (region in names(colocR)) {
        pp  <- c(pp, colocR[[region]]$result$summary[[ppH]])
    }
    names(pp)  <- names(colocR)
    pp
}



mPlot <- function(data, annotation) {
    manhattan.plot(data$chr, data$pos, data$pvalue, annotate=annotation)
}


outputResult  <- function(regions, resul, flanking, data1, data2, trait1, trait2, filePrefix) {
    print("Writing summary")
    summary  <- writeColocResult(result, trait1, trait2, paste(filePrefix, "coloc-summary.txt", sep="-"))

    print("Outputting graphics")
    for (row in 1:nrow(regions)) {
        plotRegion  <- plotRegion  <- list(chr=regions[row, "chr"], start=regions[row,"start"] - flanking, end=regions[row,"end"] + flanking, label=regions[row, "label"], pp=result[[row]]$result$summary[["PP.H4.abf"]])
        print(plotRegion$label)
        plotRegionManhattan(plotRegion, data1, data2, trait1, trait2, filePrefix, paste(filePrefix, plotRegion$label, "manhattan.pdf", sep="-"))
        plotRegionPriors(plotRegion$label, result, data1, filePrefix, paste(filePrefix, plotRegion$label, "shared-causal-variant-priors.pdf", sep="-"))
    }

    print("Saving result workspace")
    save(result, summary, file=paste(filePrefix, "result.RData.gz", sep="-"), compress=TRUE)


}


plotRegionPriors  <- function(region, result, data, type, filename) {
    r  <- result[[region]]$result$results
    subset  <- data[data$marker %in% r$snp, ]

    row.names(subset)  <- subset$marker
    row.names(r)  <- r$snp

    positions  <- subset[, "pos", drop=FALSE]
    r <- merge(r, positions, by=0)
    row.names(r) <- r$snp

     plot <- ggplot(data=r, aes(x=pos, y=SNP.PP.H4, group=1)) +
         ggtitle(paste(type, region, "Prior Prob Shared Causal Variant", sep=" - ")) +
         geom_point(aes(colour = cut(SNP.PP.H4, c(-Inf, .1, Inf)))) +
         scale_color_manual(name = "SNP.PP.H4",
                     values = c("(-Inf,0.1]" = "black",
                                  "(0.1, Inf]" = "red"),
                     labels = c("<= 10%", "> 10%")) +
         geom_text(aes(label=ifelse(SNP.PP.H4>0.1,as.character(snp),'')),hjust=0,vjust=0)


    pdf(filename)
    print(plot)
    dev.off()

}

plotRegionManhattan  <- function(region, data1, data2, trait1, trait2, type, filename) {

    subset1 <- data1[with(data1, chr==region$chr & pos >= (region$start) & pos < (region$end)), ,drop=FALSE]
    subset1  <- subset1[with(subset1, order(pos)), ]

    overlap <- intersect(subset1$marker, data2$marker)

    subset1.filtered <- subset1[subset1$marker %in% overlap, ,drop=FALSE]
    row.names(subset1.filtered)  <- subset1.filtered$marker
    positions  <- subset1.filtered[,"pos", drop=FALSE]


    subset2 <- data2[data2$marker %in% overlap, ,drop=FALSE]
    row.names(subset2)  <- subset2$marker

    subset2 <- merge(subset2, positions, by=0)
    row.names(subset2)  <- subset2$marker


    ylim.upper  <- max(max(subset1$neg_log10_pvalue), max(subset2$neg_log10_pvalue))

    ## pushViewport(viewport(layout = grid.layout(3, 1)))

    plot1 <- ggplot(data=subset1, aes(x=pos, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(type, region$label, trait1, paste("N", nrow(subset1), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
        ylim(0, ylim.upper) +
        ylab("-log10 p") +
        geom_point(aes(colour = cut(neg_log10_pvalue, c(-Inf, 5, Inf)))) +
        scale_color_manual(name = "-log10 p",
                     values = c("(-Inf,5]" = "black",
                                  "(5, Inf]" = "red"),
                     labels = c("<= 1e-5", "> 1e-5"))


     plot1.filtered <- ggplot(data=subset1.filtered, aes(x=pos, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(type, region$label, trait1, "filtered", paste("N", nrow(subset1.filtered), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
         ylim(0, ylim.upper) +
         ylab("-log10 p") +
        geom_point(aes(colour = cut(neg_log10_pvalue, c(-Inf, 5, Inf)))) +
        scale_color_manual(name = "-log10 p",
                     values = c("(-Inf,5]" = "black",
                                  "(5, Inf]" = "red"),
                     labels = c("<= 1e-5", "> 1e-5"))


      plot2 <- ggplot(data=subset2, aes(x=pos, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(type, region$label, trait2, paste("N", nrow(subset2), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
          ylim(0, ylim.upper) +
          ylab("-log10 p") +
        geom_point(aes(colour = cut(neg_log10_pvalue, c(-Inf, 5, Inf)))) +
        scale_color_manual(name = "-log10 p",
                     values = c("(-Inf,5]" = "black",
                                  "(5, Inf]" = "red"),
                     labels = c("<= 1e-5", "> 1e-5"))

    pdf(filename, width=12, height=12)
    grid.arrange(plot1, plot1.filtered, plot2, ncol=1, nrow=3)
    ##print(plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    ##print(plot1.filtered, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    ## print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    dev.off()
    ##print(plot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

    ##plot(subset1$pos, subset1$neg_log10_pvalue, col="red",
     ##    main=paste(region$label, trait1), xlab="Position", ylab="-log10 pvalue")

#    lines(subset2$pos, subset$neg_log10_pvalue, type = "o", col = "blue")



}
