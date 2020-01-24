library(dplyr)
library(coloc)
library(ggplot2)
library(gridExtra)

subEnvVars  <- function(text) {
    DATA_DIR <- Sys.getenv("DATA_DIR")
    WORKING_DIR  <- Sys.getenv("WORKING_DIR")

    newText  <- sub("\\$DATA_DIR", DATA_DIR, text)
    sub("\\$WORKING_DIR", WORKING_DIR, newText)
}

verifyDirectory  <- function(dirName) {
    dir.create(dirName)
}


loadData <- function(config) {
    config$file <- subEnvVars(config$file)
    print(paste("Loading", config$file, sep=" "))
    sep = "\t"
    if ("sep" %in% names(config)) {
        sep = config$sep
    }
    data  <- read.table(config$file, sep=sep, header=T, stringsAsFactors=F)
    if ("transform_pvalue" %in% names(config)) { # when pvalue is small, may come out of zero b/c of limits of python/perl/postgres
        if (config$transform_pvalue) {
            data$pvalue = 10 ^ (-1 * data$neg_log10_pvalue)
        }
    }
    data

}

buildDetails <- function(config) {
    list(name=config$name, ncases=config$ncases, ncontrols=config$ncontrols, type="cc", nsamples=config$ncases + config$ncontrols)
}



customColoc  <- function(regions, trait1, trait2, config) {
    ## subset
    result <- NULL

    details1  <- buildDetails(trait1$config)
    details2  <- buildDetails(trait2$config)

    if (!is.null(config$log)) {
        sink(file = config$log, append = FALSE, type = c("output", "message"), split = FALSE)
    }

    for (row in 1:nrow(regions)) {

        label  <- regions[["label"]][row]

        print(paste(row, label, sep=": "))
        print(paste("Span Length:", regions$end[row] -  regions$start[row], sep=" "))

        if (regions$start[row]  == regions$end[row]) {
            print("1bp region; skipping")
            result[[label]][["result"]] <- NULL
            result[[label]][["span"]]  <- list(chr=regions$chr[row], start=regions$start[row], end=regions$end[row])
            result[[label]][["nvariants"]] <- list(trait1=1)
            result[[label]][["phenotype"]] <- regions$phenotype[row]
            result[[label]][["summary"]]  <- NULL
            next
        }

        subset1 <- trait1$data[with(trait1$data, chr==regions$chr[row] & position >= (regions$start[row]) & position < (regions$end[row])), ,drop=FALSE]
        overlap <- intersect(subset1$marker, trait2$data$marker)
        subset2 <- trait2$data[trait2$data$marker %in% overlap, ,drop=FALSE]
        row.names(subset2)  <- subset2$marker

        result[[label]][["nvariants"]] <- list(trait1=nrow(subset1),trait2=nrow(subset2),overlap=length(overlap))
        print(paste("NVARIANTS (1,2,overlap):", nrow(subset1), nrow(subset2), length(overlap)))

        subset1 <- subset1[subset1$marker %in% overlap, ,drop=FALSE]
        row.names(subset1)  <- subset1$marker

        adj  <- adjustDirectionality(subset1, subset2)

        input1 <- list(snp=adj$subset1$marker, pvalues=adj$subset1$pvalue, type=details1$type, N=details1$nsamples, s=details1$ncases/details1$nsamples)
        input2 <- list(snp=adj$subset2$marker,  pvalues=adj$subset2$pvalue, type=details2$type, N=details2$nsamples, s=details2$ncases/details2$nsamples)

        input1[["MAF"]]  <- adj$subset1$frequency
        input2[["MAF"]]  <- adj$subset2$frequency

        input1[["beta"]]  <-  adj$subset1$beta
        input1[["betavar"]]  <- adj$subset1$variance
        input2[["beta"]] <-  adj$subset2$beta
        input2[["betavar"]] <-  adj$subset2$variance

        r <- coloc.abf(input1, input2)

        print(summary(r))

        result[[label]][["result"]] <- r
        result[[label]][["span"]]  <- list(chr=regions$chr[row], start=regions$start[row], end=regions$end[row])
        result[[label]][["summary"]] <- summary(r)
        result[[label]][["phenotype"]] <- regions$phenotype[row]
    }
    # print(str(result))
    print("Done.")

    if (!is.null(config$log)) {
        sink()
    }
    result
}



adjustDirectionality  <- function(data1, data2) {
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
            }
        }

    list(subset1=data1, subset2=data2)
}

writeColocResult  <- function(colocR, trait1, trait2, fileName) {
    print(fileName)

    cnames  <- c("region_label", "phenotype", "chr", "start", "end", "length",
                 paste("nvariants", trait1, sep="_"), paste("nvariants", trait2, sep="_"), paste("ppH0", "no_causal", sep="_"),
                 paste("ppH1", trait1, "only", sep="_"), paste("ppH2", trait2, "only", sep="_"),
                 paste("ppH3", "more_than_one_shared_causal", sep="_"), paste("ppH4", "one_shared_causal", sep="_"))
    df <- NULL
    for (label in names(colocR)) {

        if (colocR[[label]]$nvariants$trait1 == 1) {
            print(paste("WARNING: no result for", label, "as region span was 1bp."))
            next
        }
        row  <-  list(label, colocR[[label]]$phenotype, colocR[[label]]$span$chr, colocR[[label]]$span$start, colocR[[label]]$span$end,
                      colocR[[label]]$span$end - colocR[[label]]$span$start,
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
    manhattan.plot(data$chr, data$position, data$pvalue, annotate=annotation)
}


saveResult  <- function(regions, result, trait1, trait2, config) {
    print("Writing summary")
    summary  <- writeColocResult(result, trait1$config$name, trait2$config$name, paste(config$output_path, paste(config$comparison, "result_summary.txt", sep="_"), sep="/"))

    print("Saving result workspace")
    fileName  <- paste(config$comparison, "result.RData.gz", sep="_")
    save(result, summary, file=paste(config$output_path, fileName, sep="/"), compress=TRUE)
    summary
}

generateGraphics  <- function(regions, result, trait1, trait2, config) {
    print("Outputting graphics")
    for (row in 1:nrow(regions)) {
        plotRegion  <- plotRegion  <- list(chr=regions[row, "chr"], start=regions[row,"start"], end=regions[row,"end"], label=regions[row, "label"], pp=result[[row]]$result$summary[["PP.H4.abf"]])
        print(plotRegion$label)
        ##region, trait1, trait2, type, filename) {
        fileName  <- paste(config$comparison, plotRegion$label, "manhattan.pdf", sep="_")
        plotRegionManhattan(plotRegion, trait1, trait2, config$comparison, paste(config$output_path, fileName, sep="/"))
        fileName  <- paste(config$comparison, plotRegion$label, "shared_cvps.pdf", sep="_")
        plotRegionPriors(plotRegion$label, result, data1, filePrefix, paste(config$output_path, fileName, sep="/"))
    }
}


plotRegionPriors  <- function(region, result, data, comparison, filename) {
    r  <- result[[region]]$result$results
    subset  <- data[data$marker %in% r$snp, ]

    row.names(subset)  <- subset$marker
    row.names(r)  <- r$snp

    positions  <- subset[, "pos", drop=FALSE]
    r <- merge(r, positions, by=0)
    row.names(r) <- r$snp

     plot <- ggplot(data=r, aes(x=pos, y=SNP.PP.H4, group=1)) +
         ggtitle(paste(comparison, region, "Prior Prob Shared Causal Variant", sep=" - ")) +
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

plotRegionManhattan  <- function(region, trait1, trait2, comparison, filename) {



    subset1 <- trait1$data[with(trait1$data, chr==region$chr & position >= (region$start) & position < (region$end)), ,drop=FALSE]
    subset1  <- subset1[with(subset1, order(position)), ]
    print("subset1")
    print(dim(subset1))

    overlap <- intersect(subset1$marker, trait2$data$marker)

    subset1.filtered <- subset1[subset1$marker %in% overlap, ,drop=FALSE]
    row.names(subset1.filtered)  <- subset1.filtered$marker
    positions  <- subset1.filtered[,"position", drop=FALSE]


    subset2 <- trait2$data[trait2$data$marker %in% overlap, ,drop=FALSE]
    row.names(subset2)  <- subset2$marker

    subset2 <- merge(subset2, positions, by=0)
    row.names(subset2)  <- subset2$marker
    print("subset2")
    print(dim(subset2))

    ylim.upper  <- max(max(subset1$neg_log10_pvalue), max(subset2$neg_log10_pvalue))

    ## pushViewport(viewport(layout = grid.layout(3, 1)))

    plot1 <- ggplot(data=subset1, aes(x=position, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(comparison, region$label, trait1, paste("N", nrow(subset1), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
        ylim(0, ylim.upper) +
        ylab("-log10 p") +
        geom_point(aes(colour = cut(neg_log10_pvalue, c(-Inf, 5, Inf)))) +
        scale_color_manual(name = "-log10 p",
                     values = c("(-Inf,5]" = "black",
                                  "(5, Inf]" = "red"),
                     labels = c("<= 1e-5", "> 1e-5"))

        print("plot1")

     plot1.filtered <- ggplot(data=subset1.filtered, aes(x=position, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(comparison, region$label, trait1, "filtered", paste("N", nrow(subset1.filtered), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
         ylim(0, ylim.upper) +
         ylab("-log10 p") +
        geom_point(aes(colour = cut(neg_log10_pvalue, c(-Inf, 5, Inf)))) +
        scale_color_manual(name = "-log10 p",
                     values = c("(-Inf,5]" = "black",
                                  "(5, Inf]" = "red"),
                     labels = c("<= 1e-5", "> 1e-5"))


      plot2 <- ggplot(data=subset2, aes(x=position, y=neg_log10_pvalue, group=1)) +
        ggtitle(paste(comparison, region$label, trait2, paste("N", nrow(subset2), sep="="), paste("pp", region$pp, sep="="), sep=" - ")) +
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
