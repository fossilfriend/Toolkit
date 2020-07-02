library(dplyr)
library(coloc)
library(ggplot2)
library(gridExtra)
library(reshape)
library(mvc)


subEnvVars  <- function(text) {
    DATA_DIR <- Sys.getenv("DATA_DIR")
    WORKING_DIR  <- Sys.getenv("WORKING_DIR")

    newText  <- sub("\\$DATA_DIR", DATA_DIR, text)
    sub("\\$WORKING_DIR", WORKING_DIR, newText)
}

verifyDirectory  <- function(dirName) {
    dir.create(dirName)
}

  adjFreq <- function(x) {
    if (x > 0.5) {
      return (1.0 - x)
	}
    return(x)
	      }

    adjBetaSign <- function(x) {
      if (x > 0.5) {
        return(-1)
	  }
      return (1)
		 }

loadData <- function(config) {
    config$file <- subEnvVars(config$file)
    print(paste("Loading", config$file, sep=" "))
    sep = "\t"
    if ("sep" %in% names(config)) {
        sep = config$sep
    }
    data  <- read.table(config$file, sep=sep, header=T, stringsAsFactors=F)
    if ("transform_pvalue" %in% names(config)) { 
      ## when pvalue is small, may come out of zero b/c of limits of python/perl/postgres
        if (config$transform_pvalue) {
            data$pvalue = 10 ^ (-1 * data$neg_log10_pvalue)
        }
    }
    ## adjust maf/beta (make sure using maf/switch beta sign if necessary)
    data$beta <- sapply(data$frequency, adjBetaSign) * data$beta
    data$frequency <- sapply(data$frequency, adjFreq)
    
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

        # adj  <- adjustDirectionality(subset1, subset2) # can't do this anymore b/c we adjusted beta when loading to supply MAF; no way to know which allele is the test allele

        input1 <- list(snp=subset1$marker, pvalues=subset1$pvalue, type=details1$type, N=details1$nsamples, s=details1$ncases/details1$nsamples)
        input2 <- list(snp=subset2$marker,  pvalues=subset2$pvalue, type=details2$type, N=details2$nsamples, s=details2$ncases/details2$nsamples)

        input1[["MAF"]]  <- subset1$frequency
        input2[["MAF"]]  <- subset2$frequency

        input1[["beta"]]  <-  subset1$beta
        input1[["betavar"]]  <- subset1$variance
        input2[["beta"]] <-  subset2$beta
        input2[["betavar"]] <-  subset2$variance

        r <- coloc.abf(input1, input2)

        print(summary(r))

        result[[label]][["result"]] <- r
        result[[label]][["span"]]  <- list(chr=regions$chr[row], start=regions$start[row], end=regions$end[row])
        result[[label]][["summary"]] <- summary(r)
        result[[label]][["phenotype"]] <- regions$phenotype[row]
        result[[label]][["coloc_abf"]] <- r
    }
    # print(str(result))
    print("Done.")

    if (!is.null(config$log)) {
        sink()
    }
    result
}

customSensitivity  <- function(result, rule, config) {
    for (i in seq(1,length(result))) {
        tad  <- names(result)[i]
        print(paste(i, tad, sep=": "))
        abf  <-  result[tad][[1]]$coloc_abf ## result data frame
        sensitivity(abf, rule)
    }    
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
        row  <-  list(label,
                      colocR[[label]]$phenotype,
                      colocR[[label]]$span$chr,
                      colocR[[label]]$span$start,
                      colocR[[label]]$span$end,
                      colocR[[label]]$span$end - colocR[[label]]$span$start,
                      colocR[[label]]$nvariants$trait1,
                      colocR[[label]]$nvariants$trait2,
                      colocR[[label]]$result$summary[["PP.H0.abf"]],
                      colocR[[label]]$result$summary[["PP.H1.abf"]],
                      colocR[[label]]$result$summary[["PP.H2.abf"]],
                      colocR[[label]]$result$summary[["PP.H3.abf"]],
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



saveResult  <- function(regions, result, trait1, trait2, config, isSensitivityAnalysis) {
    print("Writing summary")
    summary  <- writeColocResult(result, trait1$config$name, trait2$config$name, paste(config$output_path, paste(config$comparison, "result_summary.txt", sep="_"), sep="/"))

    print("Saving result workspace")
    fileName  <- paste(config$comparison, "result.RData.gz", sep="_")
    outputPath <- if (isSensitivityAnalysis) paste0(config$output_path, "/sensitivity") else config$output_path
    save(result, summary, file=paste(outputPath, fileName, sep="/"), compress=TRUE)
    summary
}

generateGraphics  <- function(result, data, config, isSensitivityAnalysis) {
    ## data should have at least marker, position fields
    print("Outputting graphics")
    trait1  <- config$conditions[[1]]$phenotype
    trait2  <- config$conditions[[2]]$phenotype
    print(paste("Comparison:", trait1, "/", trait2, sep=" "))
          
    for (i in seq(1,length(result))) {
        tad  <- names(result)[i]

        print(paste(i, tad, sep=": "))
        strokeVariant <- strsplit(tad, "_")[[1]][1]
        print(strokeVariant)

        print("EXTRACTING RESULTS")
        r  <-  result[tad][[1]]$result$results ## result data frame

        s <- result[tad][[1]]$result$summary ## result summary table

        d.pp1 <- s[3]
        d.pp2 <- s[4]
        d.pp4 <- s[6]

        print("WARNING: Subsetting data")
        data.s <- data[data$marker %in% r$snp, c("marker", "position")] # extract variants in tad region

        df  <- merge(data.s, r, by.x="marker", by.y="snp") # map position to result
        df.pvalues <- df # do it again, so we can have a dataframe for doing the manhattan plots

        df$pp1 <- exp(df$lABF.df1 - logsum(df$lABF.df1)) # calc per variant pp1
        df$pp2 <- exp(df$lABF.df2 - logsum(df$lABF.df2)) # calc per variant pp2

        ## restructure dataframes so can generate one plot per variable
        df <- melt(df[,c("marker","position","pp1","pp2","SNP.PP.H4")], id.vars=c("marker","position"))
        df.pvalues <- melt(df.pvalues[,c("marker","position","pvalues.df1", "pvalues.df2" )], id.vars=c("marker","position"))

        df.pvalues$value <- -1 * log(df.pvalues$value, 10) # -log10p

        df$variable <- sub("pp1", paste0("pp1 (", trait1, ") = ", d.pp1) ,df$variable)
        df$variable <- sub("pp2", paste0("pp2 (", trait2, ") = ", d.pp2) ,df$variable)
        df$variable <- sub("SNP.PP.H4", paste0("pp4 (Both) = ", d.pp4) ,df$variable)

        df$variable <- sub("pp1", paste0("pp1 (", trait1, ") = ", d.pp1) ,df$variable)
        df$variable <- sub("pp2", paste0("pp2 (", trait2, ") = ", d.pp2) ,df$variable)
        df$variable <- sub("SNP.PP.H4", paste0("pp4 (Both) = ", d.pp4) ,df$variable)

        df.pvalues$variable <- sub("pvalues.df1", paste0("pvalues (", trait1 ,")"), df.pvalues$variable)
        df.pvalues$variable <- sub("pvalues.df2", paste0("pvalues (", trait2 ,")"), df.pvalues$variable)

        ## find the top 10 snps
        print("WARNING: Finding top 10 variants")
        df.ord = df[order(df$value, decreasing=TRUE), ]
        snps <- unique(df.ord$marker)[1:10]

        df$label <- ifelse(df$marker %in% snps, df$marker,"")

        df$label <- ifelse(df$marker == strokeVariant, "", df$label)
        df$is_stroke_variant <- ifelse(df$marker == strokeVariant, df$marker,"")

        df.pvalues$label <- ifelse(df.pvalues$marker %in% snps, df.pvalues$marker,"")

        df.pvalues$label <- ifelse(df.pvalues$marker == strokeVariant, "", df.pvalues$label)
        df.pvalues$is_stroke_variant <- ifelse(df.pvalues$marker == strokeVariant, df.pvalues$marker, "")

        ttl <- paste0(trait1, ' & ', trait2, ' ', tad)
        chr  <- strsplit(strsplit(tad, "_")[[1]][3], ":")[[1]][1]

        pPlot <- ggplot(df, aes_string(x="position",y="value")) +
            geom_point(data=subset(df,label==""),size=1.5) +
            geom_point(data=subset(df,label!=""),col="red",size=1.5) +
            geom_text(aes_string(label="label"),hjust=-0.1,vjust=0.5,size=2.5,col="red") +
            geom_point(data=subset(df,is_stroke_variant!=""),col="blue",size=1.5) +
            geom_text(aes_string(label="is_stroke_variant"),hjust=-0.1,vjust=0.5,size=2.5,col="blue") +
            facet_grid(variable ~ .) +
            theme(legend.position="none") + ylab("Posterior probability") +
            ggtitle(ttl) + theme_bw()

        mPlot <- ggplot(df.pvalues, aes_string(x="position",y="value")) +
            geom_point(data=subset(df.pvalues,label==""),size=1.5) +
            geom_point(data=subset(df.pvalues,label!=""),col="red",size=1.5) +
            geom_text(aes_string(label="label"),hjust=-0.1,vjust=0.5,size=2.5,col="red") +
            geom_point(data=subset(df.pvalues,is_stroke_variant!=""),col="blue",size=1.5) +
            geom_text(aes_string(label="is_stroke_variant"),hjust=-0.1,vjust=0.5,size=2.5,col="blue") +
            facet_grid(variable ~ .) +
            theme(legend.position="none") + ylab("-log10 pvalue") +
            ggtitle(ttl) + theme_bw()

        fileName  <- paste0(gsub(":", "_", tad), ".pdf")
        outputPath <- if (isSensitivityAnalysis) paste0(config$output_path, "/sensitivity") else config$output_path
        ggsave(paste(outputPath, fileName , sep='/'), plot = marrangeGrob(grobs = list(PP=pPlot, M=mPlot), nrow=2, ncol=1), device="pdf")

    }
}

writeSingleVariantResults  <- function(results, data, config, isSensitivityAnalysis) {
    output <- NULL
  
    for (i in seq(1,length(results))) {
        tad  <- names(results)[i]

        
        print(paste0("Extracting single variants results for TAD: ", tad))
        marker <- strsplit(tad, "_")[[1]][1]

        r  <-  results[tad][[1]]$result$results

        s <- results[tad][[1]]$result$summary

        d.pp1 <- s[3]
        d.pp2 <- s[4]
        d.pp3 <- s[5]
        d.pp4 <- s[6]

        data.s <- data[data$marker %in% r$snp, c("marker", "pvalue", "position")]

        ##print(dim(r))
        ##print(dim(data.s))

        df  <- merge(data.s, r, by.x="marker", by.y="snp")

        ##print(head(df))


        df$pp1 <- exp(df$lABF.df1 - logsum(df$lABF.df1))
        df$pp2 <- exp(df$lABF.df2 - logsum(df$lABF.df2))
        df$tad  <- tad
        df$tad_marker  <- marker

        filename = paste(config$output_path, paste0(gsub(":", "_", tad), "_single_variant_results.txt"), sep="/")
        write.table(df, filename, quote=FALSE,row.names=FALSE, sep="\t")
        output  <- rbind(output, df)
    }

    outputPath <- if (isSensitivityAnalysis) paste0(config$output_path, "/sensitivity") else config$output_path
    filename = paste(outputPath, "complete_single_variant_results.txt", sep="/")
    write.table(output, filename, quote=FALSE,row.names=FALSE, sep="\t")
    output
}

