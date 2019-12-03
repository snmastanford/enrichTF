
setClass(Class = "TissueOpennessSpecificity",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "TissueOpennessSpecificity",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        openBedInput <- allparam[["openBedInput"]]
        sampleTxtInput <- allparam[["sampleTxtInput"]]
        bedOutput <- allparam[["bedOutput"]]
        distPdfOutput <- allparam[["distPdfOutput"]]
        heatmapPdfOutput <- allparam[["heatmapPdfOutput"]]
        sampleTxtOutput <- allparam[["sampleTxtOutput"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            bedInput0 <- getParam(prevStep,"bedOutput")
            input(.Object)$bedInput <- bedInput0
        }

        if(!is.null(bedInput)){
            input(.Object)$bedInput <- bedInput
        }

        if(!is.null(openBedInput)){
            input(.Object)$openBedInput <- openBedInput
        }else{
            input(.Object)$openBedInput <- getRefFiles("OpenRegion")
        }

        if(!is.null(sampleTxtInput)){
            input(.Object)$sampleTxtInput <- sampleTxtInput
        }else{
            input(.Object)$sampleTxtInput <- getRefFiles("SampleName")
        }


        if(is.null(bedOutput)){
            output(.Object)$bedOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "bed")
        }else{
            output(.Object)$bedOutput <- bedOutput
        }

        if(is.null(distPdfOutput)){
            output(.Object)$distPdfOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "dist.pdf")
        }else{
            output(.Object)$distPdfOutput <- distPdfOutput
        }

        if(is.null(heatmapPdfOutput)){
            output(.Object)$heatmapPdfOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "heat.pdf")
            output(.Object)$heatmapDataOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "heat.pdf.Rdata")
        }else{
            output(.Object)$heatmapPdfOutput <- heatmapPdfOutput
            output(.Object)$heatmapDataOutput <- paste0(heatmapPdfOutput,".Rdata")
        }

        if(is.null(sampleTxtOutput)){
            output(.Object)$sampleTxtOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "txt")
        }else{
            output(.Object)$sampleTxtOutput <- sampleTxtOutput
        }


        .Object
    }
)

#' @importFrom ggpubr ggarrange
#' @importFrom  heatmap3 heatmap3
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 xlab
#'
#'

setMethod(
    f = "processing",
    signature = "TissueOpennessSpecificity",
    definition = function(.Object,...){
        bedInput <- getParam(.Object,"bedInput")
        bedOutput <- getParam(.Object,"bedOutput")
        openBedInput <-getParam(.Object,"openBedInput")
        sampleTxtInput <-getParam(.Object,"sampleTxtInput")
        distPdfOutput <- getParam(.Object,"distPdfOutput")
        heatmapPdfOutput <- getParam(.Object,"heatmapPdfOutput")
        sampleTxtOutput <- getParam(.Object,"sampleTxtOutput")
        heatmapDataOutput <- getParam(.Object,"heatmapDataOutput")

        message("read input data")
        region <- import.bed(bedInput)
        openTable <- read.table(openBedInput,header = F,sep = "\t")
        spname <- read.table(sampleTxtInput, header = F, sep = "\t")

        message("transform open table into GRanges format")
        openBed <- openTable[,1:3]
        colnames(openBed) <- c("chrom", "start", "end")
        openValue <- openTable[,4:ncol(openTable)]
        colnames(openValue) <- as.character(seq_len(ncol(openValue)))
        openRanges <- as(openBed,"GRanges")
        mcols(openRanges) <- openValue

        message("find overlapped region")
        pairs <- findOverlapPairs(openRanges, region, ignore.strand = TRUE)
        openRegion <- first(pairs)
        openValue <- mcols(openRegion)
        colnames(openValue) <- as.character(seq_len(ncol(openValue)))
        openRegion <- as.data.frame(openRegion)
        write.table(openRegion[,c(1:3,6:ncol(openRegion))],file = bedOutput,
                    sep="\t", col.names = FALSE, row.names = FALSE)



        message("median openness level of each tissue sample")

        rs<-lapply(as.list(openValue), median)
        allidx<-order(unlist(rs),decreasing = TRUE)
        showspname<-spname[,3:4]
        showspname<-cbind(seq_len(nrow(showspname)),showspname,unlist(rs))
        showspname <- showspname[allidx,]
        colnames(showspname)<-c("Index","Tissue / Cell Type", "ENCODE", "Median")

        write.table(showspname, file = sampleTxtOutput,col.names = TRUE, row.names = FALSE,quote = FALSE, sep = '\t')

         message("draw distribution")

        idx <- allidx

        plt<-lapply(idx, function(x){
            v <- openValue[[x]]
            ggplot(data.frame(v=v),aes(x=v)) +
                geom_histogram(binwidth = 0.1) +
                geom_vline(xintercept = median(v)) +
                annotate("text", x = median(v),
                         y = 50, color="white",
                         size=2 ,label = paste("median:", median(v))) + xlab(spname[x,2])
        })

        plt[["nrow"]] <- ceiling(length(idx)/2)
        plt[["ncol"]] <- 2

        pdf(distPdfOutput)
        do.call(what = ggarrange,args = plt)
        dev.off()

        message("draw heatmap")

        openheat <- openValue
        rownames(openheat) <- as.character(seq_len(nrow(openValue)))
        colnames(openheat) <- as.character(spname[,3])

        heatmapData <- as.matrix(openheat)

        pdf(heatmapPdfOutput)
        heatmap3(heatmapData, useRaster = TRUE)
        dev.off()

        save(heatmapData,file = heatmapDataOutput)





        .Object
    }
)


#' @name TissueOpennessSpecificity
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Tissue-specific openness scores for the given regions
#' @description
#' Users provide region positions through a BED file.
#' This function will provide tissue-specific openness degree of the given regions.
#' Median value of Openness level,  openness distribution and clustering result (heatmap)
#' based on tissues and regions will be caculated.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' This parameter is available when the upstream step functions
#' (printMap() to check the previous functions)
#' have been sucessfully called.
#' Accepted value can be the object returned by any step function or be fed by
#' \code{\%>\%} from last step function.
#' @param bedInput \code{Character} scalar.
#' The directory of input BED file for analysis.
#' @param openBedInput \code{Character} scalar.
#' The openness level file of different tissues for analysis. The first three columns are chromosome name, start and end positions,
#' The remaining columns are the openness values for each tissue.
#' The order of tissues should be consistent with the file provided by sampleTxtInput.
#' @param sampleTxtInput \code{Character} scalar.
#' The tissue annotation of samples in the file provided by openBedInput.
#' There are 4 columns seperated by tab. The first column indicates the order.
#' The second column is the tissue annotation in detail.
#' The third column is the tissue name.
#' The forth column is the accession number in its original database such as ENCODE project.
#' @param bedOutput \code{Character} scalar.
#' The output file directory of merged BED files.
#' Default: NULL (generated based on bedInput)
#' @param distPdfOutput \code{Character} scalar.
#' The openness level distribution figure for each tissue will be provided in PDF file.
#' It is shown in descending order.
#' @param heatmapPdfOutput \code{Character} scalar.
#' The openness level hierarchical clustering heatmap base on regions and tissues will be provided in this PDF file.
#' The corresponding heatmap data will be stored in the Rdata format at the same directory.
#' @param  sampleTxtOutput \code{Character} scalar.
#' In this file, there are five columns seperated by tab.
#' The first four columns are the same with sampleTxtInput:
#' The first column is the order number.
#' The second column is the tissue detail information.
#' The third column is the tissue name.
#' The forth column is the code from source project such as ENCODE project.
#' The last column is the median value of the openness level for each tissue.
#' The table is in descending order of last column
#' @param ... Additional arguments, currently unused.
#' @details
#' We collect 201 DNase-seq or ATAC-seq samples from ENCODE and calculate their openness level value.
#' They can be downloaded and installed automatically. So users do not need to configure by themselves.
#' @return An invisible \code{\link{EnrichStep-class}}
#' object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{unzipAndMergeBed}}


#' @examples
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' tissueOpennessSpecificity(bedInput = foregroundBedPath)
#'



setGeneric("enrichTissueOpennessSpecificity",
           function(prevStep,
                    bedInput = NULL,
                    openBedInput = NULL,
                    sampleTxtInput = NULL,
                    bedOutput = NULL,
                    distPdfOutput = NULL,
                    heatmapPdfOutput = NULL,
                    sampleTxtOutput = NULL,
                    ...) standardGeneric("enrichTissueOpennessSpecificity"))



#' @rdname TissueOpennessSpecificity
#' @aliases enrichTissueOpennessSpecificity
#' @export
setMethod(
    f = "enrichTissueOpennessSpecificity",
    signature = "Step",
    definition = function(prevStep,
                          bedInput = NULL,
                          openBedInput = NULL,
                          sampleTxtInput = NULL,
                          bedOutput = NULL,
                          distPdfOutput = NULL,
                          heatmapPdfOutput = NULL,
                          sampleTxtOutput = NULL, ...){
        allpara <- c(list(Class = "TissueOpennessSpecificity",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname TissueOpennessSpecificity
#' @aliases tissueOpennessSpecificity
#' @export
tissueOpennessSpecificity <- function(bedInput,
                                      openBedInput = NULL,
                                      sampleTxtInput = NULL,
                                      bedOutput = NULL,
                                      distPdfOutput = NULL,
                                      heatmapPdfOutput = NULL,
                                      sampleTxtOutput = NULL, ...){
    allpara <- c(list(Class = "TissueOpennessSpecificity",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
