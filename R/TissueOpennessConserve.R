
setClass(Class = "TissueOpennessConserve",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "TissueOpennessConserve",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        openConserveBedInput <- allparam[["openConserveBedInput"]]
        bedOutput <- allparam[["bedOutput"]]
        distrPdfOutput <- allparam[["distrPdfOutput"]]

        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            bedInput0 <- getParam(prevStep,"bedOutput")
            input(.Object)$bedInput <- bedInput0
        }

        if(!is.null(bedInput)){
            input(.Object)$bedInput <- bedInput
        }

        if(!is.null(openConserveBedInput)){
            input(.Object)$openConserveBedInput <- openConserveBedInput
        }else{
            input(.Object)$openConserveBedInput <- getRefFiles("ConserveRegion")
        }


        if(is.null(bedOutput)){
            output(.Object)$bedOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "bed")
        }else{
            output(.Object)$bedOutput <- bedOutput
        }

        if(is.null(distrPdfOutput)){
            output(.Object)$distrPdfOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "pdf")
        }else{
            output(.Object)$distrPdfOutput <- distrPdfOutput
        }

        .Object
    }
)


#' @importFrom ggplot2 ylab
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom stats median
#'
setMethod(
    f = "processing",
    signature = "TissueOpennessConserve",
    definition = function(.Object, ...){
        bedInput <- getParam(.Object,"bedInput")
        bedOutput <- getParam(.Object,"bedOutput")
        openConserveBedInput <-getParam(.Object,"openConserveBedInput")
        distrPdfOutput <- getParam(.Object,"distrPdfOutput")


        # read input data
        region <- import.bed(bedInput)
        openTable <- read.table(openConserveBedInput,header = F,sep = "\t")

        # transform openness table into GRanges format
        openBed <- openTable[,1:3]
        colnames(openBed) <- c("chrom", "start", "end")
        openValue <- openTable[,4:ncol(openTable)]
        colnames(openValue) <- seq_len(ncol(openValue))
        openRanges <- as(openBed,"GRanges")
        mcols(openRanges) <- openValue

        # find overlapped regions
        pairs <- findOverlapPairs(openRanges, region, ignore.strand = TRUE)
        openRegion <- first(pairs)
        openRegion <- as.data.frame(openRegion)
        write.table(openRegion[,c(1:3,6:ncol(openRegion))],file = bedOutput,
                    sep="\t", col.names = FALSE, row.names = FALSE)


        # draw distribution


        pdf(distrPdfOutput)
        ggplot(openRegion, aes(X2)) + geom_histogram(binwidth = 0.01) + xlab("conserve") + ylab("count")
        dev.off()


        .Object
    }
)

#' @name TissueOpennessConserve
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Conservation score of openness across tissues for the given regions
#' @description
#' Users provide region positions through a BED file.
#' This function will provide the conservation score of openness across tissues for these given regions.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' This parameter is available when the upstream step function
#' (printMap() to check the previous functions)
#' has been sucessfully called.
#' Accepted value can be the object returned by any step function or be fed by
#' \code{\%>\%} from last step function.
#' @param bedInput \code{Character} scalar.
#' The directory of input BED file for analysis.
#' @param openConserveBedInput \code{Character} scalar.
#' The openness conservation score file for analysis. The first three columns are chromosome name, start and end positions,
#' and the remaining two columns are the region name and conservation score across tissues.
#' @param bedOutput \code{Character} scalar.
#' The output file directory of merged BED files.
#' Default: NULL (generated based on bedInput)
#' @param distrPdfOutput \code{Character} scalar.
#' The openness conservation distribution figure will be provided in a PDF file.
#' @param ... Additional arguments, currently unused.
#' @details
#' We collected 201 DNase-seq/ATAC-seq samples from ENCODE and calculated their openness levels.
#' These can be downloaded and installed automatically. So users do not need to configure by themselves.
#' @return An invisible \code{\link{EnrichStep-class}}
#' object (based on \code{\link{Step-class}}) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{unzipAndMergeBed}}


#' @examples
#'
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' tissueOpennessConserve(bedInput = foregroundBedPath)


setGeneric("enrichTissueOpennessConserve",
           function(prevStep,
                    bedInput = NULL,
                    openConserveBedInput = NULL,
                    bedOutput = NULL,
                    distrPdfOutput = NULL,
                    ...) standardGeneric("enrichTissueOpennessConserve"))



#' @rdname TissueOpennessConserve
#' @aliases enrichTissueOpennessConserve
#' @export
setMethod(
    f = "enrichTissueOpennessConserve",
    signature = "Step",
    definition = function(prevStep,
                          bedInput = NULL,
                          openConserveBedInput = NULL,
                          bedOutput = NULL,
                          distrPdfOutput = NULL, ...){
        allpara <- c(list(Class = "TissueOpennessConserve",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname TissueOpennessConserve
#' @aliases tissueOpennessConserve
#' @export
tissueOpennessConserve <- function(bedInput,
                                   openConserveBedInput = NULL,
                                   bedOutput = NULL,
                                   distrPdfOutput = NULL, ...){
    allpara <- c(list(Class = "TissueOpennessConserve",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
