#' Generic method to input data to iDA 
#' 
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
#' 
#' 

setGeneric("iDA", signature=c("object"), 
function(object, columns=NULL, ...) standardGeneric("iDA"))


#' Method for SingleCellExperiment object to input data to iDA 
#' 
#' @param object The single cell experiment object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return SingleCellExperiment object with iDA cell weights and gene weights stored in reducedDims and cluster assignemts
#' stored in rowLabels
#' @export
#' 
#' 

setMethod("iDA", "SingleCellExperiment",
          function(object, columns, ...) {
            counts <- assay(object)
            
            iDA_sce <- iDA(t(counts))
            
            reducedDims(object) <- list(iDA_cellweights = iDA_sce[2], iDA_geneweights = iDA_sce[3])
            rowLabels(object) <- iDA(iDA_clusters = iDA_sce[1])
            
            return(object)
          })


