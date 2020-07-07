#' Generic method to input data to iDA 
#' 
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
#' 
#' 

setGeneric("iDA", signature=c("object"), 
function(object, ...) standardGeneric("iDA"))

#' Set method for matrix to input data to iDA 
#' 
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
#' 
#' 

setMethod("iDA", "matrix",
          function(object, ...) {
            
            iDAoutput <- iDA_core(object)
            return(iDAoutput)
            
          })



#' Method for SingleCellExperiment object to input data to iDA 
#' 
#' @param object The single cell experiment object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @importFrom  SingleCellExperiment SingleCellExperiment
#' @return SingleCellExperiment object with iDA cell weights and gene weights stored in reducedDims and cluster assignemts
#' stored in rowLabels
#' @export
#' 
#' 

setMethod("iDA", "SingleCellExperiment",
          function(object, ...) {
            counts <- assay(object)
            
            iDA_sce <- iDA(t(counts))
            
            reducedDims(object) <- list(iDA_cellweights = iDA_sce[2], iDA_geneweights = iDA_sce[3])
            rowLabels(object) <- iDA_sce[1]
            
            return(object)
          })


#' Method for Seurat object to input data to iDA 
#' 
#' @param object The single cell experiment object to run iDA on
#' @param assay The assay to take counts from
#' @param ... Additonal arguments passed to object constructors
#' @import Seurat Seurat
#' @return Seurat object with iDA cell weights and gene weights stored in object[["iDA"]] and cluster assignemts stored in rowLabels
#' @export
#' 
#' 

setMethod("iDA", "Seurat",
          function(object, assay, ...) {
            
           if (length(object[[assay]]@scale.data) == 0){
              object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
              all.genes <- rownames(object)
              object <- ScaleData(object)
              counts <- object[[assay]]@scale.data
            } else {
              counts <- object[[assay]]@scale.data
            }
            
            
            iDA_seurat <- iDA(counts, scaled = TRUE, ...)

            object[["iDA"]] <- CreateDimReducObject(embeddings = iDA_seurat[[2]], 
                                                    key = "LD", 
                                                    loadings = iDA_seurat[[3]], 
                                                    assay = assay)
  
            object@meta.data[["iDA_clust"]] <- iDA_seurat[[1]]
        
            return(object)
          })


