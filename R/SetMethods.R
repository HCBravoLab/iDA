#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature=c("object"),
           function(object, ...) standardGeneric("iDA"))

#' Set method for matrix to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setMethod("iDA", "matrix",
          function(object, ...) {
            iDAoutput <- iDA_core(object, ...)
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
setMethod("iDA", "SingleCellExperiment",
          function(object, ...) {
              if (!('logcounts' %in% names(assays(object)))){
                    counts <- assay(object, "counts")
                    libsizes <- colSums(counts)
                    size.factors <- libsizes/mean(libsizes)
                    logcounts(object) <- log2(t(t(counts)/size.factors) + 1)
              }


              normcounts <-  logcounts(object)
              iDA_sce <- iDA(normcounts, scaled = TRUE, ...)
              reducedDims(object) <- list(iDAcellweights = iDA_sce[[2]])
              colLabels(object) <- list(iDAclusters = iDA_sce[[1]])
              rowData(object[iDA_sce[[4]],]) <- list(iDAgeneweights = iDA_sce[[3]])

              return(object)
          })

#' Method for Seurat object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param assay The assay to take counts from
#' @param ... Additional arguments passed to object constructors
#' @import Seurat Seurat
#' @return Seurat object with iDA cell weights and gene weights stored in object[["iDA"]] and cluster assignments stored in rowLabels
#' @export

setMethod("iDA", "Seurat",
          function(object, assay, ...) {

            if (length(object[[assay]]@scale.data) == 0){
              object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
              all.genes <- rownames(object)
              object <- ScaleData(object)
              counts <- object[[assay]]@scale.data
              NormCounts <- object[[assay]]@data
            } else {
              counts <- object[[assay]]@scale.data
              NormCounts <- object[[assay]]@data
            }

            iDA_seurat <- iDA(counts, scaled = TRUE, ...)
            object[["iDA"]] <- CreateDimReducObject(embeddings = iDA_seurat[[2]],
                                                    key = "LD_",
                                                    loadings = iDA_seurat[[3]],
                                                    assay = assay)

            object@meta.data[["iDA_clust"]] <- iDA_seurat[[1]]
            object@reductions$iDA@stdev <- iDA_seurat[[5]]
            return(object)
          })
