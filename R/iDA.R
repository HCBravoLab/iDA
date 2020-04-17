#'  Find embedding space using iterative LDA
#'
#'  Takes scaled data and iterates between clustering using the Louvain community detection method and embedding in LDA space, then recluster in
#'  the LDA transformed data space.
#'@param data.use A dataframe of scaled data to find embedding for. (sample x feature)
#'@param mean.low.cutoff  Bottom cutoff on x-axis for identifying variable genes
#'@param mean.high.cutoff Top cutoff on x-axis for identifying variable genes
#'@param dispersion.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param dims.use A vector of the dimensions to use in construction of the SNN
#' graph (e.g. To use the first 10 PCs, pass 1:10)
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param diag Diagonalize the within class scatter matrix (assume the features are independent
#' within each cluster)
#' @param decompose.fxn Which method to use for computing eigenvectors. Available methods are:
#' \itemize{
#' \item{svd:}{ use svd() function to find singular value decomposition [default]}
#' \item{irlba:}{use irlba() function: Fast and memory efficient methods for truncated singular value decomposition
#'  and principal components analysis of large sparse and dense matrices.}
#' }
#' @param resolution The resolution parameter to use for Louvain modularity clustering
#' @param reduction.type Type of reduction to use in the iterative step. Available methods are:
#' \itemize{
#' \item{LDA:}{ use linear discriminant analysis to find discriminants between classes (assumes equal covariance of clusters) [default]}
#' \item{QDA:}{use quadratic discriminant analysis to find discriminants between classes (allows for different covariances between clusters)
#'
#'
#' @import irlba 
#' @import igraph
#' @import FNN
#' @import RANN 
#' @import scran
#'@return n number of dataframes for each cluster's data
#'
#'@export


iDA <- function(data.use,  
                mean.low.cutoff = 0.1, 
                mean.high.cutoff = 8,
                dispersion.cutoff = 1,
                k.param = 10,
                prune.SNN = 1/15,
                nn.eps = 0,
                reduction.type = "LDA",
                dims.use = 10,
                diag = FALSE, 
                set.seed = FALSE,
                resolution = 1.0
                ){


  #find variable features
  svd_time = 0 
  var.features <- VariableGenes(data.use, dispersion.cutoff = dispersion.cutoff, mean.low.cutoff = mean.low.cutoff, mean.high.cutoff = mean.high.cutoff)

  #calculate svd for covariance matrix of variable_features

start_svd = Sys.time()
    var_data <- data.use[var.features,]
    if(!is.numeric(set.seed)){
    svd <- svdr(as.matrix(var_data), k = dims.use)
    } else if (is.numeric(set.seed)){
      set.seed(set.seed)
      svd <- svdr(as.matrix(var_data), k = dims.use)
    }
end_svd <- Sys.time()
svd_time <- svd_time + (end_svd - start_svd)
  #transform data
    transformed <- svd$v
    rownames(transformed) <- colnames(var_data)


    #calculate SNN matrix for top PC's
louvain_time <- 0
start_louvain <- Sys.time()

      snn <- getSNN(data.use = transformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN, nn.eps = 0)
      
    #cluster
      louvainClusters = getLouvain(snn, resolution = resolution, random.seed = set.seed)
      
      clusters <- c(rep(1, dim(snn)[1]))
      clusters <- cbind(clusters, louvainClusters)
end_louvain <- Sys.time()

      rownames(clusters) <- rownames(transformed)


louvain_time = louvain_time + (end_louvain - start_louvain)

      i = 1
    #start iterations
    while(sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1] < .98) {
      concordance = sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1]
      message(paste0("iteration ", i))
      message(paste0("concordance: ", concordance))

      #merge data with cluster
        currentcluster <- as.data.frame(clusters[,i + 1])
        rownames(currentcluster) <- rownames(clusters)
        merged <- merge(currentcluster, t(var_data), by = 0)

      #split by cluster
        splitclusters <- split_clusters(merged, merged[,2])

        
        #calculate within cluster scatter matrix
                   Sw <- withinclass_scattermatrix_LDA(splitclusters = splitclusters, diag = diag)
         
                 #calculate between cluster scatter matrix
                   Sb <- betweenclass_scatter_matrix(splitclusters = splitclusters)
         
        #Sw-1 %*% Sb
         start_svd = Sys.time()
                   eigenvecs <- decomposesvd(Sw, Sb, nu = length(splitclusters) - 1, set.seed = set.seed)
         end_svd = Sys.time()
         
         svd_time = svd_time + (end_svd - start_svd)
        
        
        
        
#       if (reduction.type == "LDA") {
#         #calculate within cluster scatter matrix
#           Sw <- withinclass_scattermatrix_LDA(splitclusters = splitclusters, diag = diag)
# 
#         #calculate between cluster scatter matrix
#           Sb <- betweenclass_scatter_matrix(splitclusters = splitclusters)
# 
#         #Sw-1 %*% Sb
#           
# start_svd = Sys.time()
#           eigenvecs <- decomposesvd(Sw, Sb, nu = length(splitclusters) - 1, set.seed = set.seed)
# end_svd = Sys.time()
# 
# svd_time = svd_time + (end_svd - start_svd)
        
      # } else if (reduction.type == "QDA") {
      #   Sw = withinclass_scattermatrix_QDA(splitclusters = splitclusters, diag = diag)
      #   Sb = betweenclass_scatter_matrix(splitclusters = splitclusters)
      #   invSwSb = c()
      #   k=1
      #   for (i in Sw){
      #     invSwSb[[k]] = decomposesvd(i,Sb, nu = 1)
      #     k = k + 1
      #   }
      #   eigenvecs = as.matrix(invSwSb, ncol = length(Sw))
      # }

      #transform data
        eigenvectransformed <- t(var_data) %*% eigenvecs

      #calculate SNN matrix for top LDs
start_louvain = Sys.time()
        snn_transformed <- getSNN(data.use = eigenvectransformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN, nn.eps = 0)

      #cluster
        currentclust <- getLouvain(snn_transformed, resolution = resolution, random.seed = set.seed)
        clusters <- cbind(clusters, currentclust)
        i = i + 1

end_louvain = Sys.time()

louvain_time = louvain_time + (end_louvain - start_louvain)
    }
  concordance = sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1]
  geneweights = eigenvecs
  message(paste0("final concordance: "))
  message(paste0(concordance))
  return(list(clusters[,dim(clusters)[2]], transformed, geneweights, louvain_time, svd_time))
}








