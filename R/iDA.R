#'  Find embedding space using iterative LDA
#'
#'  Takes scaled data and iterates between clustering using the Louvain community detection method and embedding in LDA space, then recluster in
#'  the LDA transformed data space.
#'  
#' @param data.use (data.frame) A dataframe of scaled data to find embedding for. (sample x feature)
#' @param scaled (boolean) An indicator of if the data has already been normalized and scaled.
#' @param mean.low.cutoff  (numeric) Bottom cutoff on mean for identifying variable genes, passed to function [`VariableGenes`]
#' @param mean.high.cutoff (numeric) Top cutoff on mean for identifying variable genes (passed to [`VariableGenes`])
#' @param dispersion.cutoff (numeric) Bottom cutoff on dispersion for identifying variable genes (passed to [`VariableGenes`])
#' @param dims.use (numeric) A vector of the dimensions to use in construction of the SNN
#' graph (e.g. To use the first 10 PCs, pass 1:10) (passed to [`getSNN`])
#' @param k.param (numeric) Defines k for the k-nearest neighbor algorithm (passed to [`getSNN`])
#' @param prune.SNN (numeric) Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything). Passed to [`getSNN`]
#' @param diag Diagonalize the within class scatter matrix (assume the features are independent
#' within each cluster)
#' @param set.seed (numeric or FALSE) seed random number generator before building KNN graph. (passed to [`getSNN`])
#'
#' @import irlba 
#' @import igraph
#' @return n number of dataframes for each cluster's data
#'
#'@export
#'

iDA_core <- function(data.use,
                     scaled = FALSE,
                     mean.low.cutoff = 0.1, 
                     mean.high.cutoff = 8,
                     dispersion.cutoff = 1,
                     k.param = 10,
                     prune.SNN = 1/15,
                     dims.use = 10,
                     diag = FALSE, 
                     set.seed = FALSE
){
  
 # if (scaled == FALSE){
    #normalize data by dividing by the sum of cell feature counts and then multiplying the cell counts by 10000
  #  data.use.norm <- Matrix::t((Matrix::t(data.use)/ Matrix::colSums(data.use))* 10000)
   # data.use.norm <- log1p(data.use.norm)
    
    #scale data with max 10
    #data.use.scaled <- scale(data.use.norm)
#  }

  #find variable features
  #  svd_time <- 0 
  var.features <- VariableGenes(data.use, dispersion.cutoff = dispersion.cutoff, mean.low.cutoff = mean.low.cutoff, mean.high.cutoff = mean.high.cutoff)

  #calculate svd for covariance matrix of variable_features
  
  #start_svd <- Sys.time()
  var_data <- data.use[var.features,]
  if(!is.numeric(set.seed)){
    svd <- svdr(as.matrix(var_data), k = dims.use)
  } else if (is.numeric(set.seed)){
    set.seed(set.seed)
    svd <- svdr(as.matrix(var_data), k = dims.use)
  }
  #end_svd <- Sys.time()
  #svd_time <- svd_time + (end_svd - start_svd)
  #transform data
  transformed <- svd$v
  rownames(transformed) <- colnames(var_data)
  
  
  #calculate SNN matrix for top PC's
  #louvain_time <- 0
  #start_louvain <- Sys.time()
  
  snn <- getSNN(data.use = transformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
  
  #cluster
  walktrapClusters <- igraph::cluster_walktrap(snn)
  
  
  #pick highet modularity 
  modularity <- c()
  for (i in 1:15){
    modularity <- c(modularity,  modularity(snn, igraph::cut_at(walktrapClusters, n = i)))
    
  }
  
  maxmodclust <- igraph::cut_at(walktrapClusters, n = which.max(modularity))
  clusters <- cbind(start = rep(1,dim(transformed)[1]), currentclust = maxmodclust)
  
  #end_louvain <- Sys.time()
  #louvain_time = louvain_time + (end_louvain - start_louvain)
  
  rownames(clusters) <- rownames(transformed)
  
  
  #start iterations
  i = 1
  while(sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1] < .98) {
    concordance <- sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1]
    message(paste0("iteration ", i))
    message(paste0("concordance: ", concordance))
    
    #merge data with cluster
    currentcluster <- as.data.frame(clusters[,i + 1])
    #rownames(currentcluster) <- rownames(clusters)
    merged <- merge(currentcluster, t(var_data), by = 0)
    
    #split by cluster
    splitclusters <- split_clusters(merged, merged[,2])
    
    #calculate within cluster scatter matrix
    Sw <- withinclass_scattermatrix_LDA(splitclusters = splitclusters, diag = diag)
    
    #calculate between cluster scatter matrix
    Sb <- betweenclass_scatter_matrix(splitclusters = splitclusters)
    
    #Sw-1 %*% Sb
    #   start_svd = Sys.time()
    eigenvecs <- decomposesvd(Sw, Sb, nu = length(splitclusters) - 1, set.seed = set.seed)
    #   end_svd = Sys.time()
    
    #   svd_time = svd_time + (end_svd - start_svd)
    
    
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
    #start_louvain = Sys.time()
    
    snn_transformed <- getSNN(data.use = eigenvectransformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
    
    #cluster
    walktrapClusters <- suppressWarnings(igraph::cluster_walktrap(snn_transformed))
    
    
    #pick highest modularity 
    modularity = c()
    for (j in 1:15){
      modularity <- c(modularity, modularity(snn_transformed, igraph::cut_at(walktrapClusters, n = j)))
    }
    
    maxmodclust <- igraph::cut_at(walktrapClusters, n = which.max(modularity))
    clusters <- cbind(clusters, maxmodclust)
    
    
    #end_louvain = Sys.time()
    
    #louvain_time = louvain_time + (end_louvain - start_louvain)
    i = i + 1
  }
  concordance <- sum(clusters[,dim(clusters)[2]-1] == clusters[,dim(clusters)[2]])/dim(clusters)[1]
  geneweights <- eigenvecs
  
  rownames(geneweights) <- var.features
  colnames(geneweights) <- paste("LD", 1:dim(geneweights)[2], sep = "")
  
  rownames(eigenvectransformed) <- rownames(transformed)
  colnames(eigenvectransformed) <- paste("LD", 1:dim(eigenvectransformed)[2], sep = "")
  
  message(paste0("final concordance: "))
  message(paste0(concordance))
  return(list(clusters[,dim(clusters)[2]], eigenvectransformed, geneweights, var.features))
}

