#'  Find embedding space using iterative LDA
#'
#'  Takes scaled data and iterates between clustering using the Louvain community detection method and embedding in LDA space, then recluster in
#'  the LDA transformed data space.
#'
#' @param data.use (data.frame) A dataframe of scaled data to find embedding for. (sample x feature)
#' @param scaled (boolean) An indicator of if the data has already been normalized and scaled.
#' @param var.Features Which method to use when finding variable features
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
#' @param c.param (numeric) Defines the number of desired clusters to be found in the embedding
#'
#' @import irlba
#' @import igraph
#' @import plyr
#' @return n number of dataframes for each cluster's data
#'
#'@export

iDA_core <- function(data.use,
                     NormCounts = NULL, 
                     scaled = FALSE,
                     var.Features = "scran",
                     mean.low.cutoff = 0.1, 
                     mean.high.cutoff = 8,
                     dispersion.cutoff = 1,
                     k.param = 10,
                     prune.SNN = 1/15,
                     dims.use = 10,
                     diag = TRUE, 
                     set.seed = FALSE, 
                     c.param = NULL,
                     cluster.method = "walktrap"
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
  if (var.Features == "scran") {
    stats <- scran::modelGeneVar(NormCounts)
    if (dim(data.use)[1] < 3000){
      var.features <- rownames(data.use)
    } else {
      var.features <- scran::getTopHVGs(stats, n = 3000)
    }
    
  } else if (var.Features == "disp") {
    if (is.null(NormCounts)) {
      # variance stabilizing transformation using Deseq2
      NormCounts <- varianceStabilizingTransformation(data.use)
    }
    var.features <- VariableGenes(NormCounts, dispersion.cutoff = dispersion.cutoff, mean.low.cutoff = mean.low.cutoff, mean.high.cutoff = mean.high.cutoff)
  }
  
  if(length(var.features) == 0){
    stop("No variable features found.")
  }
  
  message(length(var.features), " variable features found. \n")
  
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
  
  #cluster
  if(cluster.method == "louvain") {
    snn <- getSNN(data.use = transformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
    louvainClusters <- getLouvain(SNN = snn, set.seed = set.seed)
    clusters <- cbind(start = rep(1,dim(transformed)[1]), currentclust = louvainClusters)
    
  } else if (cluster.method == "kmeans"){
    kmeansclusters <- kmeans(transformed, centers = c.param)
    clusters <- cbind(start = rep(1,dim(transformed)[1]), currentclust = kmeansclusters$cluster)
    
  } else if (cluster.method == "walktrap"){
    snn <- getSNN(data.use = transformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
    if(!is.numeric(set.seed)){
      walktrapClusters <- suppressWarnings(igraph::cluster_walktrap(snn))
    } else if (is.numeric(set.seed)){
      set.seed(set.seed)
      walktrapClusters <- suppressWarnings(igraph::cluster_walktrap(snn))
    }
    
    #pick highest modularity
    if (is.null(c.param)){
      modularity <- c(0)
      for (i in 2:15){
        modularity <- c(modularity,  modularity(snn, suppressWarnings(igraph::cut_at(walktrapClusters, n = i))))
      }
      
      maxmodclust <- igraph::cut_at(walktrapClusters, n = which.max(modularity))
      clusters <- cbind(start = rep(1,dim(transformed)[1]), currentclust = maxmodclust)
    } else if (is.numeric(c.param)) {
      maxmodclust <- igraph::cut_at(walktrapClusters, n = c.param)
      clusters <- cbind(start = rep(1,dim(transformed)[1]), currentclust = maxmodclust)
    } else {
      stop("Invalid c.param")
    }
  }
  
  
  rownames(clusters) <- rownames(transformed)
  
  concordance <- adjustedRandIndex(clusters[,(dim(clusters)[2]-1)], clusters[,(dim(clusters)[2])])
  
  #start iterations
  i = 1        
  while(concordance < .98) {
    if(i > 1){
      message(paste0("iteration ", i-1))
      message(paste0("concordance: ", concordance))
    }
    
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
    
    #   start_svd = Sys.time()
    eigenvecs <- decomposesvd(Sw, Sb, nu = length(splitclusters) - 1, set.seed = set.seed)[["eigenvecs"]]
    
    #transform data
    eigenvectransformed <- t(var_data) %*% eigenvecs
    
    #calculate SNN matrix for top LDs
    if (cluster.method == "louvain") {
      snn <- getSNN(data.use = eigenvectransformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
      louvainClusters <- getLouvain(SNN = as.matrix(snn), set.seed = set.seed)
      clusters <- cbind(clusters, currentclust = louvainClusters)
      
    } else if (cluster.method == "kmeans"){
      kmeansclusters <- kmeans(eigenvectransformed, centers = c.param)
      clusters <- cbind(clusters, currentclust = kmeansclusters$cluster)
      
    } else if (cluster.method == "walktrap"){
      snn_transformed <- getSNN(data.use = eigenvectransformed, set.seed = set.seed, k.param = k.param, prune.SNN = prune.SNN)
      #cluster
      walktrapClusters <- suppressWarnings(igraph::cluster_walktrap(snn_transformed))
      
      #pick highest modularity 
      if (is.null(c.param)){
        modularity <- c(0)
        for (j in 2:15){
          modularity <- c(modularity,  modularity(snn_transformed, suppressWarnings(igraph::cut_at(walktrapClusters, n = j))))
        }
        maxmodclust <- igraph::cut_at(walktrapClusters, n = which.max(modularity))
        clusters <- cbind(clusters, currentclust = maxmodclust)
      } else if (is.numeric(c.param)) {
        maxmodclust <- igraph::cut_at(walktrapClusters, n = c.param)
        clusters <- cbind(clusters, currentclust = maxmodclust)
      } else {
        stop("Invalid c.param")
      }
    }
    
    concordance <- adjustedRandIndex(clusters[,(dim(clusters)[2]-1)], clusters[,(dim(clusters)[2])])
    
    #end_louvain = Sys.time()
    
    #louvain_time = louvain_time + (end_louvain - start_louvain)
    i = i + 1
  }
  
  geneweights <- as.data.frame(eigenvecs)
  
  rownames(geneweights) <- var.features
  colnames(geneweights) <- paste("LD", 1:dim(geneweights)[2], sep = "")
  
  rownames(eigenvectransformed) <- rownames(transformed)
  colnames(eigenvectransformed) <- paste("LD", 1:dim(eigenvectransformed)[2], sep = "")
  
  message(paste0("final concordance: "))
  message(paste0(concordance))
  return(list(clusters[,(dim(clusters)[2])], eigenvectransformed, geneweights, var.features))
  
}
