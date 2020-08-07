#' Find Variable Genes
#'
#' Identify rows in a data frame with high sclaed dispersion. This rule was
#' taken from the Seurat package.
#'
#' @param data.use (data.frame) features are rows, samples are columns. Rows must be named.
#' @param dispersion.cutoff (numeric) rows returned will have scaled dispersion higher provided cutoff
#' @param mean.low.cutoff (numeric) rows returned will have average higher than this cutoff
#' @param mean.high.cutoff (numeric) rows returned will have average lower than this cutoff
#'
#' @importFrom stats sd var
#'
#' @return (character) a list of row names with high dispersion rows
#' @export
VariableGenes <- function(data.use,
                          dispersion.cutoff,
                          mean.low.cutoff,
                          mean.high.cutoff) {
  #calculate logged means and VMR
  ExpMeans <- apply(data.use, 1, FUN = function(x) log(mean(exp(x) - 1) + 1))
  dispersions <- apply(data.use, 1, FUN = function(x) {log(var(exp(x) - 1) / mean( exp(x) - 1))})

  names(x = ExpMeans) <- names(x = dispersions) <- rownames(x = data.use)
  dispersions[is.na(x = dispersions)] <- 0
  ExpMeans[is.na(x = ExpMeans)] <- 0

  #create bins to scale mean and dispersions
  num.bin <- 20
  data.x.bin <- cut(x = ExpMeans, breaks = num.bin)

  names(x = data.x.bin) <- names(x = ExpMeans)

  mean.y <- tapply(X = dispersions, INDEX = data.x.bin, FUN = mean)
  sd.y <- tapply(X = dispersions, INDEX = data.x.bin, FUN = sd)

  #scale dispersions
  scaled.dispersions <- (dispersions - mean.y[as.numeric(x = data.x.bin)]) /
    sd.y[as.numeric(x = data.x.bin)]
  names(x = scaled.dispersions) <- names(x = ExpMeans)

  #find variable features
  var.features <- names(dispersions[scaled.dispersions > dispersion.cutoff &
                                      ExpMeans > mean.low.cutoff &
                                      ExpMeans < mean.high.cutoff])
  return(var.features)
}

#' Split dataset by cluster label
#'
#' Takes a data.frame to and the name of the column cluster identifier column and outputs n dataframes (n = number of clusters)
#' @param data (data.frame) A dataframe to be split
#' @param clusterIDcol (character) The name of the column in data with cluster identifiers
#'
#' @return n number of dataframes for each cluster's data
split_clusters <- function(data, clusterIDcol) {
  split(data , f = as.factor(clusterIDcol))
}

#' Compute each cluster's within class scatter matrix
#'
#' Takes in the output from split_clusters() and computes the within class scatter matrix
#'
#' @param splitclusters A list of dataframes with scaled data from each cluster (output from split_clusters())
#' @param diag if off diagonal entries in within class scatter matrix should be zeroed
#'
#' @return returns the within class scatter matrix
withinclass_scattermatrix_LDA <- function(splitclusters, diag = FALSE) {
  #calculate means vector for each cluster
  clustermeans <- c()
  k=1
  for (i in splitclusters) {
    clustermeans[[k]] <- colMeans(i[,3:(length(i))])
    k = k + 1
  }

  #calculate within class scatter matrix for each cluster
  wcsm <- c()
  k= 1
  for (i in splitclusters) {
    dataMatrix <- t(i[,3:(length(i))])
    wcsm[[k]] <- (t(t(dataMatrix) - clustermeans[[k]])) %*% (t(dataMatrix) - clustermeans[[k]])
    k = k + 1
  }

  #add all within class scatter matrices together
  #Sw <- matrix(0L, nrow = nrow(wcsm[[1]]), ncol = ncol(wcsm[[1]]))
  Sw <- array(0L, dim(wcsm[[1]]))
  k = 1
  list <- sapply(splitclusters, function(l) l[1])
  n_obs <- sum(lengths(list))
  for (i in wcsm) {
    Sw <- Sw + ((dim(splitclusters[[k]])[1])/n_obs) * i
    k = k + 1
  }

  if (diag == TRUE) {
    #set off-diagonal entries to 0
    Sw <- diag(diag(Sw))
  }
  return(Sw)
}

#' Compute each cluster's within class scatter matrix (for QDA)
#'
#' Takes in the output from split_clusters() and computes the within class scatter matrix
#'
#' @param splitclusters A list of dataframes with scaled data from each cluster (output from split_clusters())
#' @param diag if off diagonal entries in within class scatter matrix should be zeroed
#'
#' @return returns the within class scatter matrix
withinclass_scattermatrix_QDA <- function(splitclusters, diag = FALSE) {
  #calculate means vector for each cluster
  clustermeans <- c()
  k <- 1
  for (i in splitclusters) {
    clustermeans[[k]] <- colMeans(i[,3:(length(i))])
    k = k + 1
  }

  #calculate within class scatter matrix for each cluster
  wcsm <- c()
  k= 1
  for (i in splitclusters) {
    dataMatrix <- t(i[,4:(length(clustermeans[[k]])+3)])
    wcsm[[k]] <- (t(t(dataMatrix) - clustermeans[[k]])) %*% (t(dataMatrix) - clustermeans[[k]])
    k <- k + 1
  }

  if (diag == TRUE) {
    #set off-diagonal entries to 0
    wcsm_diag <- c()
    k <- 1
    for (i in wcsm){
      wcsm_diag[[k]] <- diag(diag(i))
      k <-  k + 1
    }
    message("within cluster scatter matrix complete")
    return(wcsm_diag)
  } else {
    message("within cluster scatter matrix complete")
    return(wcsm)
  }
}

#' Compute the between class scatter matrix
#'
#' Takes in a list of dataframes with scaled data (output from split_clusters) and returns the between class scatter matrix
#' @param splitclusters A list of dataframes (from the output of split_clusters) with scaled data from each cluster
#'
#' @return returns the between class scatter matrix
#'
betweenclass_scatter_matrix <- function(splitclusters){
  #calculate means vector for each cluster
  clustermeans <- c()
  k <- 1
  for (i in splitclusters) {
    clustermeans[[k]] <- colMeans(i[,3:(length(i))])
    k <- k + 1
  }

  #calculate overallMeans for each feature
  overallMeanVector <- c()
  for (i in 1:length(clustermeans[[1]])) {
    overallMeanVector[[i]] <- mean(sapply(clustermeans, function(l) l[[i]]))
  }

  #calculate each btsc matrix per cluster
  btsc <- c()

  for (i in 1:length(clustermeans)) {
    btsc[[i]] <- ((clustermeans[[i]] - unlist(overallMeanVector)) %*% t(clustermeans[[i]] - unlist(overallMeanVector)))
    #* length(rownames(splitclusters[[1]]))
  }

  #add all btsc's together
  Sb <- array(0L, dim(btsc[[1]]))
  k <- 1
  for (i in btsc) {
    Sb <- Sb + i
    k <- k + 1
  }
  return(Sb)
}

## TODO: This is not using set.seed properly
decomposesvd <- function(withinclust_sc_mat,
                         betweenclust_sc_mat,
                         nu = 10, set.seed = FALSE) {
  if (!is.numeric(set.seed)) {
    svd <- svd(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u[,1:nu]
    top_eigenvalues <- svd$d[1:nu]
    return(list(top_eigenvectors, top_eigenvalues))
  } else if(is.numeric(set.seed)) {
    set.seed <- set.seed
    svd <- svd(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u[,1:nu]
    top_eigenvalues <- svd$d[1:nu]
    return(list(top_eigenvectors, top_eigenvalues))
  }
}

## TODO: This is not using set.seed properly
decomposeirlba <- function(withinclust_sc_mat, betweenclust_sc_mat, nu = 10, set.seed = FALSE) {
  if(!is.numeric(set.seed)){
    svd <- irlba(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u
    return(top_eigenvectors)
  } else if (is.numeric(set.seed)){
    set.seed <- set.seed
    svd <- irlba(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u
    return(top_eigenvectors)
  }
}

#' Cluster Determination
#'
#' Calculate k-nearest neighbors and construct a shared nearest neighbor (SNN) graph.
#'
#' @param data.use (matrix) Matrix with scaled data to find nearest neighbors
#' @param k.param (numeric) Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN (numeric) Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param set.seed (numeric or FALSE) seed random number generator before building KNN graph
#'
#' @import igraph
#' @import scran
#'
#' @return The SNN graph (igraph object)
getSNN <- function(data.use,
                   k.param = 10,
                   prune.SNN = 1/15,
                   set.seed = FALSE) {
  data.use <- as.matrix(data.use)
  n.obs <- nrow(x = data.use)

  if (n.obs < k.param) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
    k.param <- n.obs - 1
  }

  ## TODO: refactor this to avoid code duplication
  if (!is.numeric(set.seed)){

    SNN_igraph <- scran::buildKNNGraph(
      data.use,
      k = k.param,
      transposed = TRUE)
    snn.matrix <- similarity(
      SNN_igraph,
      method = "jaccard")

    snn.matrix[snn.matrix < 1/15] <- 0
    rownames(x = snn.matrix) <- rownames(x = data.use)
    colnames(x = snn.matrix) <- rownames(x = data.use)

    snn.graph <- graph_from_adjacency_matrix(snn.matrix, weighted = TRUE, mode = "undirected")
    return(snn.graph)

  } else if (is.numeric(set.seed)){
    set.seed(set.seed)

    SNN_igraph <- scran::buildKNNGraph(
      data.use,
      k = k.param,
      transposed = TRUE)
    snn.matrix <- similarity(
      SNN_igraph,
      method = "jaccard")

    snn.matrix[snn.matrix < 1/15] <- 0
    rownames(x = snn.matrix) <- rownames(x = data.use)
    colnames(x = snn.matrix) <- rownames(x = data.use)

    snn.graph <- graph_from_adjacency_matrix(snn.matrix, weighted = TRUE, mode = "undirected")
    return(snn.graph)
  }
}

# #' Cluster Determination
# #'
# #' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
# #' optimization based clustering algorithm. Optimize the modularity function to
# #' determine clusters. For a full description of the algorithms, see Waltman and
# #' van Eck (2013) \emph{The European Physical Journal B}.
# #'
# #'
# #'@param SNN a matrix of shared nearest neighbors (output from getSNN)
# #'@param resolution resolution parameter for louvain clustering. Low resolution = few clusters, high resolution = many clusters
# #'@param random.seed Seed of the random number generator.
#
# #'@importFrom NetworkToolbox louvain
# getLouvain <- function(SNN, resolution = 1, random.seed = 0){
#   if (!is.numeric(set.seed)){
#     louvain_clusts <- louvain(SNN, gamma = resolution)
#     idents <- louvain_clusts$community
#     return(idents)
#   }  else if (is.numeric(set.seed)){
#     set.seed(set.seed)
#
#     louvain_clusts <- louvain(SNN, gamma = resolution)
#     idents <- louvain_clusts$community
#     return(idents)
#   }
#
# }
