
#'  Find Variable Genes 
#'
#'  Takes in dataframe to be split and an cluster identifier column and outputs n number of dataframes (n = number of clusters)
#'@param data.use A scaled dataframe to find the variable features. features are rows, samples are columns.
#'
#'@return variable genes
#'
#'@export

VariableGenes <- function(data.use){
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













#'  Split dataset by cluster
#'
#'  Takes in dataframe to be split and an cluster identifier column and outputs n number of dataframes (n = number of clusters)
#'@param data A dataframe to be split
#'@param clusterIDcol The column  in data with cluster identifiers
#'
#'@return n number of dataframes for each cluster's data
#'
#'@export
split_clusters <- function(data, clusterIDcol) {
  out <- split(data , f = as.factor(clusterIDcol))
}





#'  Compute each cluster's within class scatter matrix
#'
#'  Takes in the output from split_clusters() and computes the within class scatter matrix
#'@param splitclusters A list of dataframes with scaled data from each cluster (output from split_clusters())
#'@param diag if off diagonal entries in within class scatter matrix should be zeroed
#'
#'@return returns the within class scatter matrix
#'
#'@export



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


withinclass_scattermatrix_QDA <- function(splitclusters, diag = FALSE) {
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
    dataMatrix <- t(i[,4:(length(clustermeans[[k]])+3)])
    wcsm[[k]] <- (t(t(dataMatrix) - clustermeans[[k]])) %*% (t(dataMatrix) - clustermeans[[k]])
    k = k + 1
  }

  if (diag == TRUE) {
    #set off-diagonal entries to 0
    wcsm_diag <- c()
    k=1
    for (i in wcsm){
      wcsm_diag[[k]] <- diag(diag(i))
      k = k + 1
    }
    message("within cluster scatter matrix complete")
    return(wcsm_diag)
  } else {
    message("within cluster scatter matrix complete")
    return(wcsm)
  }

}





#'  Compute the between class scatter matrix
#'
#'  Takes in a list of dataframes with scaled data (output from split_clusters) and returns the between class scatter matrix
#'@param splitclusters A list of dataframes (from the output of split_clusters) with scaled data from each cluster
#'
#'@return returns the between class scatter matrix
#'
#'@export


betweenclass_scatter_matrix <- function(splitclusters){
  #calculate means vector for each cluster
  clustermeans <- c()
  k=1
  for (i in splitclusters) {
    clustermeans[[k]] <- colMeans(i[,3:(length(i))])
    k = k + 1
  }

  #calculate overallMeans for each feature
  overallMeanVector <- c()
  for (i in 1:length(clustermeans[[1]])) {
    overallMeanVector[[i]] = mean(sapply(clustermeans, function(l) l[[i]]))
  }

  #calculate each btsc matrix per cluster
  btsc <- c()
  for (i in 1:length(clustermeans)) {
  btsc[[i]] <- ((clustermeans[[i]] - overallMeanVector) %*% t(clustermeans[[i]] - overallMeanVector))
  #* length(rownames(splitclusters[[1]]))
  }

  #add all btsc's together
  Sb <- array(0L, dim(btsc[[1]]))
  k = 1
  for (i in btsc) {
    Sb <- Sb + i
    k = k + 1
  }
  return(Sb)
}



decomposesvd <- function(withinclust_sc_mat, betweenclust_sc_mat, nu = 10, set.seed = FALSE) {
  if(!is.numeric(set.seed)){
    svd <- svd(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u[,1:nu]
    return(top_eigenvectors)
  } else if(is.numeric(set.seed)){
    set.seed = set.seed
    svd <- svd(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u[,1:nu]
    return(top_eigenvectors)
  }
}


decomposeirlba <- function(withinclust_sc_mat, betweenclust_sc_mat, nu = 10, set.seed = FALSE) {
  if(!is.numeric(set.seed)){
    svd <- irlba(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u
    return(top_eigenvectors)
  } else if (is.numeric(set.seed)){
    set.seed = set.seed
    svd <- irlba(solve(withinclust_sc_mat) %*% betweenclust_sc_mat, nu)
    top_eigenvectors <- svd$u
    return(top_eigenvectors)
  }
}



#'  Find linear discriminants which best separate clusters
#'
#'  Takes in a dataframe of scaled data with cluster identification and finds the embedding space with linear discriminants to optimize the ratio of between cluster variance to within cluster variance
#'@param data a dataframe of scaled data
#'@param betweenclust_sc_mat The between class scatter matrix (output from betweenclust_sc_mat())
#'@param num_clusts The number of clusters
#'
#'@return returns the between class scatter matrix
#'
#'@export

LDAtransform <- function(data, clustident, diagSw = FALSE, set.seed = FALSE) {
  #split data into dataframes per cluster
  splitclusters <- split_clusters(data = data, clusterIDcol = clustident)

  #calculate within cluster scatter matrix
  wcsm <- withinclass_scattermatrix(splitclusters = splitclusters, diag = diagSw)

  #calculate between cluster scatter matrix
  bcsm <- betweenclass_scatter_matrix(splitclusters = splitclusters)

  #decompose Sw^-1 * Sb
  set.seed = set.seed
  feature_embedding <- decompose(withinclust_sc_mat = wcsm, betweenclust_sc_mat = bcsm, nu = length(unique(clustident)))

  #transform original data with feature_embedding
  LDAtransformed <- as.data.frame(data %*% feature_embedding)

  return(LDAtransformed)
}



#' @include seurat.R
NULL
#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. First calculate k-nearest neighbors
#' and construct the SNN graph. Then optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#' @param data.use Matrix with scaled data to find nearest neighbors
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#'
#' @importFrom igraph plot.igraph graph.adjlist
#' @importFrom Matrix sparseMatrix
#' @return
#'
#' @export
#'

ComputeSNN <- function(nn_ranked, prune) {
  .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn_ranked, prune)
}


WriteEdgeFile <- function(snn, filename, display_progress) {
  invisible(.Call('_Seurat_WriteEdgeFile', PACKAGE = 'Seurat', snn, filename, display_progress))
}



findNearestNeighbors <- function(data.use, k.param = 10, prune.SNN = 1/15, nn.eps = 0, set.seed = FALSE) {
  data.use <- as.matrix(x = data.use)
  n.obs <- nrow(x = data.use)

  if (n.obs < k.param) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
    k.param <- n.obs - 1
  }

  if (!is.numeric(set.seed)){
    my.knn <- nn2(
      data = data.use,
      k = k.param,
      searchtype = 'standard',
      eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    snn.matrix <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(x = snn.matrix) <- rownames(x = data.use)
    colnames(x = snn.matrix) <- rownames(x = data.use)
    return(snn.matrix)
  } else if (is.numeric(set.seed)){
    set.seed(set.seed)
    my.knn <- nn2(
      data = data.use,
      k = k.param,
      searchtype = 'standard',
      eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    snn.matrix <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(x = snn.matrix) <- rownames(x = data.use)
    colnames(x = snn.matrix) <- rownames(x = data.use)
    return(snn.matrix)
  }
}


#' @include seurat.R
NULL
#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. Optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#'@param SNN a matrix of shared nearest neighbors (output from findNearestNeighbors)
#'@param modularity Modularity function (1 = standard; 2 = alternative)
#'@param resolution resolution parameter for louvain clustering. Low resolution = few clusters, high resolution = many clusters
#'@param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM
#' algorithm)
#'@param n.start Number of random starts.
#'@param n.iter Maximal number of iterations per random start.
#'@param random.seed Seed of the random number generator.
#'@param edge.file.name Edge file to use as input for modularity optimizer jar.
#'@importFrom igraph plot.igraph graph.adjlist
#'@importFrom Matrix sparseMatrix
#'
#'



RunModularityClustering <- function(
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  edge.file.name = NULL
) {
  seurat.dir <- system.file(package = "Seurat")
  ModularityJarFile <- paste0(seurat.dir, "/java/ModularityOptimizer.jar")


  seurat.dir.base <- strsplit(x = seurat.dir, split = "/")[[1]]
  seurat.dir <- paste0(
    seurat.dir.base[0:(length(x = seurat.dir.base) - 1)],
    collapse = "/"
  )
  seurat.dir <- paste0(seurat.dir, "/")
  unique_ID <- sample(x = 10000:99999, size = 1)
  temp.file.location <- "/Users/taa/Downloads/"

  if(is.null(edge.file.name)) {
    edge_file <- paste0(temp.file.location, "_edge_", unique_ID, ".txt")
    #while (file.exists(edge_file)) {
    # unique_ID <- sample(x = 10000:99999, size = 1)
    #edge_file <- paste0(temp.file.location, "_edge_", unique_ID, ".txt")
    #}
    WriteEdgeFile(snn = SNN,
                  filename = edge_file,
                  display_progress = print.output)
  } else {
    if(!file.exists(edge.file.name)) {
      stop("Edge file provided doesn't exist")
    }
    edge_file <- edge.file.name
  }

  output_file <- paste0(temp.file.location, "_output_", unique_ID, ".txt")
  file.create(output_file)

  if (print.output) {
    print.output <- 1
  } else {
    print.output <- 0
  }

  command <- paste(
    "java -jar",
    shQuote(string = ModularityJarFile),
    shQuote(string = edge_file),
    shQuote(string = output_file),
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output
  )
  system(command, wait = TRUE)
  ident.use <- read.table(file = output_file, header = FALSE, sep = "\t")[, 1]
  idents <- data.frame(cells = rownames(SNN), ident = ident.use)

  file.remove(output_file)
  file.remove(edge_file)
  return(idents)
}
