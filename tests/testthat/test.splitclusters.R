context("Split clusters")
library(iDA)
library(datasets)
data(iris)


test_that("split_clusters has correct dataframes in list ", {
  expect_equal(length(split_clusters(data = iris, clusterIDcol = iris$Species)), length(unique(iris$Species)))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[1]]), dim(iris[which(iris$Species == "setosa"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[2]]), dim(iris[which(iris$Species == "versicolor"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[3]]), dim(iris[which(iris$Species == "virginica"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[1]])[1] + dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[2]])[1] + dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[3]])[1], dim(iris)[1])
})


