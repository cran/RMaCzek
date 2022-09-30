## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


#' @title Calculate the distance matrix between clusters.
#' @description Calculate the distance matrix for a czek_matrix with clustering result or a data set with its clustering labels.
#' @param x A data set or a matrix with class czek_matrix.
#' @param y If x is the data set, y is the cluster label.
#' @param distfun Specifies which distance function should be used.
#' @param dist_method Four linkage criteria: single, complete, average and SSD.
#'
#' @return A distance matrix.
#' @export
#'
#' @examples
#' # Clustering Result on czek_matrix
#' x = czek_matrix(iris[,-5], cluster = TRUE, num_cluster = 3)
#' dist_czek = cluster_dist(x)
#' plot(czek_matrix(dist_czek))
#'
#' # Clustering Result on a Data Set with Clustering Labels
#' dist_data = cluster_dist(x = iris[,-5], y = iris$Species)
#' plot(czek_matrix(dist_data))
#'
cluster_dist = function(x, y, distfun = dist, dist_method = "average"){
  if(inherits(x,"czek_matrix")){
    if(is.null(attr(x, "cluster"))){
      stop("Please try a clustered czek_matrix.")
    }

    dist_mat = x
    num_cluster = attr(x, "num_cluster")
    y = attr(x, "cluster_res")

    clusters = list()
    for (i in 1:num_cluster) {
      clusters[[i]] = which(y == i)
    }
    names = paste("cluster", 1:num_cluster)
  }else if(inherits(x,"matrix") | inherits(x,"data.frame")){
    if(nrow(x) != length(y)){
      stop("The lengths of x and y do not match.")
    }
    dist_mat = as.matrix(do.call(distfun, list(x)))
    names = unique(y)
    num_cluster = length(names)

    clusters = list()
    for (i in names) {
      clusters[[i]] = which(y == i)
    }
  }

  if(num_cluster <= 2){
    warning("The results would not be plotted by Czekanowski's diagram because the number of clusters is less than 3.")
  }

  cluster_dist = matrix(NA, nrow = num_cluster, ncol = num_cluster)
  for (i in 1:num_cluster) {
    for (j in 1:num_cluster) {
      if(i == j){
        cluster_dist[i,j] = 0
      }else{
        if(dist_method == "single"){
          # Single linkage algorithm
          cluster_dist[i,j] = min(dist_mat[clusters[[i]], clusters[[j]]])
        }else if(dist_method == "complete"){
          # Complete linkage algorithm
          cluster_dist[i,j] = max(dist_mat[clusters[[i]], clusters[[j]]])
        }else if(dist_method == "average"){
          # Average linkage algorithm
          cluster_dist[i,j] = sum(dist_mat[clusters[[i]], clusters[[j]]])/(length(clusters[[i]]) * length(clusters[[j]]))
        }else if(dist_method == "SSD"){
          # Sum of squares of deviations
          cluster_dist[i,j] = sum(dist_mat[clusters[[i]], clusters[[j]]]^2)/(length(clusters[[i]]) * length(clusters[[j]]))
        }else stop("Cluster distance method is not correct.")
      }
    }
  }

  dimnames(cluster_dist) = list(names, names)
  return(cluster_dist)
}









