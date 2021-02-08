#Author: Somnath Tagore, Ph.D.
#Title: Functions: Master Regulator Analysis of COVID data
#Script Name: Functions_MasterRegulators.R
#Last Updated: 02/05/2021

#' Louvain clustering over a range of resolution parameters (uses getComMembership from the MUDAN package)
#'
#' @param dat.mat Matrix of data (features X samples)
#' @param dist.mat Distance matrix
#' @param rmin Lowest resolution to try. Default of 10.
#' @param rmax Maximum resolution to try. Default of 100.
#' @param rstep Step size for resolution parameter. Default of 10.
#' @param verbose Switch to control terminal read out of progress. Default of TRUE.
#' @return A list of two lists; 'clusterings', which contains the cluster labels and 'sils', which has silhouette scores
#' @export
LouvainClustering <- function(dat.mat, dist.mat, rmin = 10, rmax = 100, rstep = 10, verbose = TRUE) {
  # iterate through resolution params
  res <- rmin
  clusterings <- list()
  sils <- list()
  while (res <= rmax) {
    if (verbose) {print(paste('Clustering with res=', res, '...', sep = ''))}
    # cluster
    set.seed(343)
    lclust <- MUDAN::getComMembership(t(dat.mat), k = res, method = igraph::cluster_walktrap, verbose = FALSE)
    clusterings[[paste('res', res, sep = '')]] <- lclust
    # evaluate
    sil.score <- cluster::silhouette(as.integer(lclust), dist.mat)
    sils[[paste('res', res, sep = '')]] <- mean(sil.score[,3])
    # iterate resolution param
    res <- res + rstep
  }
  return(list('clusterings' = clusterings, 'sils' = sils))
}

#' Stouffer integrates the given vector of data.
#' 
#' @param dat.vect Vector of data to be integrated
#' @param weight.vect Vector of weights, if specified.
#' @return Stouffer integrated value for the given data.
#' @export
StoufferIntegrate <- function(dat.vect, weight.vect) {
  if (!missing(weight.vect)) {
    s.int <- sum(dat.vect * weight.vect) / sqrt(sum(weight.vect**2))
  } else {
    s.int <- sum(dat.vect) / length(dat.vect)
  }
  return(s.int)
}

#' Master regulators by Stouffer integration.
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param clust.vect Vector of cluster labels.
#' @param weight.vect Vector of weights, if specified.
#' @return Integrated 
StoufferMRs <- function(dat.mat, clust.vect) {
  
}

#' Identifies cluster specific master regulators using the Mann-Whitney U-test.
#' Approximates p-vals using a normal distribution for n > 30.
#' 
#' @param dat.mat Matrix of data.
#' @param clust.vect Vector of cluster labels.
#' @return For each cluster, two sorted lists of log p-values, split by positive / negative
#' @export
MasterRegulators <- function(dat.mat, clust.vect) {
  # identify cluster names
  clust.names <- unique(clust.vect)
  clust.mrs <- list()
  # cluster specific mrs
  for (cn in clust.names) {
    print(paste('Identifying MRs for cluster', cn))
    # set up labels and test statistics
    clust.samps <- which(clust.vect == cn); n.1 <- length(clust.samps)
    ref.samps <- which(clust.vect != cn); n.2 <- length(ref.samps)
    u.mean <- (n.1 * n.2) / 2; u.sd <- sqrt((n.1 * n.2 * (n.1 + n.2 + 1)) / 12)
    if (n.1 < 30 | n.2 < 30) { print('WARNING: Group size <30, normal approximation may not be appropriate...') }
    # generate tests; scale; transform to p-value
    clust.wStat <- apply(dat.mat, 1, function(x) {wilcox.test(x[clust.samps], x[ref.samps])$statistic} )
    clust.zScore <- sapply(clust.wStat, function(x) {(x - u.mean) / u.sd} )
    clust.logp <- sapply(clust.zScore, function(x) {log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE)})
    # check medians
    median.dif <- apply(dat.mat, 1, function(x) {sign(median(x[clust.samps]) - median(x[ref.samps]))} ) 
    # sort and return
    mr.lists <- list('positive' = sort(clust.logp[which(median.dif == 1)]),
                     'negative' = sort(clust.logp[which(median.dif == -1)]))
    clust.mrs[[cn]] <- mr.lists
  }
  return(clust.mrs)
}

#' Identifies a vector of colors for a given number of clusters. Internal function.
#'
#' @param k Number of clusters.
#' @param offset Optional argument to shift colors along color wheel.
#' @return A vector of hues
#' @export
ClusterColors <- function(k, offset = 0) {
  hues <- seq(15, 375, length = k + 1) + offset
  return(hcl(h = hues, l = 65, c = 100)[1:k])
}

#' Generates breaks for a color scale based on quantiles.
#'
#' @param dat.mat Data matrix (features X samples).
#' @param n Number of breaks to generate. If not specified, uses first three stdevs.
#' @return Numeric vector of break values.
#' @export
QuantileBreaks <- function(xs, n) {
  if (!missing(n)) {
    breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  } else {
    breaks <- quantile(xs, c(0.003, 0.05, 0.32, 0.5, 0.68, 0.95, 0.997))
  }
  return(unique(breaks))
}

# Save the pheatmap object
save_pheatmap_pdf <- function(x, filename, width=15, height=15) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
