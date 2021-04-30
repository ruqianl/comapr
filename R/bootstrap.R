#' bootstrapDist
#'
#' Generating distribution of sample genetic distances
#'
#' It takes the crossover counts for samples in multiple groups that is returned
#' by `countCO`. It then draws samples from a group with replacement and
#' calculate the distribution of relevant statistics.
#'
#' @param co_gr GRanges or RangedSummarizedExperiment object that contains the
#' crossover counts for each marker interval across all samples.
#' Returned by \code{countCOs}
#' @param mapping_fun character default to "k" (kosambi mapping function). It
#' can be one of the mapping functions: "k","h"
#' @param B integer the number of sampling times
#' @param group_by, the prefix for each group that we need to generate
#' distributions for(only when co_gr is a GRanges object). Or the column name
#' for colData(co_gr) that contains the
#' group factor (only when co_gr is a RangedSummarizedExperiment object)

#' @importFrom BiocParallel bplapply
#' @export
#' @return lists of numeric genetic distances for multiple samples
#' @examples
#' bootsDiff <- bootstrapDist(coCount, group_by = "sampleGroup")
#' @author Ruqian Lyu
#'
bootstrapDist <- function(co_gr, B = 1000, mapping_fun = "k", group_by)
{
    if(is(co_gr,"RangedSummarizedExperiment")) {
        count_matrix <- as.matrix(assay(co_gr))

        group_idx <- lapply(unique(as.character(colData(co_gr)[, group_by])),
                        function(gp) {
                          grep(gp, colData(co_gr)[, group_by])})
        } else {
          count_matrix <- as.matrix(mcols(co_gr))
          group_idx <- lapply(group_by, function(gp){
            grep(gp,colnames(mcols(co_gr)))
    })
  }

  #group_size <- sapply(group_idx, length)
  boots_result <- BiocParallel::bplapply(seq_len(B), function(b) {
    boots_sample1 <- sample(unlist(group_idx[1]),replace = TRUE)
    boots_sample2 <- sample(unlist(group_idx[2]),replace = TRUE)

    rb_rate1 <- rowMeans(count_matrix[,boots_sample1])
    dist_1 <- switch(mapping_fun,
                     k = 25 * log((1 + 2*rb_rate1)/(1 - 2*rb_rate1)),
                     h = -50 * log(1 - 2*rb_rate1))
    rb_rate2 <- rowMeans(count_matrix[, boots_sample2])
    dist_2 <- switch(mapping_fun,
                     k = 25 * log((1 + 2*rb_rate2)/(1 - 2*rb_rate2)),
                     h = -50 * log(1 - 2*rb_rate2))
    sum(dist_1) - sum(dist_2)
  })
  unlist(boots_result)
}


#' permuteDist
#'
#' Permutation test of two sample groups
#'
#' It shuffles the group labels for the samples and calculate a difference
#' between two groups after shuffling.
#'
#' @inheritParams bootstrapDist
#' @importFrom BiocParallel bplapply

#' @export
#' @return A list of three elements. `permutes` of length B with numeric
#' differences of permuted group differences,`observed_diff` the observed
#' genetic distances of two groups, `nSample`, the number of samples in the
#' first and second group.
#' @examples
#' perms <- permuteDist(coCount, group_by = "sampleGroup")
#' @author Ruqian Lyu
#'
permuteDist <- function(co_gr,B=100,mapping_fun="k",group_by){

  if(is(co_gr,"RangedSummarizedExperiment")){
    ## group_by contains the
    count_matrix <- as.matrix(assay(co_gr))
    group_idx <- lapply(unique(as.character(colData(co_gr)[, group_by])),
                        function(gp) {
                          grep(gp, colData(co_gr)[, group_by])})

  } else {
    ## co_gr is of GRanges class
    count_matrix <- as.matrix(mcols(co_gr))
    #stopifnot(length(group_by)==2)

    group_idx <- lapply(group_by, function(gp){
    grep(gp,colnames(mcols(co_gr)))
    })
  }
  stopifnot(length(group_idx)==2)
  two_g_idx <- unlist(group_idx)
  len_1 <- length(group_idx[[1]])
  ## there is a potential issue here.
  ## when the randomly generated "samples" have co_rate >0.5.
  ## Currently perhaps just abandon these trials
  perm_result <- BiocParallel::bplapply(seq_len(B),function(b){
      g1_idx <- sample(two_g_idx,size = len_1,replace = FALSE)
      g2_idx <- setdiff(two_g_idx,g1_idx)

      rb_rate_1 <- rowMeans(count_matrix[,g1_idx])
      rb_rate_2 <- rowMeans(count_matrix[,g2_idx])
      dist_1 <- switch(mapping_fun,
                     k = 25*log((1+2*rb_rate_1)/(1-2*rb_rate_1)),
                     h = -50 * log(1 - 2 * rb_rate_1))
      dist_2 <- switch(mapping_fun,
                       k = 25*log((1+2*rb_rate_2)/(1-2*rb_rate_2)),
                       h = -50 * log(1 - 2 * rb_rate_2))
      sum(dist_1)-sum(dist_2)
      #rs

  })
  ## calculate observed diff
  g1_idx <- group_idx[[1]]
  g2_idx <- group_idx[[2]]

  rb_rate_1 <- rowMeans(count_matrix[,g1_idx])
  rb_rate_2 <- rowMeans(count_matrix[,g2_idx])
  dist_1 <- switch(mapping_fun,
                   k = 25*log((1+2*rb_rate_1)/(1-2*rb_rate_1)),
                   h = -50 * log(1 - 2 * rb_rate_1))
  dist_2 <- switch(mapping_fun,
                   k = 25*log((1+2*rb_rate_2)/(1-2*rb_rate_2)),
                   h = -50 * log(1 - 2 * rb_rate_2))
  observed_diff <- sum(dist_1)-sum(dist_2)


  list(permutes=unlist(perm_result),
       observed_diff = observed_diff,
       nSample=c(length(g1_idx),length(g2_idx)))
}


