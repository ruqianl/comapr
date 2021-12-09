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
#' for `colData(co_gr)` that contains the
#' group factor (only when co_gr is a RangedSummarizedExperiment object)

#' @importFrom BiocParallel bplapply
#' @export
#' @return lists of numeric genetic distances for multiple samples
#' @examples
#' data(coCount)
#'
#' bootsDiff <- bootstrapDist(coCount, group_by = "sampleGroup",B=10)
#' @author Ruqian Lyu
#'
bootstrapDist <- function(co_gr, B = 1000, mapping_fun = "k", group_by )
{
  .check_mapping_fun(mapping_fun)
  count_m_idex <- .get_count_matrix(co_gr_rse = co_gr, group_by = group_by )
  group_idx <- count_m_idex$gi
  count_matrix <- count_m_idex$c

  #group_size <- sapply(group_idx, length)
  result_fun <- function(){

    bpl_fun <- function(b) {
      boots_sample1 <- sample(unlist(group_idx[1]),replace = TRUE)
      boots_sample2 <- sample(unlist(group_idx[2]),replace = TRUE)

      rb_rate1 <- rowMeans(count_matrix[,boots_sample1])
      dist_1 <- .rb_to_dist(rb_rate1,mapping_fun=mapping_fun)

      rb_rate2 <- rowMeans(count_matrix[, boots_sample2])
      dist_2 <- .rb_to_dist(rb_rate2,mapping_fun=mapping_fun)

      sum(dist_1) - sum(dist_2)
    }
    bplapply(seq_len(B), function(b) bpl_fun(b) )
  }
  boots_result <- result_fun()
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
#' data(coCount)
#' perms <- permuteDist(coCount, group_by = "sampleGroup",B=10)
#' @author Ruqian Lyu
#'
permuteDist <- function(co_gr,B=100,mapping_fun="k",group_by){
  .check_mapping_fun(mapping_fun)
  count_m_idex <- .get_count_matrix(co_gr_rse = co_gr, group_by = group_by )
  group_idx <- count_m_idex$gi
  count_matrix <- count_m_idex$c

  stopifnot(length(group_idx)==2)
  two_g_idx <- unlist(group_idx)
  len_1 <- length(group_idx[[1]])
  ## there is a potential issue here.
  ## when the randomly generated "samples" have co_rate >0.5.
  ## Currently perhaps just abandon these trials
  result_fun <- function(){
    bpl_fun <- function(b){
      g1_idx <- sample(two_g_idx,size = len_1,replace = FALSE)
      g2_idx <- setdiff(two_g_idx,g1_idx)

      rb_rate_1 <- rowMeans(count_matrix[,g1_idx])
      rb_rate_2 <- rowMeans(count_matrix[,g2_idx])

      dist_1 <- .rb_to_dist(rb_rate_1,mapping_fun=mapping_fun)
      dist_2 <- .rb_to_dist(rb_rate_2,mapping_fun=mapping_fun)

      sum(dist_1)-sum(dist_2)
      #rs
    }
    bplapply(seq_len(B),function(b) bpl_fun(b))
  }
  perm_result <- result_fun()
  ## calculate observed diff
  g1_idx <- group_idx[[1]]
  g2_idx <- group_idx[[2]]

  rb_rate_1 <- rowMeans(count_matrix[,g1_idx])
  rb_rate_2 <- rowMeans(count_matrix[,g2_idx])

  dist_1 <- .rb_to_dist(rb_rate_1,mapping_fun=mapping_fun)
  dist_2 <- .rb_to_dist(rb_rate_2,mapping_fun=mapping_fun)

  observed_diff <- sum(dist_1)-sum(dist_2)


  list(permutes=unlist(perm_result),
       observed_diff = observed_diff,
       nSample=c(length(g1_idx),length(g2_idx)))
}

#' Kosambi mapping function
#' @noRd
#' @keywords internal
.k_map <- function(rb_rate){
  25*log((1+2*rb_rate)/(1-2*rb_rate))
}
#' Haldane mapping function
#' @noRd
#' @keywords internal
.h_map <- function(rb_rate){
  -50 * log(1 - 2 * rb_rate)
}

.rb_to_dist <- function(rb_rate,mapping_fun){
  switch(mapping_fun,
         k = .k_map(rb_rate),
         h = .h_map(rb_rate))
}
#' Get the count matrix from RSE or Granges object's mcols
#' @noRd
#' @keywords internal
.get_count_matrix <- function(co_gr_rse,group_by){

  if(is(co_gr_rse,"RangedSummarizedExperiment")){
    ## group_by contains the
    count_matrix <- as.matrix(assay(co_gr_rse))
    group_idx <- lapply(unique(as.character(colData(co_gr_rse)[, group_by])),
                        function(gp) {
                          grep(gp, colData(co_gr_rse)[, group_by])})

  } else {
    ## co_gr is of GRanges class
    count_matrix <- as.matrix(mcols(co_gr_rse))
    #stopifnot(length(group_by)==2)
    group_idx <- lapply(group_by, function(gp){
      grep(gp,colnames(mcols(co_gr_rse)))
    })
  }
  list(c = count_matrix, gi = group_idx)
}
