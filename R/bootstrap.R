#' bootstrapDist
#' 
#' Generating distribution of sample genetic distances 
#' 
#' It takes the crossover counts for samples in multiple groups that is returned
#' by `countCO`. It then draws samples from a group with replacement and 
#' calculate the distribution of relevant statistics.
#' 
#' @param co_gr GRanges the crossover counts for each marker interval across all
#' samples. Returned by \code{countCO}
#' @param mapping_fun character default to "k" (kosambi mapping function). It 
#' can be one of the mapping functions: "k","h"
#' @param B integer the number of sampling times
#' @param chr character the specific chromosome to test or NULL indicates all
#' chrs in the GRanges object will be summed together.
#' @param group_prefix the prefix for each group that we need to generate 
#' distributions for
#' @param BPPARAM, the bpparam backend registered for parallel computing

#' @importFrom BiocParallel bplapply
#' @export
#' @return lists of numeric genetic distances for multiple samples
#'  
#' @author Ruqian Lyu
#'
bootstrapDist <- function(co_gr,B=100,mapping_fun="k",chr=NULL,group_prefix,
                          BPPARAM=BiocParallel::bpparam()){
  if(is.null(chr)){
    chr <- unique(as.character(seqnames(co_gr)))
  }
  co_gr <- co_gr[as.character(seqnames(co_gr)) %in% chr,]
  group_idx <- sapply(group_prefix, function(gp){
    grep(gp,colnames(mcols(co_gr)))
  },USE.NAMES=T)
  ## there is a potential issue here.
  ## when the randomly generated "samples" have co_rate >0.5. 
  ## Currently perhaps just abandon these trials
  boots_result <- BiocParallel::bplapply(1:B,function(b){
    boots_dist <- sapply(group_idx, function(t){
      rs <- sample(t,size = length(t),replace = T)
      rb_rate <- rowMeans(as.matrix(mcols(co_gr))[,rs])
      dist <- switch(mapping_fun,
        k = 25*log((1+2*rb_rate)/(1-2*rb_rate)),
        h = -50 * log(1 - 2 * rb_rate))
      sum(dist)
      #rs
    })
    
  })
  do.call(rbind,boots_result)
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

#' @param group_prefix a character vector of two
#' @export
#' @return the list of length B with numeric differences of two groups
#' @author Ruqian Lyu
#' 
permuteDist <- function(co_gr,B=100,mapping_fun="k",chr=NULL,group_prefix,
                          BPPARAM=BiocParallel::bpparam()){
  stopifnot(length(group_prefix)==2)
  if(is.null(chr)){
    chr <- unique(as.character(seqnames(co_gr)))
  }
  co_gr <- co_gr[as.character(seqnames(co_gr)) %in% chr,]
  
  group_idx <- sapply(group_prefix, function(gp){
    grep(gp,colnames(mcols(co_gr)))
  },USE.NAMES=T)
  two_g_idx <- unlist(group_idx)
  len_1 <- length(group_idx[[1]])
  ## there is a potential issue here.
  ## when the randomly generated "samples" have co_rate >0.5. 
  ## Currently perhaps just abandon these trials
  perm_result <- BiocParallel::bplapply(1:B,function(b){
      g1_idx <- sample(two_g_idx,size = len_1,replace = F)
      g2_idx <- setdiff(two_g_idx,g1_idx)
      
      rb_rate_1 <- rowMeans(as.matrix(mcols(co_gr))[,g1_idx])
      rb_rate_2 <- rowMeans(as.matrix(mcols(co_gr))[,g2_idx])
      dist_1 <- switch(mapping_fun,
                     k = 25*log((1+2*rb_rate_1)/(1-2*rb_rate_1)),
                     h = -50 * log(1 - 2 * rb_rate_1))
      dist_2 <- switch(mapping_fun,
                       k = 25*log((1+2*rb_rate_2)/(1-2*rb_rate_2)),
                       h = -50 * log(1 - 2 * rb_rate_2))
      sum(dist_1)-sum(dist_2)
      #rs
  
    
  })
  unlist(perm_result)
} 


