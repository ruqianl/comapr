#' findDupSamples
#'
#' Find the duplicated samples by look at the number of matching genotypes
#' in all pair-wise samples
#'
#' @inherit countGT
#' @importFrom grid grid.text grid.rect gpar grid.circle
#' @importFrom circlize colorRamp2
#' @param threshold the frequency cut-off of number of matching genotypes out of
#' all geneotypes for determining whether the pair of samples
#' are duplicated, defaults to \code{0.99}. NAs are regarded as same genotypes for
#' two samples if they both have NA for a marker.

#' @param in_text whether text of frequencies should be displayed in the
#' heatmap cells
#' @export
#' @author Ruqian Lyu
#' @return
#' The paris of duplicated samples.
#'
#' @examples
#' data(snp_geno)
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[,1] <- or_geno[,5]
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.)
#' dups <- findDupSamples(cr_geno)
findDupSamples <- function(geno, threshold = 0.99,
                           in_text = FALSE){


  geno <- as.matrix(geno)
  n <- seq_len( ncol(geno) )

  #  Make combinations
  id <- expand.grid( n , n )

  #  Get result
  out <- matrix( colSums( (is.na(geno[ , id[,1] ]) & is.na(geno[ , id[,2] ])) |
                            (geno[ , id[,1] ] == geno[ , id[,2] ]),
                          na.rm =TRUE) , ncol = length(n) )
  dimnames(out)[1] <- dimnames(geno)[2]
  dimnames(out)[2] <- dimnames(geno)[2]

  out_freq <- out / dim(geno)[1]
  l_out_freq <- out_freq
  l_out_freq[lower.tri(out_freq,diag=TRUE)] <- 0


  dups_index <- which(l_out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(l_out_freq)[x])
  #if(plot) dups <- list(htmap,dups)
  return(dups)
}
