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
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[,1] <- or_geno[,5]
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' dups <- findDupSamples(cr_geno)
findDupSamples <- function(geno, threshold = 0.99, 
                           in_text = FALSE){
  
  # similarity.matrix<-apply(gt_matrix,2,function(x)colSums(identical(x,gt_matrix[,1])))
  # diag(similarity.matrix)<-0
  #
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

  # if(plot){
  #   col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  #   if(in_text){
  #     htmap <- Heatmap(out_freq,
  #                      column_title  = "Samples' pair-wise frequencies of having
  #           same genotype across all markers",
  #                      col = col_fun,
  #                      cell_fun = function(j, i, x, y, width, height, fill) {
  #                        grid.rect(x = x, y = y, width = width, height = height,
  #                                  gp = gpar(col = "grey", fill = NA))
  #                        if(i == j) {
  #                          grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
  #                                    gp = gpar(fontsize = 10))
  #                        } else if(i > j) {
  #                          grid.rect(x = x, y = y,  
  #                                    width = width, 
  #                                    height = height,
  #                                    gp = gpar(col = "grey",
  #                                              fill = col_fun(out_freq[i, j])))
  #                        } else {
  #                          grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
  #                                    gp = gpar(fontsize = 10))
  #                        }
  #                      }, cluster_rows = FALSE, cluster_columns = FALSE,
  #                      name = "Sample pair-wise\ncorrelation")
  #     
  #   } else {
  #     htmap <- Heatmap(out_freq,
  #                      column_title  = "Samples' pair-wise frequencies of having
  #           same genotype across all markers",
  #                      col = col_fun,cluster_rows = FALSE,
  #                      cluster_columns = FALSE,
  #                      name = "Sample pair-wise\ncorrelation",
  #                      row_names_gp = gpar(fontsize = 10),
  #                      column_names_gp = gpar(fontsize = 10))
  #     
  #   }
  #   
  # }
  
  l_out_freq <- out_freq
  l_out_freq[lower.tri(out_freq,diag=TRUE)] <- 0
  
  
  dups_index <- which(l_out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(l_out_freq)[x])
  #if(plot) dups <- list(htmap,dups)
  return(dups)
}
