
#' findDupMarkers, deprecate
#'
#' Find the duplicated markers by look at the number of matching genotypes across 
#' samples
#' 
#' @inheritParams countGT
#' @param threshold the frequency cut-off for determining whether the pair of 
#' markers are duplicated, defaults to \code{0.9}
#' @param in_text whether text of frequencies should be displayed in the
#' heatmap cells
#' @param plot whether a frequency heatmap plot should be generated. If the 
#' number of markers is large, do not plot.
#' @return
#' The paris of duplicated markers.
#'
#' @export
#' @author Ruqian Lyu
#' @keywords internal
#' @examples
#'   or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[,1] <- or_geno[,5]
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                      alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' dups <- findDupMarkers(cr_geno,plot = TRUE)


findDupMarkers <- function(geno, threshold = 0.99, plot =FALSE,
                           in_text = TRUE){
  
  n <- seq_len( nrow(geno) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  out <- matrix( rowSums( geno[  id[,1], ] == geno[ id[,2], ],
                          na.rm =TRUE) , nrow = length(n) )
  dimnames(out)[1] <- dimnames(geno)[1]
  dimnames(out)[2] <- dimnames(geno)[1]
  
  out_freq <- out / dim(geno)[2]
  
  if (plot) {
    
    #  heatmap(out_freq,scale = NULL,margins = c(7, 7),xlab ="Markers Pair-wise frequencies of having\nsame genotype across all samples")
    
    col_fun <- colorRamp2(c(0, 0.5, 0.8,1), c("blue", "white", "red","darkred"))
    
    if (in_text) {
      htmap <- Heatmap(out_freq,
                       column_title  = "Frequencies of pair-wise markers for having
                     same genotype across all samples",
                       col = col_fun,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.rect(x = x, y = y, width = width, height = height,
                                   gp = gpar(col = "grey", fill = NA))
                         if(i == j) {
                           grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                                     gp = gpar(fontsize = 10))
                         } else if(i > j) {
                           grid.rect(x = x, y = y,  width = width, height = height,
                                     gp = gpar(col = "grey",
                                               fill = col_fun(out_freq[i, j])))
                         } else {
                           grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                                     gp = gpar(fontsize = 10))
                         }
                       }, cluster_rows = FALSE, cluster_columns = FALSE,
                       name = "Marker pair-wise\ncorrelation",
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))
    } else {
      htmap <- Heatmap(out_freq,
                       column_title  = "Frequencies of pair-wise markers for having
                     same genotype across all samples",
                       col = col_fun,cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       name = "Marker pair-wise\ncorrelation",
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))
      
    }
    
  }
  
  l_out_freq <- out_freq
  l_out_freq[lower.tri(out_freq,diag=TRUE)] <- 0
  
  
  dups_index <- which(l_out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(l_out_freq)[x])
  if(plot) dups <- list(htmap,dups)
  return(dups)
  
}