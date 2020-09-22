

#' findDupSamples
#'
#' Find the duplicated samples by look at the number of matching genotypes
#' in all pair-wise samples
#'
#' @inherit countGT
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.text grid.rect gpar grid.circle
#' @importFrom circlize colorRamp2
#' @param threshold the frequency cut-off of number of matching genotypes out of
#' all geneotypes for determining whether the pair of samples
#' are duplicated, defaults to \code{0.99}. NAs are regarded as same genotypes for
#' two samples if they both have NA for a marker.
#'
#' @param plot whether a frequency plot should be generated for paired-samples
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
#' dups <- findDupSamples(cr_geno,plot = TRUE)
findDupSamples <- function(geno, threshold = 0.99, plot =TRUE,
                           in_text =TRUE){

  # similarity.matrix<-apply(gt_matrix,2,function(x)colSums(identical(x,gt_matrix[,1])))
  # diag(similarity.matrix)<-0
  #
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



  if(plot){
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    if(in_text){
       htmap <- Heatmap(out_freq,
            column_title  = "Samples' pair-wise frequencies of having
            same genotype across all markers",
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
            name = "Sample pair-wise\ncorrelation")

    } else {
      htmap <- Heatmap(out_freq,
                       column_title  = "Samples' pair-wise frequencies of having
            same genotype across all markers",
                       col = col_fun,cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       name = "Sample pair-wise\ncorrelation",
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


#' findDupMarkers
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
#'
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

#' getDistortedMarkers
#'
#' Marker disortation detection using chisq-test
#'
#' @details
#' We expect the genotypes to appear with the frequenceis of 1:1 homo_alt:hets.
#' We use chisq-test for finding markers that have genotypes among samples that
#' are significantly diffferent from the 1:1 ratio and report them
#'
#' @param p the expected geneotype ratio in a numeric vector, defauls to c(0.5,0.5)
#' @importFrom stats chisq.test p.adjust
#' @inheritParams countGT
#' @author Ruqian Lyu
#'
#' @export

getDistortedMarkers <- function(geno, p = c(0.5,0.5)){
 geno.table <- sapply(rownames(geno), function(marker){
   list(Het = ifelse(is.na(table(geno[marker,],useNA = "no")["Het"]),
                     0,
                     table(geno[marker,],useNA = "no")["Het"]),
        Homo_alt = ifelse(is.na(table(geno[marker,],useNA = "no")["Homo_alt"]),
                          0,
                          table(geno[marker,],useNA = "no")["Homo_alt"]))
 })
 geno.table <- data.frame(Markers = colnames(geno.table),
                          No.Hets =  as.character(unlist(geno.table["Het",])),
                          No.Homo_alt = as.character(unlist(geno.table["Homo_alt",])),
                          stringsAsFactors = FALSE)

 pvals <- sapply(as.character(geno.table$Markers), function(marker){
   ctest <- chisq.test(as.numeric(geno.table[geno.table$Markers==marker,2:3]),
                       p = p)
   ctest$p.value
 })

 #names(pvals) == geno.table$Markers
 geno.table$Pvals <- pvals
 geno.table$Adj.pvals <- p.adjust(pvals, method = "BH")

 return(geno.table[order(geno.table$Adj.pvals),])
}


#' plotGTFreq
#'
#' Function to plot the genotypes for all samples faceted by genotype

#' @importFrom plotly plot_ly subplot
#' @importFrom reshape2 melt
#' @inheritParams countGT
#' @author Ruqian Lyu
#' @param interactive, it determines whether an interactive
#' plot will be generated.
#' @param color_set, the RColorBrewer::brewer.pal color set names 
#' @export
#' @examples
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[1,] <- rep("Fail",dim(or_geno)[2])
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' ft_gt <- filterGT(cr_geno)
#' plotGTFreq(ft_gt)

plotGTFreq <- function(geno, interactive = FALSE,color_set="Set1"){

  geno.table <- sapply(colnames(geno), function(sample){
    list(Het = ifelse(is.na(table(geno[,sample],useNA = "no")["Het"]),
                      0,
                      table(geno[,sample],useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(geno[,sample],useNA = "no")["Homo_alt"])
                           ,0,
                           table(geno[,sample],useNA = "no")["Homo_alt"]))
  })
  geno.table <- data.frame(samples = colnames(geno.table),
                           Freq.Hets =  as.numeric(unlist(geno.table["Het",]))/nrow(geno),
                           Freq.Homo_alt = as.numeric(unlist(geno.table["Homo_alt",])/nrow(geno)),
                           stringsAsFactors = FALSE)
  pltdf <- melt(geno.table)

  if(interactive){
    ply1 <- plot_ly(pltdf, x=~.data$samples,y=~.data$value,type = "scatter",
                    color = ~.data$variable,mode = "markers",colors = color_set)
    return(ply1)
  } else {

    stplt1 <- ggplot(data = pltdf)+
      geom_point(mapping = aes_string(x = "samples", y = "value",
                               color = "variable"))+
      theme_classic()+
      ylab("Genotype Frequecies for each sample")+labs(color ="Genotype")+
      theme(axis.text.x = element_text(angle =-90))
    return(stplt1)

  }

}




