
#' countGT
#' 
#' count how many samples have genotypes calls across markers and count how 
#' many markers that each individual has called genotypes for. This function 
#' helps identify poor samples or poor markers for filtering. It can also generate
#' plots that help identify outlier samples/markers
#'
#' @importFrom plotly plot_ly subplot
#' @importFrom gridExtra grid.arrange

#' @author Ruqian Lyu
#'
#' @param geno the genotype data.frame of markers by samples from output of
#' function \code{correctGT}
#'
#' @param plot, it determines whether a plot will be generated, defaults to TRUE
#' @param interactive, it determines whether an interactive plot will be generated
#' @export
#'
#' @return
#' A list of two elements including \code{n_markers} and \code{n_samples}

countGT <- function(geno, plot =TRUE,interactive=FALSE){
  if(plot){
    if(interactive){
      ply1 <- plot_ly(data = data.frame(marker_index =
                                          seq(1:length(rowSums(!is.na(geno)))),
                                        No.Samples =  rowSums(!is.na(geno)),
                                        marker_ID = rownames(geno)),
                      x = ~marker_index, y = ~No.Samples, hoverinfo="text",
                      text = ~paste('</br> Marker ID: ',marker_ID),
                      name = "No. samples by marker",
                      mode = "markers",
                      type = "scatter")
      
      ply2 <- plotly::plot_ly(data = 
                                data.frame(sample_index = 
                                             seq(1:length(colSums(!is.na(geno)))),
                                           No.Markers =  colSums(!is.na(geno)),
                                           sample_ID = colnames(geno)),hoverinfo="text",
                              x = ~sample_index, y = ~No.Markers,
                              text = ~paste('</br> Sample ID: ',sample_ID), 
                              name = "No. markers by sample",type = "scatter",
                              mode = "markers")
      
      p <- subplot(ply1,ply2)
      return(list(ply = p,n_samples = rowSums(!is.na(geno)),
                  n_markers = colSums(!is.na(geno))))
      
    } else {
      
      p1 <- ggplot()+
        geom_point(mapping = aes(x = seq(1:length(rowSums(!is.na(geno)))),
                                 y = rowSums(!is.na(geno))))+
        theme_classic()+
        ylab("Number of samples")+xlab("markers index")+
        ggtitle("No. samples by marker")
      
      p2 <- ggplot()+
        geom_point(mapping = aes(x = seq(1:length(colSums(!is.na(geno)))),
                                 y = colSums(!is.na(geno))))+
        theme_classic()+
        ylab("Number of markers")+
        xlab("samples index")+
        ggtitle("No. markers by sample")
      p <- grid.arrange(p1, p2, nrow = 1)
    }
    return(list(plot = p,
                n_samples = rowSums(!is.na(geno)),
                n_markers = colSums(!is.na(geno))))
  }
  
  return(list(n_samples = rowSums(!is.na(geno)),
              n_markers = colSums(!is.na(geno))))
}

#' filterGT
#'
#' Filter markers or samples that have too many missing values
#'
#' @inheritParams countGT
#'
#' @param min_markers the minimum number of markers for a sample to be kept
#' @param min_samples the minimum number of samples for a marker to be kept
#'
#' @details
#' This function takes the \code{geno} data.frame and filter the data.frame by
#' the provided cut-offs.
#'
#' @author Ruqian Lyu
#'
#' @return
#' A filtered genotype matrix
#' @export
#'
filterGT <- function(geno, min_markers = 5, min_samples = 3){
  
  gt_counts <- countGT(geno,plot = FALSE)
  keep_markers <- gt_counts$n_samples >= min_samples
  keep_samples <- gt_counts$n_markers >= min_markers
  
  message(paste0( "filter out ",sum(keep_markers==FALSE)," marker(s)"))
  message(paste0( "filter out ",sum(keep_samples==FALSE)," sample(s)"))
  
  return(geno[keep_markers, keep_samples])
  
}

#' Plot markers with missing genotypes, deprecate
#'
#' @author Ruqian Lyu
#'
#' @param geno the genotype data.frame of markers by samples
#' @param plot_wg
#' whether to plot all markers across whole genome or just the markers that are
#' ever missing across all samples
#' @param missing
#' the label in the matrix that is used for encoding the missing or failed data
#' @param plot_type
#' whether a 'dot' plot or a 'bar' plot should be drawn
#'
#' @details
#' This functions plots the missing markers in a 'dot' plot or 'bar' plot.
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#'
#' @keywords internal
#' @return
#' a ggplot2 object for markers with missing genotype across samples
#' @export
#'
#' @examples
#' or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' plotMissingGT(cr_geno)


plotMissingGT <- function(geno, missing = "Fail", plot_wg = FALSE,
                          plot_type = "dot"){
  
  stopifnot(plot_type == "dot" | plot_type == "bar")
  
  
  mis_matrix <- apply(geno, 2, function(es){
    is.na(es) | es == missing
  })
  
  plot_df <- melt(mis_matrix)
  
  if(plot_wg){
    switch (plot_type,
            dot = ggplot(data = plot_df)+
              geom_point(mapping =
                           aes_string(x = "Var1", colour = "value", y = "Var2"))+
              xlab("markers")+ylab("samples")+
              labs(colour = "is_missing")+
              scale_color_manual(values = c("TRUE" = "red",
                                            "FALSE"= "lightgrey"))+
              theme_classic()+theme(axis.text.x = element_text(angle = -90)),
            
            bar =  ggplot(data = plot_df)+
              geom_bar(mapping = aes_string( fill = "value", x = "Var1"))+
              xlab("markers")+ylab("samples")+
              labs(fill = "is_missing")+
              scale_fill_manual(values = c("TRUE" = "red",
                                           "FALSE"= "lightgrey"))+
              theme_classic()+theme(axis.text.x = element_text(angle = -90))
            
    )
    
    
  } else {
    remain_m <- plot_df %>% group_by(.data$Var1) %>%
      summarise(no_missing = sum(.data$value)) %>% filter(.data$no_missing >0)
    
    plot_df <- plot_df[plot_df$Var1 %in% remain_m$Var1,]
    switch (plot_type,
            dot = ggplot(data = plot_df)+
              geom_point(mapping =  aes_string(x = "Var1", colour = "value", 
                                               y = "Var2"))+
              xlab("markers")+ylab("samples")+
              labs(colour = "is_missing")+
              scale_color_manual(values = c("TRUE" = "red",
                                            "FALSE"= "lightgrey"))+
              theme_classic()+theme(axis.text.x = element_text(angle = -90)),
            
            bar =  ggplot(data = plot_df)+
              geom_bar(mapping = aes_string( fill = "value", x = "Var1"))+
              xlab("markers")+ylab("samples")+
              labs(fill = "is_missing")+
              scale_fill_manual(values = c("TRUE" = "red",
                                           "FALSE"= "lightgrey"))+
              theme_classic()+theme(axis.text.x = element_text(angle = -90))
            
    )
    
    
  }
}
