
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
  n_samples <- rowSums(!is.na(geno))
  n_markers <- colSums(!is.na(geno))
  
  pl_df <- data.frame(counts=c(n_samples,n_markers),
                      type = c(rep("by_marker", each = length(n_samples)),
                               rep("by_sample", each= length(n_markers))))
                      # name = c(names(n_samples),
                      #          names(n_markers)))
  type_colors <- c("by_marker"="#003f5c",
                   "by_sample"="#ffa600")
  if(plot){
    if(interactive){
      p <- plot_ly(data =pl_df,
                   x = ~type,y = ~counts, hoverinfo="text", jitter = 0.3,
                   boxpoints = "outliers",color = ~type,
                   text = ~paste('</br> counts ID: ',counts),
                      type = c("box"))

    #
      # name = c("by_sample"="No. samples by marker",
      #          "by_marker"="No. markers by sample")
      return(list(ply = p,n_samples = rowSums(!is.na(geno)),
                  n_markers = colSums(!is.na(geno))))
      
    } else {
      
      p <- ggplot(data = pl_df,mapping = aes(y = counts,x = type))+
        geom_boxplot(size=1.5,aes(fill=type),alpha = 0.5)+
        geom_jitter(position = "jitter",aes(color=type),size=2,alpha=.8 )+
        theme_classic(base_size=22)+facet_wrap(.~type,scales = "free")+
        ylab("Number of samples/markers")+
#        ggtitle("Genotype counts\nfor markers/samples")+
        scale_color_manual(values=type_colors)+
        scale_fill_manual(values=type_colors)
    
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
#' The filtered genotype matrix
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
setGeneric("filterGT")

#' @rdname filterGT
setMethod("filterGT",signature = c(geno ="matrix",min_markers = "numeric",
                                   min_samples = "numeric"),
          function(geno ,min_markers,
                   min_samples){
            filterGT(geno ,min_markers,
                     min_samples)
  
})

#' @rdname filterGT
setMethod("filterGT",signature = c(geno ="GRanges",min_markers = "numeric",
                                   min_samples = "numeric"),
          function(geno ,min_markers,
                   min_samples){
            
            gt_counts <- countGT(mcols(geno),plot = FALSE)
            keep_markers <- gt_counts$n_samples >= min_samples
            keep_samples <- gt_counts$n_markers >= min_markers
            geno <- geno[keep_markers,]
            
            mcols(geno) <- mcols(geno)[,keep_samples]
            message(paste0( "filter out ",sum(keep_markers==FALSE)," marker(s)"))
            message(paste0( "filter out ",sum(keep_samples==FALSE)," sample(s)"))
            
            return(geno)
            
          })

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
#' @noRd
#' @return
#' a ggplot2 object for markers with missing genotype across samples
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
