#' plotGeneticDist
#' 
#' Plotting the calculated genetic distanced for each bin or marker interval 
#' supplied by the GRanges object
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom scales unit_format
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb fetchExtendedChromInfoFromUCSC genome seqinfo  
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width 
#' @importFrom GenomicRanges tileGenome binnedAverage mcols
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels
#' @importFrom S4Vectors Rle
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @param gr GRanges object with genetic distances calculated for listed
#' intervals
#' @param bin TRUE or FALSE, indicating whether the supplied GRange objecct 
#' is for binned interval
#' @param chr the specific chrs selected to plot
#' @param cumulative TRUE or FALSE, indicating whether it plots the bin-wise
#' genetic distances or the cumulative distances

#' @return 
#' ggplot2 plot
#' @author Ruqian Lyu
#' 
#' 
plotGeneticDist <- function(gr,bin=TRUE,chr=NULL,cumulative=FALSE){
  col_to_plot <- colnames(GenomicRanges::mcols(gr))
  sample_group_colors <- RColorBrewer::brewer.pal(ifelse(length(col_to_plot)>2,
                                                         length(col_to_plot),3),
                                                  name="Set1")
  names(sample_group_colors)[seq_along(col_to_plot)] <- col_to_plot
  
  if(cumulative){
   
    GenomicRanges::mcols(gr) <- apply(mcols(gr),2 ,function(x,seq = as.character(seqnames(gr))) {
      temp_df <- data.frame(x=x,seq=seq) %>% dplyr::group_by(seq) %>% 
        dplyr::mutate(cum = cumsum(x))
      temp_df$cum
      })
  } 
  plot_df <- data.frame(gr)
  ## in case colnames for groups (| -> .) was escaped 
  colnames(plot_df)[(ncol(plot_df)-length(col_to_plot)+1):ncol(plot_df)] <- col_to_plot
  
  plot_df <- plot_df %>% dplyr::mutate(x_tick = 0.5*(.data$start + .data$end))
#                                       bin_dist = mcols(gr)[,1]) 
  
  plot_df <- plot_df %>% tidyr::pivot_longer(cols=col_to_plot,
                                             names_to = "SampleGroup",
                                             values_to = "bin_dist") 
  x_tick <- bin_dist <- end <- SampleGroup <- NULL
    if(is.null(chr)){
      p <-  plot_df %>%
        ggplot()+
        geom_step(mapping = aes(x = x_tick ,y = bin_dist,
                              color=SampleGroup),size =1)
  } else {
      p <- plot_df %>% dplyr::filter(seqnames %in% chr) %>% 
        ggplot()+
        geom_step(mapping = aes(x = end,y=bin_dist,
                                color=SampleGroup),size =1)
  }
    
    p <- p + scale_x_continuous(labels = scales::unit_format(unit = "M",
                                                        scale = 1e-6))+
      facet_wrap(.~seqnames,ncol=1,scales = "free")+
      theme_classic(base_size=18)+
      xlab("Chromosome positions")+
      scale_color_manual(values = sample_group_colors)
    
    if(cumulative){
      p + ylab("cumulative centiMorgans")
    } else {
      p + ylab("centiMorgans")
    }
    
  # } else {
  #   ## plot cumulative distances for chr 
  #   if(is.null(chr)){
  #     p <-  plot_df %>% dplyr::group_by(seqnames,SampleGroup)%>% 
  #       dplyr::mutate(cum_dist = cumsum(bin_dist)) %>%
  #       ggplot()+
  #       geom_step(mapping = aes(x = x_tick ,y=cum_dist,
  #                               color=SampleGroup),size =1)
  #   } else {
  #     p <- plot_df %>% dplyr::group_by(seqnames) %>% 
  #       dplyr::filter(seqnames %in% chr )%>%
  #       dplyr::mutate(cum_dist = cumsum(bin_dist)) %>%
  #       ggplot()+
  #       geom_step(mapping = aes(x = end,y=cum_dist,
  #                               color=SampleGroup),size =1)
  #   }
  #   


}
