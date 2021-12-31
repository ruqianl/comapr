#' plotGeneticDist
#'
#' Plotting the calculated genetic distanced for each bin or marker interval
#' supplied by the GRanges object
#'
#' @importFrom dplyr filter left_join
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom scales unit_format
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb genome seqinfo
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width
#' @importFrom GenomicRanges tileGenome binnedAverage mcols
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels
#' @importFrom S4Vectors Rle
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @param gr GRanges object with genetic distances calculated for marker
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
#' @examples
#' data(coCount)
#' dist_se <- calGeneticDist(coCount)
#' plotGeneticDist(SummarizedExperiment::rowRanges(dist_se))

plotGeneticDist <- function(gr,bin=TRUE,chr=NULL,cumulative=FALSE){
  col_to_plot <- colnames(mcols(gr))
  sample_group_colors <- brewer.pal(ifelse(length(col_to_plot)>2,
                                                         length(col_to_plot),3),
                                                  name="Set1")
  names(sample_group_colors)[seq_along(col_to_plot)] <- col_to_plot

  sample_group_colors <- sample_group_colors[col_to_plot]

  if(cumulative){
    mcols(gr) <- apply(mcols(gr),2 ,
                       function(x,
                                seq = as.character(seqnames(gr))){
      temp_df <- data.frame(x=x,seq=seq) %>% group_by(seq) %>%
        mutate(cum = cumsum(x))
      temp_df$cum
      })
  }
  plot_df <- data.frame(gr)
  ## in case colnames for groups (| -> .) was escaped
  rename_cols <- (ncol(plot_df)-length(col_to_plot)+1):ncol(plot_df)
  colnames(plot_df)[rename_cols] <- col_to_plot

  plot_df <- plot_df %>% mutate(x_tick = 0.5*(.data$start + .data$end))

  plot_df <- plot_df %>% pivot_longer(cols=col_to_plot,
                                             names_to = "SampleGroup",
                                             values_to = "bin_dist")
  x_tick <- bin_dist <- end <- SampleGroup <- NULL
    if(is.null(chr)){
      p <-  plot_df %>%
        ggplot()+
        geom_step(mapping = aes(x = x_tick ,y = bin_dist,
                              color=SampleGroup),size =1)
  } else {
      p <- plot_df %>% filter(seqnames %in% chr) %>%
        ggplot()+
        geom_step(mapping = aes(x = end,y=bin_dist,
                                color=SampleGroup),size =1)
  }

    p <- p + scale_x_continuous(labels = unit_format(unit = "M",
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
}

#' Plot cumulative genetic distances across the genome
#'
#' This function takes the calculated genetic distances for all
#' marker intervals across all chromosomes provided and plot the
#' cumulative genetic distances
#' @importFrom GenomicRanges mcols
#' @importFrom dplyr select mutate
#' @export
#' @param  gr, GRanges object with genetic distances calculated for marker
#' intervals
#' @return A ggplot object
#' @examples
#' data(coCount)
#' dist_se <- calGeneticDist(coCount)
#' plotWholeGenome(SummarizedExperiment::rowRanges(dist_se))
plotWholeGenome <- function(gr){

  end <- chr_len <- tot <- all_of <- x_tick <- bin_dist <- SampleGroup <- NULL
  BPcum <- NULL

  mcols(gr) <- apply(mcols(gr), 2, function(x) {
    temp_df <- data.frame(x = x) %>%
      mutate(cum = cumsum(x))
    temp_df$cum
  })

  gr_df <- data.frame(gr,check.names = FALSE)
  col_to_plot <- gsub("\\|",".", colnames(mcols(gr)))
  don <- gr_df %>%

    # Compute chromosome size
    group_by(seqnames) %>%
    summarise(chr_len=max(end)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(gr_df, ., by=c("seqnames"="seqnames")) %>%

    # Add a cumulative position of each SNP
    mutate( BPcum=end+tot,
                   x_tick = (0.5 * (start + end) + tot))

    plot_df <- don %>% pivot_longer(cols = all_of(col_to_plot),
                                         names_to = "SampleGroup",
                                         values_to = "bin_dist")
    p <- plot_df %>% ggplot() + geom_step(mapping = aes(x = x_tick,
                                                        y = bin_dist,
                                                        color = SampleGroup),
                                          size = 1)+

      theme_classic(base_size = 18) +
      xlab("Chromosome positions \n cumulative whole genome") +
      ylab("cumulative centiMorgans")

    axisdf <-  don %>% group_by(seqnames) %>%
      summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

    p <- p+scale_x_continuous(labels = axisdf$seqnames, breaks= axisdf$center)+
      theme_classic(base_size = 15) +
      xlab("Chromosome positions \n cumulative whole genome") +
      ylab("cumulative centiMorgans")
    p
}
