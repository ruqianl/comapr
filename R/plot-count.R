
#' plotCount
#'
#' Plot the number of COs per sample group or per chromosome
#'
#' @param co_count
#' GRange or RangedSummarizedExperiment object, returned by \code{countCO}
#' @param group_by, the column name in `colData(co_count)` that specify the
#' grouping factor. Or the character vector contains the unique prefix of
#' sample names that are used for defining different sample groups. If missing
#' all samples are assumed to be from one group
#' @param by_chr, whether it should plot each chromosome separately
#' @param plot_type, determins what type the plot will be, choose from "error_bar"
#'  or "hist". Only relevant when by_chr=TRUE
#' @return ggplot object
#' @importFrom plyr mapvalues
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10,
#'                              minCellSNP = 0)
#' s1_count <- countCOs(s1_rse_state)
#' plotCount(s1_count)
#' @export
#'

setGeneric("plotCount",
           function(co_count,
                    by_chr=FALSE,
                    group_by  = "sampleGroup",
                    plot_type="error_bar")
             standardGeneric("plotCount"))


#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'RangedSummarizedExperiment',
                                    by_chr='missing',
                                    group_by='missing',
                                    plot_type='ANY'),
          function(co_count, by_chr = FALSE, group_by, plot_type){
            co_count$sampleGroup <- "all"
            plot_count(co_count = co_count,group_by  = "sampleGroup",
                       by_chr = by_chr,plot_type=plot_type)


          })

#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'RangedSummarizedExperiment',
                                    by_chr='missing',
                                    group_by='character',
                                    plot_type='ANY'),
          function(co_count, by_chr = FALSE, group_by, plot_type){


              plot_count(co_count = co_count,group_by  = group_by,
                         by_chr = FALSE,plot_type=plot_type)



          })

#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'RangedSummarizedExperiment',
                                    by_chr='logical',
                                    group_by='character',
                                    plot_type='ANY'),
          function(co_count, by_chr = FALSE , group_by, plot_type){
            plot_count(co_count = co_count,group_by  = group_by,
                       by_chr = by_chr,plot_type=plot_type)

          })


#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'RangedSummarizedExperiment',
                                    by_chr='logical',
                                    group_by='missing',
                                    plot_type='ANY'),
          function(co_count, by_chr = FALSE, group_by, plot_type){
            co_count$sampleGroup <- "all"
            plot_count(co_count = co_count,group_by  = "sampleGroup",
                       by_chr = by_chr,plot_type=plot_type)

          })




#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'GRanges',
                                    by_chr='logical',
                                    group_by='missing',
                                    plot_type='ANY'),
          function(co_count, by_chr = FALSE, group_by, plot_type){
           tmp_counts <- SummarizedExperiment(rowRanges = granges(co_count),
                                 colData = data.frame(sampleGroup=rep("all",ncol(mcols(co_count)))),
                                 assays=list(co_count=mcols(co_count)))
            plot_count(co_count = tmp_counts,group_by  = "sampleGroup",
                       by_chr = by_chr,plot_type=plot_type)

          })

#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'GRanges',
                                    by_chr='missing',
                                    group_by='missing',
                                    plot_type='ANY'),
          function(co_count, by_chr, group_by, plot_type){
            tmp_counts <-  SummarizedExperiment(rowRanges = granges(co_count),
                                                colData = data.frame(
                                                  sampleGroup=rep("all",
                                                                  ncol(mcols(co_count)))),
                                                assays=list(co_count=mcols(co_count)))

            plot_count(co_count = tmp_counts,group_by  = "sampleGroup",
                       by_chr = FALSE,plot_type=plot_type)

          })
#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'GRanges',
                                    by_chr='missing',
                                    group_by='character',
                                    plot_type='ANY'),
          function(co_count, by_chr, group_by, plot_type){
            sampleGroup <- rep("all",ncol(mcols(co_count)))
            for(group_prefix in group_by){
              sampleGroup[grep(group_prefix,
                               colnames(mcols(co_count))) ] <- group_prefix
            }

            tmp_counts <- SummarizedExperiment(rowRanges = granges(co_count),
                                                colData = data.frame(sampleGroup=sampleGroup),
                                                assays=list(co_count=mcols(co_count)))
            plot_count(co_count = tmp_counts,group_by  = "sampleGroup",
                       by_chr = FALSE,plot_type=plot_type)

          })
#'@rdname plotCount
setMethod("plotCount",signature = c(co_count = 'GRanges',
                                    by_chr='logical',
                                    group_by='character',
                                    plot_type='ANY'),
          function(co_count, by_chr, group_by, plot_type="error_bar"){
            if(length(plot_type)==0){plot_type <- "error_bar"}

            sampleGroup <- rep("all",ncol(mcols(co_count)))
            for(group_prefix in group_by){
              sampleGroup[grep(group_prefix,
                               colnames(mcols(co_count))) ] <- group_prefix
            }

            tmp_counts <- SummarizedExperiment(rowRanges = granges(co_count),
                                               colData = data.frame(sampleGroup=sampleGroup),
                                               assays=list(co_count=mcols(co_count)))
            plot_count(co_count = tmp_counts,group_by  = "sampleGroup",
                       by_chr = by_chr,plot_type=plot_type)

          })
#'@noRd
plot_count <- function(co_count, by_chr=FALSE, group_by  = "sampleGroup",
                      plot_type=c("error_bar","hist") ){

  stopifnot(group_by %in% colnames(colData(co_count)))
  col_to_plot <- unique(colData(co_count)[,group_by])
  sample_group_colors <- RColorBrewer::brewer.pal(ifelse(length(col_to_plot)> 2,
                                                         length(col_to_plot),3),
                                                  name = "Set1")

  names(sample_group_colors)[seq_along(col_to_plot)] <- col_to_plot

  chr <- BC <- COs <- sampleGroup <- NULL
  ChrCOs <- meanCOs <- sd <- lower <- upper <- NULL
  if(by_chr){

    tmp <- assay(co_count)
    tmp$chr <- GenomicRanges::seqnames(co_count)
    tmp <- data.frame(tmp,check.names = FALSE) %>%
      tidyr::pivot_longer(cols = colnames(co_count),
                          values_to = "COs", names_to = "BC")
    tmp$sampleGroup <- plyr::mapvalues(tmp$BC, from = colnames(co_count),
                                       to = colData(co_count)[,group_by])

    stopifnot(plot_type %in% c("error_bar","hist"))

    if(plot_type == "error_bar"){
      suppressMessages(
        p <- tmp %>% dplyr::group_by(chr,BC) %>%
          summarise(ChrCOs = sum(COs),
                    sampleGroup = unique(sampleGroup)) %>%
          dplyr::group_by(chr,sampleGroup) %>%
          dplyr::summarise(meanCOs = mean(ChrCOs),
                           lower = meanCOs-sd(ChrCOs)/sqrt(length(ChrCOs)),
                           upper = meanCOs+sd(ChrCOs)/sqrt(length(ChrCOs))) %>%
          ggplot(mapping = aes(x=chr,y = meanCOs,
                               color=chr))+geom_point()+
          geom_errorbar(mapping = aes(ymin = lower,
                                      ymax = upper))+
          facet_wrap(.~sampleGroup)+
          theme_classic(base_size = 18))
    } else {
      suppressMessages(
        p <- tmp %>% dplyr::group_by(chr,BC) %>%
          summarise(ChrCOs = sum(COs),
                    sampleGroup = unique(sampleGroup)) %>%
          dplyr::group_by(chr,sampleGroup) %>%
          dplyr::mutate(ChrCOs = as.character(ChrCOs)) %>%
          ggplot()+geom_bar(mapping = aes(x=ChrCOs),
                            stat ="count")+facet_wrap(.~chr))
    }
  } else {

    plt_df <- data.frame(COs = colSums(as.matrix(assay(co_count))),

                         sampleGroup = colData(co_count)[,group_by])
    p <- ggplot(data = plt_df,
                mapping = aes(x = sampleGroup, color = sampleGroup,
                              y = COs))+geom_boxplot()+
      geom_jitter(width = 0.2,
                  height = 0.0)+
      theme_classic(base_size = 16)+
      scale_color_manual(values = sample_group_colors)
  }
  p

}
