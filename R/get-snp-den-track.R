#' getSNPDensityTrack
#'
#' Generate the SNP density DataTrack (from `Gviz`) for selected chromosome
#' @inheritParams getCellAFTrack
#' @param log, whether the histogram of SNP density should be plotted on log
#' scale (log10)
#' @param plot_type, the DataTrack plot type, default to be `hist`
#'
#' @author Ruqian Lyu
#' @export
#' @return DataTrack object plotting the SNP density histogram
#'
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' snp_track <- getSNPDensityTrack(chrom ="chr1",
#'                                path_loc = demo_path,
#'                                sampleName = "s1")
getSNPDensityTrack <- function(chrom = "chr1",
                               sampleName = "s1",
                               path_loc = ".",
                               nwindow = 80,
                               plot_type = "hist",
                               log = TRUE){

  snp_anno <- read.table(file=paste0(path_loc, sampleName,"_",chrom, "_snpAnnot.txt"),
                         header=TRUE)
  aggregation_fun <- ifelse(log, function(x) { log10(sum(x)+1)}, sum)
  snp_track <- DataTrack( GRanges(seqnames = chrom,
                                  IRanges(start=snp_anno$POS,
                                          width = 1,
                                          genome = "mm10")),
                          name = "SNPs denstiy",
                          data = rep(1,nrow(snp_anno)),
                          window = nwindow,
                          aggregation  = aggregation_fun,
                          type=plot_type)

  snp_track
}
