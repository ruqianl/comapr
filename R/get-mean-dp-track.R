#' getMeanDPTrack
#'
#' Generate the mean DP (Depth) DataTrack (from Gviz) for cells
#' @inheritParams getCellAFTrack
#' @param log, whether the histogram of SNP density should be plotted on log
#' scale (log10)
#' @param plot_type, the DataTrack plot type, default to be `hist`
#' @param selectedBarcodes, the selected cell barcodes which should be the
#' barcodes that have been called crossovers for. If not supplied then all cells
#' are counted.
#'
#' @author Ruqian Lyu
#' @export
#' @return DataTrack object plotting the mean DP histogram for windowed
#' chromosomes
#'
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' meanDP_track <- getMeanDPTrack(chrom ="chr1",
#'                                path_loc = demo_path,
#'                                sampleName = "s1",
#'                                barcodeFile = paste0(demo_path,"s1_barcodes.txt"))
#'
getMeanDPTrack <- function( chrom = "chr1",
                            path_loc,
                            nwindow = 80,
                            sampleName,
                            barcodeFile,
                            plot_type = "hist",
                            selectedBarcodes = NULL,
                            snp_track = NULL,
                            log =TRUE){
  initial_barcodes <- read.table(file = barcodeFile)


  dpMM <- Matrix::readMM(file =paste0(path_loc,sampleName,"_",
                                      chrom,"_totalCount.mtx"))
  if(!is.null(selectedBarcodes))
  dpMM <-dpMM[,match(selectedBarcodes, initial_barcodes$V1)]

  meanDP <- Matrix::rowMeans(dpMM)

  if(is.null(snp_track )){
    snp_anno <- read.table(file=paste0(path_loc, sampleName,"_",chrom, "_snpAnnot.txt"),
                           header=TRUE)
    snp_pos <- snp_anno$POS

  } else {
    stopifnot(unlist(dimnames(seqinfo(snp_track))) == chrom)
    snp_pos <- snp_track@range@ranges@start
  }
  aggregation_fun <- ifelse(log, function(x) { log10(sum(x)+1)}, sum)

  meanDP_track <- DataTrack( GRanges(seqnames = chrom,
                                     IRanges(start= snp_pos,
                                             width = 1,
                                             genome = "mm10")),
                             name = "meanDP across cells",
                             data = meanDP,window =100,
                             aggregation = aggregation_fun,
                             type=plot_type,
                             background.panel = "#FFFEDB",
                             background.title = "lightblue")
  meanDP_track
}
