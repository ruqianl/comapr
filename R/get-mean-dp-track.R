#' getMeanDPTrack
#'
#' Generate the mean DP (Depth) DataTrack (from Gviz) for cells
#' @inheritParams getCellAFTrack
#' @importFrom Matrix rowMeans
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


  dpMM <- readMM(file =paste0(path_loc,sampleName,"_",
                                      chrom,"_totalCount.mtx"))
  if(!is.null(selectedBarcodes))
  dpMM <- dpMM[,match(selectedBarcodes, initial_barcodes$V1)]

  meanDP <- rowMeans(dpMM)

  snp_pos <- .get_snp_pos(snp_track = snp_track, path_loc = path_loc,
               sampleName = sampleName, chrom = chrom)
  aggregation_fun <- ifelse(log, function(x) { log10(sum(x)+1)}, sum)

  meanDP_track <- DataTrack(GRanges(seqnames = chrom,
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



#' getCellDPTrack
#' Generates the DataTrack for plotting DP of a selected cell
#'
#' It plots the total allele counts for the selected cell.
#' @inheritParams getCellAFTrack
#' @inheritParams getSNPDensityTrack
#' @param chrom, the chromosome
#' @param path_loc, the path prefix to the output files from sscocaller
#' including "*_totalCount.mtx"
#'
#' @param nwindow, the number of windows for binning the chromosome
#' @param barcodeFile, the barcode file containing the list of cell barcodes
#'  used as the input file for sscocaller
#' @param cellBarcode, the selected cell barcode
#' @param snp_track, the SNP position track which is used for obtaining the SNP
#' chromosome locations. It could be omitted and the SNP positions will be
#' acquired from the "*_snpAnnot.txt" file.
#' @param chunk, A integer scalar indicating the chunk size to use, i.e., number
#' of rows to read at any one time.

#' @return The DataTrack object defined in \code{\link[Gviz]{DataTrack}}
#' @importFrom  Gviz DataTrack
#' @author Ruqian Lyu
#' @export
#' @examples
#' demo_path <- system.file("extdata",package = "comapr")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10,
#'                              minCellSNP = 0)
#' s1_counts <- countCOs(s1_rse_state)
#' dp_co_tracks <- getCellDPTrack(chrom ="chr1",
#'                                path_loc = demo_path,
#'                                sampleName = "s1",
#'                                barcodeFile = file.path(demo_path,
#'                                                     "s1_barcodes.txt"),
#'                                cellBarcode = "BC1")
#'
getCellDPTrack <-  function(chrom = "chr1",
                            path_loc = "./output/firstBatch/WC_522/",
                            sampleName = "WC_522",
                            nwindow = 80,
                            barcodeFile,
                            cellBarcode,
                            snp_track = NULL,
                            chunk = 1000L,
                            log =TRUE,
                            plot_type = "hist"){
  stopifnot(length(cellBarcode)==1)
  stopifnot(file.exists(barcodeFile))
  initial_barcodes <- read.table(file = barcodeFile)
  whichCell <- match(cellBarcode, initial_barcodes$V1)

  dpMM <- readColMM(file = file.path(path_loc,paste0(sampleName,"_",
                                  chrom,"_totalCount.mtx")),
                    which.col = whichCell,
                    chunk = chunk)

  dpMM <- dpMM[,whichCell]

  snp_pos <- .get_snp_pos(snp_track = snp_track, path_loc = path_loc,
                          sampleName = sampleName, chrom = chrom)
  aggregation_fun <- ifelse(log, function(x) { log10(mean(x)+1)}, "mean")

  dp_track <- DataTrack(GRanges(seqnames = chrom,
                                      IRanges(start=snp_pos,
                                              width = 1,
                                              genome = "mm10")),
                              name = paste0(cellBarcode, " AF"),
                              data = dpMM,
                              window = nwindow,
                              na.rm=TRUE,
                              aggregation = aggregation_fun,
                              type=plot_type)

  dp_track
}

#' get_snp_pos from track or snp_annot_file
#' @return a list of SNP positions for selected chromosome
#' @keywords internal
#' @noRd
.get_snp_pos <- function(snp_track = NULL,path_loc,sampleName,chrom){

  if(is.null(snp_track)){
    snp_anno <- read.table(file=file.path(path_loc,paste0(sampleName,"_",
                                       chrom, "_snpAnnot.txt")),
                           header=TRUE)
    snp_pos <- snp_anno$POS

  } else {
    stopifnot(unlist(dimnames(seqinfo(snp_track))) == chrom)
    snp_pos <- snp_track@range@ranges@start
  }
  snp_pos
}
