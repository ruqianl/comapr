#' getCellCORange
#'
#' It finds the crossover intervals for a selected cell
#'
#' @param co_count, `GRanges` or `RangedSummarizedExperiment` object,
#' @param cellBarcode, the selected cell's barcode
#'
#' @return GRange object containing the crossover intervals for the selected
#' cell
#' @importFrom SummarizedExperiment assay assay<- rowRanges
#' @importFrom GenomicRanges mcols mcols<- reduce
#' @author Ruqian Lyu
#' @export
#' @examples
#'   demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10,
#'                              minCellSNP = 0)
#' s1_counts <- countCOs(s1_rse_state)
#'
#' co_ranges <- getCellCORange(cellBarcode = "BC1",
#'                             co_count = s1_counts)
getCellCORange <- function(co_count, cellBarcode){
  stopifnot(class(co_count) %in% c("GRanges","RangedSummarizedExperiment"))

  if(is(co_count,"GRanges")){
    cell1_co <- mcols(co_count)[,cellBarcode]
    cell1_coRange <- granges(co_count)[cell1_co!=0]
    cell1_co <- cell1_co[cell1_co!=0]
  } else {
    cell1_co <- assay(co_count)[,cellBarcode]
    cell1_coRange <- SummarizedExperiment::rowRanges(co_count)[cell1_co!=0]
    cell1_co <- cell1_co[cell1_co!=0]

    mcols(cell1_coRange) <- cell1_co

  }

  co_range_cell1 <- GenomicRanges::reduce(cell1_coRange,min.gapwidth=2)
  co_range_cell1
}

#' getCellAFTrack
#' Generates the DataTracks for plotting AF and crossover regions
#'
#' It plots the raw alternative allele frequencies and highlight the crossover
#' regions for the selected cell.
#' @param chrom, the chromosome
#' @param path_loc, the path prefix to the output files from sscocaller including
#' "*_totalCount.mtx" and "_altCount.mtx"
#' @param sampleName, the sample name, which is the prefix of sscocaller's output
#' files
#'
#' @param nwindow, the number of windows for binning the chromosome
#' @param barcodeFile, the barcode file containing the list of cell barcodes used
#' as the input file for sscocaller
#' @param cellBarcode, the selected cell barcode
#' @param snp_track, the SNP position track which is used for obtaining the SNP
#' chromosome locations. It could be omitted and the SNP positions will be acquired
#' from the "*_snpAnnot.txt" file.
#' @param co_count, `GRange` or `RangedSummarizedExperiment` object,
#' returned by \code{countCO} that contains the crossover intervals and the number
#' of crossovers in each cell.
#' @param chunk, An integer scalar indicating the chunk size to use,
#' i.e., number of rows to read at any one time.
#'
#' @return The DataTrack object defined in \code{\link[Gviz]{DataTrack}}
#' @importFrom  Gviz DataTrack
#' @author Ruqian Lyu
#' @export
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10,
#'                              minCellSNP = 0)
#' s1_counts <- countCOs(s1_rse_state)
#' af_co_tracks <- getCellAFTrack(chrom ="chr1",
#'                                path_loc = demo_path,
#'                                sampleName = "s1",
#'                                barcodeFile = paste0(demo_path,
#'                                                     "s1_barcodes.txt"),
#'                                cellBarcode = "BC1",
#'                                co_count = s1_counts)
#'
getCellAFTrack <-  function(chrom = "chr1",
                            path_loc = "./output/firstBatch/WC_522/",
                            sampleName = "WC_522",
                            nwindow = 80,
                            barcodeFile,
                            cellBarcode,
                            co_count,
                            snp_track = NULL,
                            chunk = 1000L){
  stopifnot(length(cellBarcode)==1)
  initial_barcodes <- read.table(file = barcodeFile)
  whichCell <- match(cellBarcode, initial_barcodes$V1)

  dpMM <- readColMM(file = paste0(path_loc, sampleName,"_",
                                  chrom,"_totalCount.mtx"),which.col = whichCell,
                    chunk = chunk)

  dpMM <- dpMM[,whichCell]

  altMM <- readColMM(file = paste0(path_loc, sampleName,"_",
                                   chrom, "_altCount.mtx"),which.col = whichCell,
                     chunk = chunk)


  altMM <- altMM[,whichCell]
  af_data <- altMM/dpMM
  keep_snp <- !is.na(af_data)

  if(is.null(snp_track)){
     snp_anno <- read.table(file=paste0(path_loc,sampleName,"_",
                                        chrom, "_snpAnnot.txt"),
                            header=TRUE)
     snp_pos <- snp_anno$POS

  } else {
    stopifnot(unlist(dimnames(seqinfo(snp_track))) == chrom)
    snp_pos <- snp_track@range@ranges@start
  }

  af_track <- Gviz::DataTrack(GRanges(seqnames = chrom,
                                IRanges(start=snp_pos[keep_snp],
                                        width = 1,
                                        genome = "mm10")),
                       name = paste0(cellBarcode, " AF"),
                       data = af_data[keep_snp],
                       window = nwindow,
                       na.rm=TRUE,
                       aggregation = "mean",
                       type="p",
                       ylim=0:1)
  co_range_cell1 <- getCellCORange(co_count, cellBarcode = cellBarcode)
  co_range_cell1[seqnames(co_range_cell1)==chrom,]
  # ht <- HighlightTrack(trackList = af_track,
  #                      co_range_cell1[seqnames(co_range_cell1)=="chr1",],
  #                      chromosome = "chr1" )
 list(af_track=af_track, co_range =co_range_cell1)
}

