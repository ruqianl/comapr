#' getAFTracks
#'
#'
#' Generate the raw alternative allele frequencies tracks for all cells in the
#' columns of provided `co_count`
#' @inheritParams getCellAFTrack
#'
#' @author Ruqian Lyu
#' @export
#' @return a list object, in which each element is a list of two items with the
#' cell's alternative allele frequency DataTrack and the called crossover ranges.
#' @examples
#'   demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10,
#'                              minCellSNP = 0)
#' s1_counts <- countCOs(s1_rse_state)
#'
#' af_co_tracks <- getAFTracks(chrom ="chr1",
#'                                path_loc = demo_path,
#'                                sampleName = "s1",
#'                                barcodeFile = paste0(demo_path,
#'                                                     "s1_barcodes.txt"),
#'                                co_count = s1_counts)
#'
getAFTracks <-  function(chrom = "chr1",
                            path_loc = "./output/firstBatch/WC_522/",
                            sampleName = "WC_522",
                            nwindow = 80,
                            barcodeFile,
                            co_count,
                            snp_track = NULL){
  initial_barcodes <- read.table(file = barcodeFile)
  whichCells <- match(colnames(co_count), initial_barcodes$V1)

  dpMM <- Matrix::readMM(file = paste0(path_loc, sampleName,"_",
                                  chrom,"_totalCount.mtx"))

  dpMM <-dpMM[,whichCells]
  altMM <- Matrix::readMM(file = paste0(path_loc, sampleName,"_",
                                   chrom, "_altCount.mtx"))

  altMM <- altMM[,whichCells]
  af_data <- altMM/dpMM

  lapply(colnames(co_count), function(cell){
    ithCell <- match(cell,colnames(co_count))
    cell_af <- af_data[,ithCell]
    keep_snp <- !is.na(cell_af)
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
                                  name = paste0(cell, " AF"),
                                  data = cell_af[keep_snp],
                                  window = nwindow,
                                  na.rm=TRUE,
                                  aggregation = "mean",
                                  type="p")
      co_range_cell1 <- getCellCORange(co_count, cellBarcode = cell)
      co_range_cell1[seqnames(co_range_cell1)==chrom,]
      # ht <- HighlightTrack(trackList = af_track,
      #                      co_range_cell1[seqnames(co_range_cell1)=="chr1",],
      #                      chromosome = "chr1" )
      list(af_track=af_track, co_range =co_range_cell1)

  })

}
