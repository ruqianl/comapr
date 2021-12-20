#' getAFTracks
#'
#'
#' Generate the raw alternative allele frequencies tracks for all cells in the
#' columns of provided `co_count`
#'
#' @inheritParams getCellAFTrack
#' @importFrom Matrix readMM
#' @importFrom Gviz DataTrack
#' @importFrom GenomicRanges start
#' @author Ruqian Lyu
#' @export
#' @return a list object, in which each element is a list of two items with the
#' cell's alternative allele frequency DataTrack and the called crossover ranges.
#' @examples
#' demo_path <- system.file("extdata",package = "comapr")
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
#'                                barcodeFile = file.path(demo_path,
#'                                                     "s1_barcodes.txt"),
#'                                co_count = s1_counts)
#'
getAFTracks <-  function(chrom = "chr1",
                            path_loc = "./output/firstBatch/WC_522/",
                            sampleName = "WC_522",
                            nwindow = 80,
                            barcodeFile,
                            co_count,
                            chunk = 1000L,
                            snp_track = NULL){
  stopifnot(file.exists(barcodeFile))
  initial_barcodes <- read.table(file = barcodeFile)
  whichCells <- match(colnames(co_count), initial_barcodes$V1)

  dpMM <- readMM(file = file.path(path_loc, paste0(sampleName,"_",
                                  chrom,"_totalCount.mtx")))

  dpMM <- dpMM[,whichCells]
  altMM <- readMM(file = file.path(path_loc, paste0(sampleName,"_",
                                   chrom, "_altCount.mtx")))

  altMM <- altMM[,whichCells]

  af_data <- altMM/dpMM

  lapply(colnames(co_count), function(cell){
    ithCell <- match(cell,colnames(co_count))
    cell_af <- af_data[,ithCell]
    keep_snp <- !is.na(cell_af)
    if(is.null(snp_track)){
      snp_anno <- read.table(file=file.path(path_loc, paste0(sampleName,"_",
                                         chrom, "_snpAnnot.txt")),
                             header=TRUE)
      snp_pos <- snp_anno$POS

    } else {
      stopifnot(unlist(dimnames(seqinfo(snp_track))) == chrom)
      snp_pos <- start(snp_track)
    }

      af_track <- DataTrack(GRanges(seqnames = chrom,
                                          IRanges(start=snp_pos[keep_snp],
                                                  width = 1,
                                                  genome = "mm10")),
                                  name = paste0(cell, " AF"),
                                  data = cell_af[keep_snp],
                                  window = nwindow,
                                  na.rm = TRUE,
                                  aggregation = "mean",
                                  type = "p")
      co_range_cell1 <- getCellCORange(co_count, cellBarcode = cell)
      co_range_cell1[seqnames(co_range_cell1) == chrom,]

      list(af_track = af_track, co_range = co_range_cell1)

  })

}
