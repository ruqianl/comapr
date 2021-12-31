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
                         snp_track = NULL){
  whichCells <- .get_cell_idx(barcodeFile = barcodeFile, 
                              cellBarcode = colnames(co_count))
  
  dpMM <- .get_cells_mm(path_loc = path_loc,sampleName = sampleName, 
                        chrom = chrom, whichCell = whichCells, type = "total")
  altMM <- .get_cells_mm(path_loc = path_loc,sampleName = sampleName, 
                         chrom = chrom, whichCell = whichCells, type = "alt")
  
  af_data <- altMM/dpMM

  lapply(colnames(co_count), function(cell){
    ithCell <- match(cell,colnames(co_count))
    cell_af <- af_data[,ithCell]
    keep_snp <- !is.na(cell_af)
    snp_pos  <- .get_snp_pos(snp_track = snp_track, path_loc = path_loc,
                 sampleName = sampleName, chrom = chrom)
    
    .get_af_track_and_co(chrom = chrom, snp_pos = snp_pos, 
                           keep_snp = keep_snp, cellBarcode = cell, 
                           af_data = cell_af,
                           nwindow = nwindow ,co_count = co_count)

  })

}
