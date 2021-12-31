
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
 
  af_data <- .get_af_data(chrom = chrom, path_loc = path_loc, 
                          sampleName = sampleName, barcodeFile = barcodeFile,
                          cellBarcode = cellBarcode, co_count = co_count,
                          chunk = chunk)
  keep_snp <- !is.na(af_data)

  snp_pos <- .get_snp_pos(snp_track = snp_track, path_loc = path_loc,
                          sampleName = sampleName, chrom = chrom)

  .get_af_track_and_co(chrom = chrom, snp_pos = snp_pos, 
                       keep_snp = keep_snp, cellBarcode = cellBarcode, 
                       af_data = af_data,
                       nwindow = nwindow ,co_count = co_count)
}
#' Generate DataTrack and CO range in a list 
#' @noRd
#' @keywords internal
.get_af_track_and_co <- function(chrom, snp_pos, keep_snp, cellBarcode, af_data,
                                 nwindow,co_count){
  
  af_track <- DataTrack(GRanges(seqnames = chrom,
                                IRanges(start=snp_pos[keep_snp],
                                        width = 1)),
                        name = paste0(cellBarcode, " AF"),
                        data = af_data[keep_snp],
                        window = nwindow,
                        na.rm=TRUE,
                        aggregation = "mean",
                        type="p",
                        ylim=0:1)
  co_range_cell1 <- getCellCORange(co_count, cellBarcode = cellBarcode)
  co_range_cell1[seqnames(co_range_cell1)==chrom,]
  
  list(af_track = af_track, co_range = co_range_cell1)
}
#' 
#' Get cell(s) AF data
#' @noRd
#' @keywords internal
.get_af_data <- function( chrom, path_loc, sampleName, barcodeFile,
                          cellBarcode, co_count, chunk){
  
  whichCell <- .get_cell_idx(barcodeFile = barcodeFile,
                             cellBarcode = cellBarcode )

  dpMM <- .get_cells_mm(path_loc = path_loc, sampleName = sampleName,
                           chrom = chrom, whichCell = whichCell,
                           chunk = chunk,
                           type = "total") 
  altMM <- .get_cells_mm(path_loc = path_loc, sampleName = sampleName,
                             chrom = chrom, whichCell = whichCell,
                             chunk = chunk,
                             type = "alt") 
  af_data <- altMM/dpMM
  af_data
}

#' Find cell(s) index
#' @noRd
#' @keywords internal
#' 
.get_cell_idx <- function(barcodeFile,
                          cellBarcode){
  stopifnot(file.exists(barcodeFile))
  initial_barcodes <- read.table(file = barcodeFile)
  whichCell <- match(cellBarcode, initial_barcodes$V1)
  whichCell
}


#' Get cells total DP or alt counts
#' @noRd
#' @keywords internal
.get_cells_mm <- function(path_loc, sampleName,
                          chrom, whichCell,chunk = 1000L,
                          type = "total"){
  stopifnot(type %in% c("total","alt"))
  
  if(type == "total"){
    if(length(whichCell)>1){
      mm <- readMM(file = file.path(path_loc, paste0(sampleName, "_",
                                    chrom, "_totalCount.mtx")))
    } else {
      mm <- readColMM(file = file.path(path_loc, paste0(sampleName, "_",
                                  chrom, "_totalCount.mtx")),
                    which.col = whichCell,
                    chunk = chunk)
      }
    
  } else {
    if(length(whichCell)>1){
      mm <- readMM(file = file.path(path_loc, paste0(sampleName,"_",
                                    chrom,"_altCount.mtx")))
    }
    else{
      mm <- readColMM(file = file.path(path_loc, paste0(sampleName,"_",
                                  chrom,"_altCount.mtx")),
                    which.col = whichCell,
                    chunk = chunk)
    }

  }
  mm[,whichCell]
}
