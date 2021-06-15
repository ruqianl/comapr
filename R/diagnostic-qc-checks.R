#' perCellChrQC
#'
#' A function that parses output ('_viSegInfo.txt' )from `sscocaller`
#' \url{https://gitlab.svi.edu.au/biocellgen-public/sscocaller} and
#' generate cell cell (per chr) summary statistics
#' @param path, the path to the files, with name patterns *{chrom}_vi.mtx,
#' *{chrom}_viSegInfo.txt, end with slash
#' @param barcodeFile, defaults to NULL, it is assumed to be in the same d
#' irectory as the other files and with name sampleName_barcodes.txt
#' @param chroms, the character vectors of chromosomes to parse. Multiple
#' chromosomes' results will be concated together.
#' @param sampleName, the name of the sample to parse which is used as prefix
#' for finding relevant files for the underlying sample
#' @param doPlot, whether a plot should returned, default to TRUE
#'
#' @return a list object that contains the data.frame with summarised statistics
#' per chr per cell and a plot (if doPlot)
#'
#'
#' @examples
#' demo_path <-system.file("extdata",package = "comapr")
#' pcQC <- perCellChrQC(sampleName="s1",chroms=c("chr1"),
#' path=paste0(demo_path,"/"),
#' barcodeFile=NULL)
#' @export
#' @author Ruqian Lyu
perCellChrQC <- function(sampleName,chroms = c("chr1","chr7","chr15"),
                         path,barcodeFile = NULL,doPlot = TRUE){
  if (is.null(barcodeFile)) {
    barcodeFile <- paste0(path, sampleName, "_barcodes.txt")
  }
  stopifnot(file.exists(barcodeFile))

  barcodes <- read.table(file = barcodeFile, stringsAsFactors = FALSE,
                         col.names = "barcodes")

  segInfo_list <- bplapply(chroms, function(chr) {

    segInfo <- read.table(file = paste0(path, sampleName,
                                        "_", chr, "_viSegInfo.txt"),
                          stringsAsFactors = FALSE,
                          col.names = c("ithSperm", "Seg_start", "Seg_end",
                                        "logllRatio", "nSNP", "State"))
    segInfo$Chrom <- chr
    segInfo

  })
  segInfo_chrs <- do.call(rbind,segInfo_list)
  rm(segInfo_list)

  if(grepl("chr",segInfo_chrs$Chrom[1])){
    chr_levels <- paste0("chr",seq_len(23))
  } else {
    chr_levels <- paste0(seq_len(23))
  }

  plt_df <- segInfo_chrs %>% dplyr::group_by(.data$ithSperm,.data$Chrom) %>%
    dplyr::summarise(totalSNP = sum(.data$nSNP),
                     nCORaw = length(.data$nSNP)-1) %>%
    dplyr::mutate(Chrom = factor(.data$Chrom,levels =  chr_levels))
  plt_df$barcode <- barcodes$barcodes[as.numeric(gsub("ithSperm",
                                                      "",plt_df$ithSperm))+1]
  plt_df$ithSperm <- NULL
  if(doPlot){
    plt <-
      ggplot(plt_df)+geom_point(mapping = aes(x = .data$nCORaw,
                                              y = .data$totalSNP))+
      scale_y_log10()+scale_x_log10()+facet_wrap(.~Chrom,ncol=5)+theme_classic()
    list(cellQC = plt_df,
         plot = plt)
  } else {
    plt_df
  }
}

#' countBinState
#'
#' Bins the chromosome into supplied number of bins and find the state of
#' the chromosome bins across all sperms
#'
#' This function is used for checking whether chromosome segregation pattern
#' obeys the expected ratio.
#'
#' @param chr, character, the chromosome to check
#' @param snpAnno, data.frame, the SNP annotation for the supplied chromosome
#' @param viState, dgTMatrix/Matrix, the viterbi state matrix, output from
#' `sscocaller`
#' @param genomeRange, GRanges object with seqlengths information for the genome
#' @param ntile, integer, how many tiles the chromosome is binned into
#'
#' @return a data.frame that contains chromosome bin segregation ratio
#'
#'
#' @examples
#' \donttest{
#' chrom_info <- GenomeInfoDb::getChromInfoFromUCSC("mm10")
#' seq_length <- chrom_info$size
#' names(seq_length) <- chrom_info$chrom
#'
#' dna_mm10_gr <- GenomicRanges::GRanges(
#'   seqnames = Rle(names(seq_length)),
#'   ranges = IRanges(1, end = seq_length, names = names(seq_length)),
#'   seqlengths = seq_length)
#'
#' GenomeInfoDb::genome(dna_mm10_gr) <- "mm10"
#' demo_path <- system.file("extdata",package = "comapr")
#' sampleName <- "s1"
#' chr <- "chr1"
#' vi_mtx <- Matrix::readMM(file = paste0(demo_path,"/", sampleName, "_",
#'                                        chr, "_vi.mtx"))
#'
#' snpAnno <- read.table(file = paste0(demo_path,"/", sampleName,
#'                                     "_", chr, "_snpAnnot.txt"),
#'                                     stringsAsFactors = FALSE,
#'                       header = TRUE)
#'
#' countBinState(chr = "chr1",snpAnno = snpAnno,
#' viState = vi_mtx,genomeRange = dna_mm10_gr,
#' ntile = 1)
#' }
#' @author Ruqian Lyu
#' @export
countBinState <- function(chr,
                          snpAnno,
                          viState,
                          genomeRange,
                          ntile = 5){

  tiles <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(genomeRange)[chr],
                                     ntile = ntile)
  binned_dna_mm10_gr <- unlist(tiles)

  vi_gr <- GenomicRanges::GRanges(seqnames = chr,
                                  ranges = IRanges::IRanges(start = snpAnno$POS,
                                                            width = 1))
  hits <- findOverlaps(vi_gr,binned_dna_mm10_gr)


  vi_bins <- suppressMessages(do.call(rbind,(apply(viState,2, function(x){
    data.frame(state = x, bins = hits@to) %>% dplyr::group_by(.data$bins) %>%
      dplyr::summarise(nstate2 = sum(.data$state ==2),
                       nstate1 = sum(.data$state ==1))
  }))))
  vi_bins %>% dplyr::mutate(binState =
                              dplyr::if_else(.data$nstate2>.data$nstate1,
                                             "state2",
                                             "state1")) %>%
    dplyr::group_by(.data$bins) %>%
    dplyr::summarise(binSegRatio = sum(.data$binState=="state2")/dplyr::n(),
                     nstate2 =sum(.data$binState=="state2"),
                     nstate1 = sum(.data$binState=="state1")) %>%
    dplyr::mutate(Chrom = chr)
}

#' perSegChrQC
#'
#' Plots the summary statistics of segments that are generated by `sscocaller`
#' \url{https://gitlab.svi.edu.au/biocellgen-public/sscocaller}
#' which have been detected by finding consequtive viter states along the list
#' of SNP markers.
#'
#' It provides guidance in filtering out close double crossovers that are not
#' likely biological but due to technical reasons as well as crossovers that are
#' supported by fewer number of SNPs at the ends of the chromosomes.
#'
#' @param sampleName, the name of the sample to parse which is used as prefix
#' for finding relevant files for the underlying sample
#' @param chroms, the vector of chromosomes
#' @param path, the path to the files, with name patterns *{chrom}_vi.mtx,
#' *{chrom}_viSegInfo.txt, end with slash
#' @param barcodeFile, defaults to NULL, it is assumed to be in the same
#' directory as the other files and with name sampleName_barcodes.txt
#' @param maxRawCO, if a cell has more than `maxRawCO` number of raw crossovers
#' called across a chromosome, the cell is filtered out#'
#' @author Ruqian Lyu
#' @return Histogram plots for statistics summarized across all Viterbi state
#' segments
#' @examples
#' demo_path <- system.file("extdata",package = "comapr")
#' s1_rse_qc <- perSegChrQC(sampleName="s1",
#'                             chroms=c("chr1"),
#'                             path=paste0(demo_path,"/"), maxRawCO=10)
#'
#' @export
#'
perSegChrQC <- function(sampleName,chroms = c("chr1","chr7","chr15"),
                        path,barcodeFile = NULL,maxRawCO=10){
  logllRatio <- nSNP <- bpDist <- NULL
  if (is.null(barcodeFile)) {
    barcodeFile <- paste0(path, sampleName, "_barcodes.txt")
  }
  stopifnot(file.exists(barcodeFile))

  barcodes <- read.table(file = barcodeFile,
                         stringsAsFactors = FALSE,
                         col.names = "barcodes")

  segInfo_list <- bplapply(chroms, function(chr) {

    # snpAnno <- read.table(file = paste0(path, sampleName,
    #   "_", chr, "_snpAnnot.txt"), stringsAsFactors = F,
    #   header = T)
    segInfo <- read.table(file = paste0(path, sampleName,
                                        "_", chr, "_viSegInfo.txt"),
                          stringsAsFactors = FALSE,
                          col.names = c("ithSperm", "Seg_start", "Seg_end",
                                        "logllRatio", "nSNP", "State"))
    segInfo$Chrom <- chr
    segInfo
  })

  segInfo_chrs <- do.call(rbind,segInfo_list)
  rm(segInfo_list)

  rmCells <-
    segInfo_chrs %>% dplyr::group_by(.data$ithSperm,.data$Chrom) %>%
    dplyr::summarise(nCORaw = length(.data$nSNP)-1) %>%
    dplyr::filter(.data$nCORaw >= maxRawCO) %>%
    dplyr::select(.data$ithSperm)
  rmCells <- unique(rmCells$ithSperm)

  p <- segInfo_chrs %>% dplyr::filter(!.data$ithSperm %in% rmCells)%>%
    dplyr::mutate(bpDist = (.data$Seg_end -.data$Seg_start)) %>%
    ggplot()+geom_point(mapping = aes(x =nSNP,y = bpDist,
                                      color = log10(logllRatio)))+
    scale_x_log10()+scale_y_log10()+scale_color_viridis_c()
 p
}
