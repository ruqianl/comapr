
#' countCOs
#'
#' Count number of COs within each marker interval
#' COs identified in the interval overlapping missing markers are distributed
#' according to marker interval base-pair sizes. Genotypes encoded with "0" are
#' treated as missing value.
#'
#' @importFrom dplyr group_by lag
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom BiocParallel bplapply
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRanges seqnames mcols ranges seqinfo<- mcols<-
#' @importFrom GenomicRanges gaps granges
#' @importFrom IRanges mergeByOverlaps
#' @importFrom IRanges IRanges ranges start width
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment

#' @export
#' @param geno GRanges object or RangedSummarizedExperiment object with genotype
#' matrix that has SNP positions in the rows and cells/samples in the columns

#' @return
#' GRanges object or RangedSummarizedExperiment with markers-intervals as rows
#' and samples in columns, values as the number of COs estimated for each marker
#' interval
#' @examples
#' data(twoSamples)
#' cocount <- countCOs(twoSamples)
#' @author Ruqian Lyu

setGeneric("countCOs",function(geno)
  standardGeneric("countCOs"))

#'@rdname countCOs
setMethod("countCOs",signature = c(geno="GRanges"),function(geno){
  countCOs_gr(geno)
})

#'@rdname countCOs
setMethod("countCOs",signature = c(geno="RangedSummarizedExperiment"),
          function(geno){
            geno_gr <- rowRanges(geno)
            mcols(geno_gr) <- as.matrix( assays(geno)[["vi_state"]])
            colnames(mcols(geno_gr)) <- as.character(colData(geno)$barcodes)
            result_gr <- countCOs_gr(geno_gr)
            stopifnot(colData(geno)$barcodes==colnames(mcols(result_gr)))
            SummarizedExperiment(assays=list(co_count=mcols(result_gr)),
                                 colData = colData(geno),
                                 rowRanges =granges(result_gr) )
})

#'countCOs_gr
#'
#'count number of crossovers for ranges in the input GRanges object
#'
#'@keywords internal
#'@noRd
#'
countCOs_gr <- function(geno) {
  stopifnot(!is.null(seqnames(geno)))

  crossover_counts <- bplapply(mcols(geno),
    FUN = function(sid_geno, snp_gr = geno[, 0]) {
      gps_snp_gr <- gaps(snp_gr)

      temp_df <-
        data.frame(
          chr = as.character(seqnames(snp_gr)),
          Pos = start(ranges(snp_gr)),
          GT = as.character(sid_geno),
          stringsAsFactors = FALSE
        )
      to_re <-
        temp_df %>% filter((!is.na(.data$GT)) &
                                    (!.data$GT ==
                                       "0")) %>%
        group_by(.data$chr) %>%
        mutate(CO = (.data$GT != lag(.data$GT)),
                      Prev = lag(.data$Pos)) %>%
        filter(.data$CO) %>% mutate(coid = cumsum(.data$CO))

      if (nrow(to_re) == 0) {
        ## no crossovers
        mcols(gps_snp_gr)$crossovers <- 0
      } else {
        co_gr <- GRanges(
          seqnames = to_re$chr,
          ranges = IRanges(start = to_re$Prev,
                           end = to_re$Pos -
                             1),
          coid = to_re$coid
        )

        myHits <- findOverlaps(gps_snp_gr, co_gr)

        co_gap_length <- width(ranges(gps_snp_gr))[myHits@from] + 1
        len_prop <- (co_gap_length)/(width(ranges(co_gr))[myHits@to])

        mcols(gps_snp_gr)$crossovers <- 0
        mcols(gps_snp_gr)$crossovers[myHits@from] <- len_prop
      }
      #  colnames(mcols(gps_snp_gr)) <- names(sid_geno)
      gps_snp_gr
    }
  )
  final_df <- lapply(crossover_counts, mcols)
  final_df <- do.call(cbind,final_df)
  colnames(final_df) <- names(crossover_counts)

  co_gr <- crossover_counts[[1]]
  mcols(co_gr) <- final_df
  sort(co_gr)
  co_gr[width(ranges(co_gr)) != 0 &
          rowSums(as.matrix(final_df))>0, ]

}


