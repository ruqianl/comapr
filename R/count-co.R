
#' countCOs
#' 
#' Count number of COs within each marker interval
#' COs identified in the interval overlapping missing markers are distributed
#' according to marker interval base-pair sizes. Genotypes encoded with "0" are
#' treated as missing value.
#' 
#' @importFrom dplyr group_by
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom BiocParallel bplapply
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRanges seqnames mcols ranges seqinfo<- mcols<-
#' @importFrom IRanges mergeByOverlaps
#' @importFrom IRanges IRanges ranges start width
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @param geno GRanges object or RangedSummarizedExperiment object with genotype
#' matrix that has SNP positions in the rows and cells/samples in the columns

#' @return 
#' GRanges object or RangedSummarizedExperiment with markers-intervals as rows 
#' and samples in columns, values as the number of COs estimated for each marker 
#' interval

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
            mcols(geno_gr) <- as.matrix( assay(geno))
            colnames(mcols(geno_gr)) <- as.character(colData(geno)$barcodes)
            result_gr <- countCOs_gr(geno_gr)
            stopifnot(colData(geno)$barcodes==colnames(mcols(result_gr)))
            SummarizedExperiment(assay=list(co_count=mcols(result_gr)),
                                 colData = colData(geno),
                                 rowRanges =granges(result_gr) )
})

countCOs_gr <- function(geno){
  stopifnot(!is.null(seqnames(geno)))

    crossover_counts <- bplapply(GenomicRanges::mcols(geno),
                          FUN=function(sid_geno,snp_gr=geno[,0]){

              temp_df <- data.frame(chr= as.character(GenomicRanges::seqnames(snp_gr)),
                                  Pos= IRanges::start(GenomicRanges::ranges(snp_gr)),
                                  GT = as.character(sid_geno),
                                  stringsAsFactors = F)
              to_re <- temp_df %>% dplyr::filter((!is.na(.data$GT)) & 
                                                   (!.data$GT=="0")) %>%
                dplyr::group_by(.data$chr) %>%
                dplyr::mutate(CO = (.data$GT!= dplyr::lag(.data$GT)),
                              Prev = dplyr::lag(.data$Pos)) %>%
                dplyr::filter(.data$CO) %>% dplyr::mutate(coid = cumsum(.data$CO))

              if(nrow(to_re)==0){
                re_df <- data.frame(chr=temp_df$chr,Pos = temp_df$Pos,
                                    crossovers = 0)
                colnames(re_df) <- c(colnames(temp_df)[1:2],
                                     names(sid_geno))
              } else {
                co_gr <- GenomicRanges::GRanges(
                  seqnames = to_re$chr,
                  ranges = IRanges(start = to_re$Prev,
                                   end = to_re$Pos),
                  coid = to_re$coid
                )
                mapped_marker_state <- IRanges::mergeByOverlaps(snp_gr,co_gr)
                mapped_marker_state <- as.data.frame(mapped_marker_state)

                mapped_marker_state <- mapped_marker_state %>%
                  dplyr::group_by(.data$snp_gr.seqnames,.data$coid) %>%
                  mutate(snp_gr.prev = dplyr::lag(.data$snp_gr.start,
                                                  default = dplyr::first(.data$snp_gr.start))) %>%
                  mutate(len_prop = (.data$snp_gr.start-.data$snp_gr.prev)/(unique(.data$co_gr.width)-1))
#                message(paste0("\t sec \t",re_df$Pos))
                
                re_df <- data.frame(chr=temp_df$chr,Pos = temp_df$Pos,
                                    crossovers = 0)
                ## keep non zero counts (TEST)
                ## Finds the first matched row with the same Pos. In case this Pos is
                ## the start SNP for the next interval
                mapped_marker_state <-
                  mapped_marker_state[mapped_marker_state$snp_gr.start !=
                                        mapped_marker_state$snp_gr.prev,]
                
                re_df$crossovers[match(paste0(mapped_marker_state$co_gr.seqnames,
                                              "_",mapped_marker_state$snp_gr.start),
                                      paste0(re_df$chr,"_", re_df$Pos))] <- mapped_marker_state$len_prop
                
               colnames(re_df) <- c(colnames(re_df)[1:2],names(sid_geno))

              }
#              message(paste0("\t third \t",as.character(re_df$chr),"_",re_df$Pos))
#              message(paste0("\t fourth \t",start(ranges(snp_gr))))
              rownames(re_df) <- paste0(as.character(re_df$chr),"_",re_df$Pos)
              re_df[,3,drop =F]
    })
    
        final_df <- do.call(cbind,crossover_counts)
        colnames(final_df) <- names(crossover_counts)
        
        gr <- data.frame(seqnames = sapply( strsplit(rownames(final_df),"_"), `[[`,1),
                  end = as.numeric(sapply(strsplit(rownames(final_df),"_"), `[[`,2))) %>% 
          dplyr::group_by(.data$seqnames) %>%
          mutate(start = dplyr::lag(.data$end,default = dplyr::first(.data$end)))
        
        co_gr <- GenomicRanges::GRanges(
          seqnames = gr$seqnames,
          ranges = IRanges(start =gr$start,
                           end = gr$end-1))
        GenomicRanges::mcols(co_gr) <- final_df
        
       
        co_gr <- GenomeInfoDb::sortSeqlevels(co_gr)
        sort(co_gr)
        co_gr[IRanges::width(ranges(co_gr))!=0,]

}

