
#' countCOs
#' 
#' Count number of COs within each marker interval
#' COs identified in the interval overlapping missing markers are distributed
#' according to marker interval base-pair sizes
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
#' @export
#' @param geno GRanges object with SNP positions and genotypes codes
#' @param BPPARAM the registered backend for parallel computing
#' across samples

#' @return 
#' data.frame with markers as rows and samples in columns, values as the
#' number of COs estimated for each marker interval

#' @author Ruqian Lyu
countCOs <- function(geno,BPPARAM=BiocParallel::bpparam()){
  stopifnot(!is.null(seqnames(geno)))

    
    #   snp_wg_gr <-  GenomicRanges::GRanges(
    # seqnames = sapply(strsplit(rownames(geno),"_"),`[[`,1),
    # ranges = IRanges::IRanges(start = 
    #                             as.numeric(sapply(strsplit(rownames(geno),"_"),`[[`,2)),
    #                           width = 1))
  
    crossover_counts <- bplapply(GenomicRanges::mcols(geno),
                          FUN=function(sid_geno,snp_gr=geno[,0]){

              temp_df <- data.frame(chr= as.character(GenomicRanges::seqnames(snp_gr)),
                                  Pos= IRanges::start(GenomicRanges::ranges(snp_gr)),
                                  GT = as.character(sid_geno),
                                  stringsAsFactors = F)
              to_re <- temp_df %>% dplyr::filter(!is.na(.data$GT)) %>%
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
                re_df$crossovers[match(mapped_marker_state$snp_gr.start,
                                       re_df$Pos)] <- mapped_marker_state$len_prop
                
               colnames(re_df) <- c(colnames(re_df)[1:2],names(sid_geno))

              }
#              message(paste0("\t third \t",as.character(re_df$chr),"_",re_df$Pos))
#              message(paste0("\t fourth \t",start(ranges(snp_gr))))
              rownames(re_df) <- paste0(as.character(re_df$chr),"_",re_df$Pos)
              re_df[,3,drop =F]
    },BPPARAM = BPPARAM)
    
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
