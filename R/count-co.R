
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
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges mergeByOverlaps
#' @importFrom IRanges IRanges
#' 
#' @param geno corrected genotype matrix returned by \code{correctGT}
#' @return 
#' data.frame with markers as rows and samples in columns, values as the
#' number of COs estimated for each marker interval

#' @author Ruqian Lyu
countCOs <- function(geno ){
  snp_wg_gr <-  GenomicRanges::GRanges(
    seqnames = sapply(strsplit(rownames(geno),"_"),`[[`,1),
    ranges = IRanges::IRanges(start = 
                                as.numeric(sapply(strsplit(rownames(geno),"_"),`[[`,2)),
                              width = 1)
  )
  foreach(sid = colnames(geno),.combine = "cbind",
          .packages = c("GenomicRanges","dplyr")) %dopar% {
            
            marker_names <- rownames(geno)
            gts <- geno[,sid]
            #gts <- na.omit(gts)
            re_df <- data.frame(chr=sapply(strsplit(marker_names,"_"),`[[`,1),
                                Pos=as.numeric(sapply(strsplit(marker_names,"_"),`[[`,2)),
                                GT = gts,
                                stringsAsFactors = F)
            to_re <- re_df %>% filter(!is.na(GT)) %>% dplyr::group_by(chr) %>%
              dplyr::mutate(CO = GT!= dplyr::lag(GT),
                            Prev = dplyr::lag(Pos)) %>% 
              filter(CO) %>% dplyr::mutate(coid = cumsum(CO))
            
            if(nrow(to_re)==0){
              re_df <- data.frame(chr=re_df$chr,pos = re_df$pos,
                                  crossovers = 0)
              colnames(re_df) <- c(colnames(re_df)[1:2],sid)
              
            } else {
              co_gr <- GenomicRanges::GRanges(
                seqnames = to_re$chr,
                ranges = IRanges(start = to_re$Prev,
                                 end = to_re$Pos),
                coid = to_re$coid
              )
              mapped_marker_state <- mergeByOverlaps(snp_wg_gr,co_gr)
              mapped_marker_state <- as.data.frame(mapped_marker_state)
              
              mapped_marker_state <- mapped_marker_state %>%  
                group_by(snp_wg_gr.seqnames,coid) %>% 
                mutate(snp_wg_gr.prev = dplyr::lag(snp_wg_gr.start,
                                                   default = dplyr::first(snp_wg_gr.start))) %>%
                mutate(len_prop = (snp_wg_gr.start-snp_wg_gr.prev)/(unique(co_gr.width)-1))
              
              re_df <- data.frame(chr=re_df$chr,pos = re_df$Pos,
                                  crossovers = 0)
              ## keep non zero counts (TEST)
              ## Finds the first matched row with the same Pos. In case this Pos is
              ## the start SNP for the next interval
              mapped_marker_state <- 
                mapped_marker_state[mapped_marker_state$snp_wg_gr.start != 
                                      mapped_marker_state$snp_wg_gr.prev,]
              re_df$crossovers[match(mapped_marker_state$snp_wg_gr.start,
                                     re_df$pos)] <- mapped_marker_state$len_prop
              colnames(re_df) <- c(colnames(re_df)[1:2],sid)
              
              
            }
            #   message(sid)
            rownames(re_df) <- paste0(re_df$chr,"_",re_df$pos)
            re_df[,3,drop =F]
          }
}
