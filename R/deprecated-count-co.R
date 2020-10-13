#' Count crossovers - deprecated 
#' count the number of crossovers for each sample, per chromosome, cumulatively
#'
#' @param s_gt genotypes of a sample along markers across chromosomes
#' @param chrs a vector of characters indicating the chromosome locations of these markers
#' @param chrPos a vector of characters indicating the base pair position of these markers
#' @param type "counts", or "bool" to specify whether this function should return
#' vector of cumulative counts or vector of TRUE/FALSE
#' @keywords  internal
#' @return
#' A vector of TRUE/FALSE indicating whether crossover happens between each marker
#' interval
#'
#' @author Ruqian Lyu

count_cos <- function(s_gt, interval_df, chrs, chrPos, type = "bool"){
  
  stopifnot(length(s_gt)==length(chrs))
  stopifnot(length(s_gt)==length(chrPos))
  stopifnot(type == "bool" | type == "counts")
  #  is_rb <- NULL
  interval_df$gt <- s_gt
  
  interval_df$gt_before <- interval_df$gt
  interval_df$chr_before <-  interval_df$chrs
  
  interval_df$gt_before[2:nrow(interval_df)] <- interval_df$gt[1:(nrow(interval_df)-1)]
  interval_df$chr_before[2:nrow(interval_df)] <- interval_df$chrs[1:(nrow(interval_df)-1)]
  
  interval_df$is_rb <- (interval_df$gt != interval_df$gt_before &
                          interval_df$chrs == interval_df$chr_before)
  
  interval_df <- data.frame(interval_df)
  
  
  if(type == "bool")
  {
    rownames(interval_df) <- paste0(interval_df$chrs,"_",interval_df$interval_s,"_",
                                    interval_df$interval_e)
    return (interval_df[,"is_rb",drop =FALSE])
    
  } else {
    interval_df <- interval_df %>% group_by(chrs) %>%
      mutate(cum_rb = cumsum(!is.na(.data$is_rb) & .data$is_rb))
    
    interval_df <- data.frame(interval_df)
    
    rownames(interval_df) <- paste0(interval_df$chrs,"_",interval_df$interval_s,"_",
                                    interval_df$interval_e)
    return (interval_df[,"cum_rb",drop =FALSE])
  }
}
