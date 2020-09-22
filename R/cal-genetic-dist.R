#' calGeneticDist
#' 
#' Calculate genetic map length, 
#'
#' Given whether crossover happens in each marker interval, calculate the
#' recombination fraction in samples and then derive the Haldane or Kosambi
#' genetic distances via mapping functions
#'
#' @param co_geno
#' data.frame, returned by \code{detectCO}
#' @examples
#' or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                      alt = snp_geno$FVB.NJ..i.,
#'                      chr = snp_geno$CHR)
#' co_geno <- detectCO(cr_geno,
#'                     chrs = snp_geno$CHR,
#'                     chrPos = snp_geno$POS)
#' co_geno[3,] <- rep(TRUE,dim(co_geno)[2])
#'
#'
#' rb_rate <- calGeneticMap(co_geno)
#' @importFrom plotly summarise
#' @importFrom rlang .data
#' @return data.frame
#' data.frame for all markers with Haldane and Kosambi morgans calculated
#' @export

calGeneticDist <- function(co_geno){
  stopifnot(!is.null(rownames(co_geno)))
  stopifnot(!is.null(colnames(co_geno)))
  ## NEED TO UPDATE To accommadate the new countCOs function **##
  rb_geno <- data.frame("t_counts" = rowSums(co_geno,na.rm = TRUE),
                        "total_na" = rowSums(is.na(co_geno)),
                        "f_counts" = rowSums(co_geno==FALSE,
                                             na.rm = TRUE))
  
  rb_geno_add <- data.frame( total_calls = (rb_geno$f_counts +
                                              rb_geno$t_counts),
                             
                             total_samples = dim(co_geno)[2])
  rb_geno <- cbind(rb_geno,rb_geno_add)
  
  rownames(rb_geno) <- rownames(co_geno)
  
  #stopifnot(rownames(co_geno)==rownames(rb_geno))
  
  # rb_geno <- cbind(co_geno,rb_geno)
  # rb_geno <- data.frame(rb_geno)
  # rownames(rb_geno) <- rownames(co_geno)
  
  # gt_matrix_dst <-
  #   gt_matrix_co_by_marker %>%  group_by(interval_ID) %>% summarise(
  #     t_counts = sum(Cross_over == TRUE, na.rm = TRUE),
  #     f_counts = sum(Cross_over == FALSE, na.rm = TRUE),
  #     total_calls = (f_counts +
  #                      t_counts),
  #     total_na = sum(is.na(Cross_over)),
  #     total_samples = length(Cross_over)
  #   )
  
  if (any(rb_geno$total_calls == 0)) {
    message(paste0(
      sum(rb_geno$total_calls == 0),
      " marker(s) do not have calls (failed) across all samples, they will be removed"
    ))
  }
  
  rb_geno <- rb_geno[rb_geno$total_calls != 0, ]
  #
  #   rb_geno_f <- rb_geno %>%  mutate(na_rate = .data$total_na / .data$total_samples,
  #            pointEst = .data$t_counts / .data$total_calls)
  
  rb_geno_f <- data.frame(na_rate = rb_geno$total_na / rb_geno$total_samples,
                          pointEst = rb_geno$t_counts / rb_geno$total_calls)
  
  rownames(rb_geno_f) <- rownames(rb_geno)
  rb_geno_f <- cbind(rb_geno, rb_geno_f)
  
  if (any(rb_geno_f$pointEst >= 0.5)) {
    warning(
      paste0(
        sum(rb_geno_f$pointEst >= 0.5),
        " markers have cross-over fraction larger or equal to 0.5,
        please check whether they are informative: (these markers are removed)",
        paste0(rownames(rb_geno_f[rb_geno_f$pointEst >= 0.5,]), collapse = ",")
      )
    )
    
  }
  
  rb_geno_f <- rb_geno_f[rb_geno_f$pointEst < 0.5, ]
  
  rb_geno_e <- data.frame(haldane = -0.5 * log(1 - 2 * rb_geno_f$pointEst),
                          kosambi = 0.25 * log((1 + 2 * rb_geno_f$pointEst) / (1 - 2 * rb_geno_f$pointEst)))
  rb_geno_e <- cbind(rb_geno_f,rb_geno_e)
  
  return(rb_geno_e)
}