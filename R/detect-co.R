

#' detectCO, deprecated
#'
#' detect crossovers between every two markers by detecting if there is a change
#' of genotypes
#' @importFrom dplyr group_by
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#'
#' @param geno
#' corrected genotype matrix returned by \code{correct_gt}
#' @param prefix
#' what prefix to add to sample IDs
#' @param chrPos
#' the base-pair positions of SNP markers
#' @param type
#' whether return boolean value indicating whether crossover happens between this
#' marker and its preceding marker or returne (type=count) how many crossovers
#' has been detected from all preceding intevals for each chromosome.
#'
#' @param chrs the chromosomes for markers
#' @examples
#' or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                      alt = snp_geno$FVB.NJ..i.,
#'                      chr = snp_geno$CHR)
#' co_geno <- detectCO(cr_geno,
#'                     chrs = snp_geno$CHR,
#'                     chrPos = snp_geno$POS)
#' count_geno <- detectCO(cr_geno,
#'                        chrs = snp_geno$CHR,
#'                        chrPos = snp_geno$POS,
#'                        type="counts")

#' @return
#' a data.frame of marker intervals by samples with values indicating whether
#' a crossover is detected
#'
#' @export
#' @author Ruqian Lyu
#'
detectCO <-function(geno, prefix = "Sample_",
                    chrs, chrPos, type = "bool"){
  row_names <- rownames(geno)
  colnames(geno) <- paste0(prefix,colnames(geno))
  #
  #   gt_matrix_co_counts <- apply(gt_matrix,2, count_cos,
  #                                chrs = chrs,
  #                                type = "counts")
  
  interval_df <- data.frame(chrs,chrPos,
                            stringsAsFactors = FALSE)
  
  interval_df <- interval_df %>% group_by(chrs) %>%
    mutate(interval_e = chrPos,
           interval_s = c("firstM",
                          chrPos[1:(length(chrPos)-1)]))
  
  interval_df <- data.frame(interval_df,row.names = paste0(chrs,"_",chrPos))
  
  gt_matrix_co <- apply(geno,2, count_cos,
                        chrs = chrs,
                        chrPos = chrPos,
                        interval_df = interval_df,type = type)
  
  cols_names <- names(gt_matrix_co)
  
  gt_matrix_co <-  Reduce(cbind,gt_matrix_co)
  colnames(gt_matrix_co) <- cols_names
  
  stopifnot(nrow(geno) == nrow(gt_matrix_co))
  
  # rownames(gt_matrix_co) <- row_names
  
  return(gt_matrix_co)
}
