#' getDistortedMarkers
#'
#' Marker segregation distortion detection using chisq-test
#'
#' @details
#' We expect the genotypes to appear with the frequencies of 1:1 homo_alt:hets.
#' We usechisq.test for finding markers that have genotypes among samples that
#' are significantly different from the 1:1 ratio and report them
#' @param p the expected geneotype ratio in a numeric vector, defaults to c(0.5,0.5)
#' @param adj.method Methods to adjust for multiple comparisons, defaults to "BH"
#' @importFrom stats chisq.test p.adjust
#' @inheritParams countGT
#' @author Ruqian Lyu
#' @return data.frame with each row representing one SNP marker and columns
#' containing the chisq.test results
#' @examples
#'
#' corrected_geno <- correctGT(gt_matrix = GenomicRanges::mcols(snp_geno_gr),
#' ref = parents_geno$ref,alt = parents_geno$alt,fail = "Fail",
#' wrong_label = "Homo_ref")
#' GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
#' corrected_geno <- filterGT(snp_geno_gr, min_markers = 30,min_samples = 2)
#' pvalues <- getDistortedMarkers(GenomicRanges::mcols(corrected_geno))
#' @export

getDistortedMarkers <- function(geno, p = c(0.5,0.5),adj.method="BH"){
  geno_table <- vapply(seq_len(nrow(geno)), function(nthRow){
    list(Het = ifelse(is.na(table(unlist(geno[nthRow,]),useNA = "no")["Het"]),
                      0,
                      table(unlist(geno[nthRow,]),useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(unlist(geno[nthRow,]),useNA = "no")["Homo_alt"]),
                           0,
                           table(unlist(geno[nthRow,]),useNA = "no")["Homo_alt"]))
    },FUN.VALUE = vector(mode = "list",length = 2))

  geno_table <- data.frame(Markers = seq_len(ncol(geno_table)),
                           No.Hets =  as.character(unlist(geno_table["Het",])),
                           No.Homo_alt = as.character(unlist(geno_table["Homo_alt",])),
                           stringsAsFactors = FALSE)

  pvals <- vapply(seq_len(nrow(geno_table)), function(marker){
    ctest <- chisq.test(as.numeric(geno_table[marker,2:3]),
                        p = p)
    ctest$p.value
  }, FUN.VALUE = numeric(1))

  #names(pvals) == geno_table$Markers
  geno_table$Pvals <- pvals
  geno_table$Adj.pvals <- p.adjust(pvals, method = adj.method)

  return(geno_table[order(geno_table$Adj.pvals),])
}
