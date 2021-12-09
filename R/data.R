#' Markers by genotype results for a group of samples
#'
#' @source Statistics Canada. Table 001-0008 - Production and farm value of
#'  maple products, annual. \url{http://www5.statcan.gc.ca/cansim/}
#' @format A data frame with columns:
#' \describe{
#'  \item{C57BL.6J}{genotype of reference mouse train across markers}
#'  \item{FVB.NJ..i.}{genotype of alternative mouse train across markers}
#'  \item{POS}{SNP marker base-pair location }
#'  \item{CHR}{SNP marker chromosome location}
#'  \item{X100}{a mouse sample }
#'  \item{X101}{a mouse sample }
#'  \item{X102}{a mouse sample }
#'  \item{X103}{a mouse sample }
#'  \item{X104}{a mouse sample }
#'  \item{X105}{a mouse sample }
#'  \item{X106}{a mouse sample }
#'  \item{X107}{a mouse sample }
#'  \item{X108}{a mouse sample }
#'  \item{X109}{a mouse sample }
#'  \item{X110}{a mouse sample }
#'  \item{X111}{a mouse sample }
#'  \item{X112}{a mouse sample }
#'  \item{X113}{a mouse sample }
#'  \item{X92}{a mouse sample }
#'  \item{X93}{a mouse sample }
#'  \item{X94}{a mouse sample }
#'  \item{X95}{a mouse sample }
#'  \item{X96}{a mouse sample }
#'  \item{X97}{a mouse sample }
#'  \item{X98}{a mouse sample }
#'  \item{X99}{a mouse sample }
#'  \item{rsID}{the SNP ID}
#' }
#' @usage data(snp_geno)
"snp_geno"

#' Markers by genotype results for a group of samples
#'
#' @source TBD
#' @format A GRanges object:
#' \describe{
#'  \item{X100}{a mouse sample }
#'  \item{X101}{a mouse sample }
#'  \item{X102}{a mouse sample }
#'  \item{X103}{a mouse sample }
#'  \item{X104}{a mouse sample }
#'  \item{X105}{a mouse sample }
#'  \item{X106}{a mouse sample }
#'  \item{X107}{a mouse sample }
#'  \item{X108}{a mouse sample }
#'  \item{X109}{a mouse sample }
#'  \item{X110}{a mouse sample }
#'  \item{X111}{a mouse sample }
#'  \item{X112}{a mouse sample }
#'  \item{X113}{a mouse sample }
#'  \item{X92}{a mouse sample }
#'  \item{X93}{a mouse sample }
#'  \item{X94}{a mouse sample }
#'  \item{X95}{a mouse sample }
#'  \item{X96}{a mouse sample }
#'  \item{X97}{a mouse sample }
#'  \item{X98}{a mouse sample }
#'  \item{X99}{a mouse sample }
#'  \item{rsID}{the SNP ID}
#' }
#' @usage data(snp_geno_gr)
#'
"snp_geno_gr"

#' Parents' genotype for F1 samples in `snp_geno`
#' @format A data.frame:
#' \describe{
#'  \item{C57BL.6J}{genotype of reference mouse train across markers}
#'  \item{FVB.NJ..i.}{genotype of alternative mouse train across markers}
#'}
#' @usage data(parents_geno)
"parents_geno"


#' RangedSummarizedExperiment object containing the Viterbi states SNP markers
#' for samples from two groups. `colData(twoSamples)` contains the sample group
#' factor.
#' @usage data(twoSamples)
"twoSamples"


#' RangedSummarizedExperiment object containing the crossover counts across
#' samples for the list of SNP marker intervals
#'
#' @usage data(coCount)
"coCount"


