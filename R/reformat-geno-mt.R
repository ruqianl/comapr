
## Functions for formatting input genotype data.frame including `.label_gt`,
## `.correct_gt`, `.correct_ref` `.change_missing` and a simple imputation
## for missing value method `fill_failed`

#' @import methods
NULL

#' `label_gt` for changing genotypes in alleles format to labels
#'
#' It turns a vector of Genotypes to a vector of Labels consist of
#' `Homo_ref`, `Homo_alt`, and `Het` given the known genotypes for reference
#' and alternative strains.
#' @param s_gt
#' s_gt, a vector of genotypes for one sample across markers
#' @param ref
#' ref, a vector of genotypes for reference strain across markers
#' @param alt
#' alt, a vector of genotypes for alternative strain across markers
#' @param failed
#' what was used for encoding failed genotype calling such as "Fail" in example
#' @return
#' a vector of labels \code{Homo_ref, Homo_alt, Het} indicating the progeny's
#' genotypes across markers
#' @details
#' This function takes the a sample's genotype across each SNP marker in alleles
#' and compare with genotypes of in-bred reference and alternative strains to.
#' If the sample's genotype for a particular SNP marker is the same with the
#' reference strain, it is labelled as Homo_ref \code{homogeneous reference} for
#' a particular SNP marker; if the sample's genotype is the same with the
#' alternative strain it is labelled as Homo_alt \code{homogeneous alternative}
#' for a particular SNP marker; if the sample's genotype is heterozygous then it
#' is labeled as Het \code{heterozygous} for this particular genotypes. If it
#' does not fall in any of the three cases, it is labelled as the string
#' specified by the argument `missing`.
#'
#' Note that the wrong/failed genotype is labelled as the string in `missing`
#' after this function.
#' If there is a different label for failed genotype, provide the label using
#' the `missing` argument.
#'
#' @keywords internal
#' @author Ruqian Lyu

.label_gt <- function(s_gt,ref,alt,failed = "Fail"){
  stopifnot(length(s_gt) == length(ref))
  stopifnot(length(ref) == length(alt))

  ## initialise the vector of GT with `missing`
  tem <- rep(failed,length(s_gt))

  ## if the genotype is the same as reference:
  tem[s_gt == ref] <- "Homo_ref"
  ## if the genotype is the same as alternative:
  tem[s_gt == alt] <- "Homo_alt"
  ## het1 GC
  het1 <- paste0(strtrim(ref,1),strtrim(alt,1))
  ## het2 CG
  het2 <- paste0(strtrim(alt,1),strtrim(ref,1))

  ## if the genotype is the same as het1 or het2:
  tem[s_gt == het1 | s_gt == het2] <- "Het"

  ## the wrong GTs are GTs that can not be converted to Home_ref, Home_alt, Het
  ## or Fail
  ## the wrong GT will be converted to `Fail`

  return(tem)
}



#' change SNPs with genotype 'Fail' to \code{NA}
#'
#' @param s_gt, a column of labelled genotypes
#' @param missing, the string used for encoding missing values default to
#' \code{Fail}
#' @return
#' a vector of genotypes with Fail substituted by `NA`
#' @details

#' calculation.
#'
#' @keywords internal
#' @author  Ruqian Lyu

.change_missing <- function(s_gt, missing = "Fail"){

  if(! any(s_gt %in% missing)){
    return(s_gt)
  } else {

    missing_index  <-  which(s_gt %in% missing)

    s_gt[missing_index] <-  NA
  }
  return(s_gt)
}


#' Infer the genotype of failed SNPs

#' If we have a \code{Fail} in the genotype data and the \code{Fail} in a block
#' of either Home_alt, or Het, we fill in the \code{Fail}s using values of the
#' ones adjacent to it, otherwise they remain as "Fail" to indicate missing
#' values.
#'
#' @param s_gt, a column of labelled genotypes
#' @param fail, the string that is used for encoding failed genotype results,
#' default to \code{Fail}
#' @param chr, the factor vector indicating which chromosomes the markers are
#' on, default to \code{NULL} which means the input markers are all on the same
#' chromosome.
#'
#' @return
#' a vector of genotypes with Failed genotype imputed or changed to `NA` if not
#' imputable
#' @keywords internal
#' @author  Ruqian Lyu

fill_fail <- function(s_gt,fail = "Fail",chr = NULL){
  stopifnot(length(s_gt) >=3)

  if(! any(s_gt == fail)){
    return(s_gt)
  } else if(is.null(chr)) {
    fail_index  <-  grep(fail,s_gt)
    before_index <- fail_index -1
    after_index <- fail_index +1
    if(any(before_index <1)){
      lo_i <- fail_index[which(before_index<1)]
      fail_index <- fail_index[-which(before_index<1)]
      before_index <- fail_index -1
      after_index <- fail_index + 1
    }
    if(any(after_index >= length(s_gt))){
      lo_i <- fail_index[which(after_index > length(s_gt))]
      fail_index <- fail_index[-which(after_index > length(s_gt))]
      before_index <- fail_index -1
      after_index <- fail_index + 1
    }

    able_infer_index <- fail_index[which(s_gt[before_index] == s_gt[after_index]
                                         & s_gt[before_index] != fail)]
    s_gt[able_infer_index] <- s_gt[able_infer_index-1]

    return(s_gt)

  } else { # chr is not null, we should do this for each group of markers
    stopifnot(length(chr) == length(s_gt))
    chrs <- levels(chr)
    names(s_gt) <- paste0("m",seq_len(length(s_gt)))
    to_return <- lapply(chrs, function(s_gt_chr){
      fill_fail(s_gt[chr==s_gt_chr],chr =NULL)
    })
    to_return <- unlist(to_return)
    to_return <- to_return[ paste0("m",seq_len(length(s_gt)))]
    names(to_return) <- NULL
    return(to_return)
  }
}


#' correctGT
#'
#' function for formatting and correction genotypes of markers
#'
#' This function changes genotype in alleles to genotype labels, change Homo_ref
#' to Hets/Fail, infer Failed genotype, and change "Failed" to NA for
#' counting crossover later
#'
#'
#' @param gt_matrix
#' the input genotype matrix of markers by samples with rownames as marker IDs
#' and column names as sample IDs
#' @param ref
#'  a vector of genotypes of the inbred reference strain
#' @param alt
#'  a vector of genotypes of the inbred alternative strain
#' @param failed
#' what was used for encoding failed genotype calling such as "Fail" in example
#' data \code{snp_geno}
#' @param wrong_label
#' what would be considered a wrong genotype label for example Homo_ref which
#' should not be in the possible genotypes of BC1F1 samples
#'
#'
#' @details
#' This function changes genotype in alleles to labels by calling internal
#' functions \code{lable_gt}, and changes the wrong genotype Homo_ref to Fail by
#' calling \code{.change_missing}.
#'
#' @return
#' a genotype data.frame of sample genotypes with dimension as the input
#' `gt_matrix` with genotypes converted to labels and failed calls are changed
#' to NA.
#'
#' @examples
#' data(snp_geno_gr)
#' data(parents_geno)
#' snp_gt_crt <- correctGT(gt_matrix = GenomicRanges::mcols(snp_geno_gr),
#'                       ref = parents_geno$ref,
#'                       alt = parents_geno$alt,
#'                       fail = "Fail",
#'                       wrong_label = "Homo_ref")
#' @export
#'
#' @author Ruqian Lyu

correctGT <- function(gt_matrix, ref, alt,
                      failed = 'Fail', wrong_label = "Homo_ref")
{
  stopifnot(length(ref) == length(alt))
  stopifnot(length(ref) == dim(gt_matrix)[1])

  ### change  GT to labels, non matched GTs are changed to Fail
  row_names <- rownames(gt_matrix)

  gt_matrix <- apply(as.matrix(gt_matrix),2, .label_gt,
                     ref = ref,
                     alt = alt,
                     failed = failed)


  ## change SNPs with genotype `fail` to NA
  gt_matrix <- apply(as.matrix(gt_matrix),2, .change_missing,
                     missing = failed)
  ## Homo_ref calls are changed to Fail since it is not right
  gt_matrix <- apply(as.matrix(gt_matrix),2, .change_missing,
                     missing = wrong_label)
  rownames(gt_matrix) <- row_names

  gt_matrix

}

