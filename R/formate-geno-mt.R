
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

label_gt <- function(s_gt,ref,alt,missing = "Fail"){
  stopifnot(length(s_gt) == length(ref))
  stopifnot(length(ref) == length(alt))
  
  ## initialise the vector of GT with `missing`
  tem <- rep(missing,length(s_gt))
  
  #tem <- rep("missing",length(s_gt))
  ## keep the Fail in
  ## tem[s_gt == missing] <- missing
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
  
  ## the wrong GTs are GTs that are neither Home_ref, Home_alt, Het or Fail
  ## the wrong GT will be converted to `fail`
  
  return(tem)
}



#' Function for correcting wrong genotypes
#'
#' If we see home_ref markers in the samples, then something is wrong.
#' We can change all home_ref genotypes to Hets count as a rough correction or
#' just change it to Fail.
#'
#' @param s_gt
#' a vector of genotype labels as returned by \code{label_gt}
#' @param change_to
#' By default, the Home_ref is changed to \code{Fail}, but it can be changed to
#' any labels such as \code{Het} if we believe the correct genotype for Homo_ref
#' is Het.
#' @keywords internal
#' @return
#' a vector of genotype labels with Home_ref corrected to \code{Het} or 
#' \code{Fail} as specified by the argument change_to.

correct_ref <- function(s_gt, change_to = 'Fail'){
  
  if(! any(s_gt == "Homo_ref")){
    # if for this sample, there is no Homo_ref called across all marker,
    # nothing needs to be changed.
    return(s_gt)
    
  } else {
    # otherwise, the home_ref is changed to Het or NA.
    ref_index <- grep("ref",s_gt)
    s_gt[ref_index] <- change_to
    
    return(s_gt)
  }
}



#' change SNPs with genotype 'Fail' to \code{NA}
#' @param s_gt, a column of labelled genotypes
#' @param missing, the string used for encoding missing values default to 
#' \code{Fail}
#' @return
#' a vector of genotypes with Fail substituted by `NA`
#' @details
#' After changing genotypes in alleles to genotype labels and correct Homo_ref 
#' genotypes, the `missing` genotypes are changes to `NA` for downstream 
#' calculation.
#'
#' @keywords internal
#' @author  Ruqian Lyu

change_missing <- function(s_gt, missing = "Fail"){
  
  if(! any(s_gt == missing)){
    return(s_gt)
  } else {
    
    missing_index  <-  grep(missing,s_gt)
    
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
      #      s_gt[lo_i] <- "missing"
      fail_index <- fail_index[-which(before_index<1)]
      before_index <- fail_index -1
      after_index <- fail_index + 1
    }
    if(any(after_index >= length(s_gt))){
      lo_i <- fail_index[which(after_index > length(s_gt))]
      #      s_gt[lo_i] <- "missing"
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
    names(s_gt) <- paste0("m",seq(1:length(s_gt)))
    to_return <- lapply(chrs, function(s_gt_chr){
      fill_fail(s_gt[chr==s_gt_chr],chr =NULL)
    })
    to_return <- unlist(to_return)
    to_return <- to_return[ paste0("m",seq(1:length(s_gt)))]
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
#' @param chr
#' the vector of chromosomes that the markers are from
#' @param ref
#'  a vector of genotypes of reference strain
#' @param alt
#'  a vector of genotypes of alternative strain
#' @param ref_change_to
#' change homo_reference genotype calls to, default to "Het". The other option is
#' "Fail".
#' @param infer_fail
#' whether the Failed genotype should be imputed. The ways of imputing Failed
#' genotype is described in \code{\link{fill_fail}},default to FALSE.
#' @param fail
#' what was used for encoding failed genotype calling such as "Fail" in example
#' data \code{snp_geno}
#'
#' @details
#' This function changes genotype in alleles to labels by calling \code{lable_gt},
#' corrects the Homo_reference geneotype calls to Hets or Fail by calling
#' \code{correct_ref}, infer the failed genotype calls by calling \code{fill_fail}
#' and change missing values to NA by calling \code{change_missing}.
#'
#' @return
#' a genotype data.frame of sample dimension as the input gt_matrix with genotypes
#' converted to labels and failed calls are changed to NA.
#'
#' @export
#'
#' @author Ruqian Lyu

correctGT <- function(gt_matrix, ref, alt, chr, ref_change_to = "Fail",
                      infer_fail = FALSE, fail = 'Fail')
{
  stopifnot(ref_change_to %in% c('Het','Fail'))
  stopifnot(length(ref) == length(alt))
  stopifnot(length(ref) == dim(gt_matrix)[1])
  
  ### change  GT to labels, non matched GTs are changed to Fail
  row_names <- rownames(gt_matrix)
  gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                     ref = ref,
                     alt = alt)
  ## Homo_ref calls are changed to Fail since it is not likely 
  gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref, 
                     change_to = ref_change_to )
  
  ### infer and fill in Fail or put NA in Fail
  ### Do not need to do this actually 
  if(infer_fail){
    gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail,
                       fail = "Fail",chr = as.factor(chr))
  }
  
  ## change SNPs with genotype `fail` to NA
  gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing,
                     missing = fail)
  
  rownames(gt_matrix) <- row_names
  
  gt_matrix
  
}
