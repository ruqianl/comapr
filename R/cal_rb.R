

#' @import methods
NULL

#' `label_gt` for changing genotypes in alleles format to labels
#'
#' It turns a vector of Genotypes to a vector of Labels consist of
#' `Homo_ref`, `Homo_alt`, and `Het` given the known genotypes for reference
#' and alternative strains.
#'
#'
#' @param s_gt
#' s_gt, a vector of genotypes for one sample across markers
#'
#' @param ref
#' ref, a vector of genotypes for reference strain across markers
#'
#' @param alt
#' alt, a vector of genotypes for alternative strain across markers
#'
#' @return
#' a vector of labels \code{Homo_ref, Homo_alt, Het} indicating the progeny's
#' genotypes across markers
#'
#'
#' @details
#' This function takes the a sample's genotype across each SNP marker in alleles
#' and compare with genotypes of in-bred reference and alternative strains to. 
#' If the sample's genotype for a particular SNP marker is the same with the 
#' reference strain, it is labelled as Homo_ref \code{homogeneous reference} for
#' a particular SNP marker; if the sample's genotype is the same with the 
#' alternative strain it is labelled as Homo_alt \code{homogeneous alternative} 
#' for a particular SNP marker; if the sample's genotype is heterozygous then it
#' is labeled as Het \code{heterozygous} for this particular genotypes. If it 
#' does not fall in any of the three cases, it is labelled as the string specified
#' by the argument `missing`.
#'
#' Note that the wrong/failed genotype is labelled as the string in `missing` 
#' after this function.
#' If there is a different label for failed genotype, provide the label using
#' the `missing` argument.
#'
#' @keywords internal
#' @author Ruqian Lyu

label_gt <- function(s_gt,ref,alt,missing ="Fail"){
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
#'
#' @param change_to
#' By default, the Home_ref is changed to \code{Fail}, but it can be changed to
#' any labels such as \code{Het} if we believe the correct genotype for homo_ref
#' is Het.
#'
#' @keywords internal
#'
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
#'
#' @param s_gt, a column of labelled genotypes
#' @param missing, the string used for encoding missing values default to \code{Fail}
#' @return
#' a vector of genotypes with Fail substituted by `NA`
#' @details
#' After changing genotypes in alleles to genotype labels and correct Homo_ref genotypes,
#' the `missing` genotypes are changes to `NA` for downstream calculation.
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

#' If we have a \code{Fail} in the genotype data and the \code{Fail} in a block of either Home_alt,
#' or Het, we fill in the \code{FAIL} using values of the ones adjacent to it,
#' otherwise they remain as "Fail" to indicate missing values.
#'
#' @param s_gt, a column of labelled genotypes
#' @param fail, the string that is used for encoding failed genotype results, default to \code{Fail}
#' @param chr, the factor vector indicating which chromosomes the markers are on, default to \code{NULL}
#' which means the markers are all on the same chromosome.
#'
#' @return
#' a vector of genotypes with Failed genotype imputed or changed to NA if not imputable
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
#    notable_infer_index <- fail_index[-which(s_gt[before_index] == s_gt[after_index] & s_gt[before_index] != fail)]
#    s_gt[notable_infer_index] <- "missing"
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
  ## NOT USE THIS 
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

#' Calculate genetic map length
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

calGeneticMap <- function(co_geno){
  stopifnot(!is.null(rownames(co_geno)))
  stopifnot(!is.null(colnames(co_geno)))

  #gt_matrix_co
  # rb_geno <- co_geno %>% mutate(t_counts = rowSums(co_geno,na.rm = TRUE),
  #                    total_na = rowSums(is.na(co_geno)),
  #                    f_counts = rowSums(co_geno==FALSE,na.rm = TRUE),
  #                    total_calls = (f_counts + t_counts),
  #                    total_samples = length(co_geno))

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

#' Plot markers with missing genotypes
#'
#' @author Ruqian Lyu
#'
#' @param geno the genotype data.frame of markers by samples
#' @param plot_wg
#' whether to plot all markers across whole genome or just the markers that are
#' ever missing across all samples
#' @param missing
#' the label in the matrix that is used for encoding the missing or failed data
#' @param plot_type
#' whether a 'dot' plot or a 'bar' plot should be drawn
#'
#' @details
#' This functions plots the missing markers in a 'dot' plot or 'bar' plot.
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#'
#' @return
#' a ggplot2 object for markers with missing genotype across samples
#' @export
#'
#' @examples
#' or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' plotMissingGT(cr_geno)


plotMissingGT <- function(geno, missing = "Fail", plot_wg = FALSE,
                          plot_type = "dot"){

  stopifnot(plot_type == "dot" | plot_type == "bar")


  mis_matrix <- apply(geno, 2, function(es){
    is.na(es) | es == missing
    })

  plot_df <- melt(mis_matrix)

  if(plot_wg){
    switch (plot_type,
      dot = ggplot(data = plot_df)+
        geom_point(mapping =
                     aes_string(x = "Var1", colour = "value", y = "Var2"))+
      xlab("markers")+ylab("samples")+
    labs(colour = "is_missing")+
      scale_color_manual(values = c("TRUE" = "red",
                                    "FALSE"= "lightgrey"))+
      theme_classic()+theme(axis.text.x = element_text(angle = -90)),

      bar =  ggplot(data = plot_df)+
      geom_bar(mapping = aes_string( fill = "value", x = "Var1"))+
      xlab("markers")+ylab("samples")+
      labs(fill = "is_missing")+
      scale_fill_manual(values = c("TRUE" = "red",
                                   "FALSE"= "lightgrey"))+
      theme_classic()+theme(axis.text.x = element_text(angle = -90))

    )


      } else {
      remain_m <- plot_df %>% group_by(.data$Var1) %>%
        summarise(no_missing = sum(.data$value)) %>% filter(.data$no_missing >0)

      plot_df <- plot_df[plot_df$Var1 %in% remain_m$Var1,]
      switch (plot_type,
              dot = ggplot(data = plot_df)+
                geom_point(mapping =  aes_string(x = "Var1", colour = "value", 
                                                 y = "Var2"))+
                xlab("markers")+ylab("samples")+
                labs(colour = "is_missing")+
                scale_color_manual(values = c("TRUE" = "red",
                                              "FALSE"= "lightgrey"))+
                theme_classic()+theme(axis.text.x = element_text(angle = -90)),

              bar =  ggplot(data = plot_df)+
                geom_bar(mapping = aes_string( fill = "value", x = "Var1"))+
                xlab("markers")+ylab("samples")+
                labs(fill = "is_missing")+
                scale_fill_manual(values = c("TRUE" = "red",
                                             "FALSE"= "lightgrey"))+
                theme_classic()+theme(axis.text.x = element_text(angle = -90))

      )


  }
}

#' countGT
#' count how many samples have genotypes calls across markers
#' and count how many markers that each individual has called genotypes for
#'
#' @importFrom plotly plot_ly subplot
#' @importFrom gridExtra grid.arrange

#' @author Ruqian Lyu
#'
#' @param geno the genotype data.frame of markers by samples from output of
#' function \code{correctGT}
#'
#' @param plot, it determines whether a plot will be generated, defaults to TRUE
#' @param interactive, it determines whether an interactive plot will be generated
#' @export
#'
#' @return
#' A list of two elements including \code{n_markers} and \code{n_samples}

countGT <- function(geno, plot =TRUE,interactive=FALSE){
  if(plot){
    if(interactive){
      ply1 <- plot_ly(data = data.frame(marker_index =
                                          seq(1:length(rowSums(!is.na(geno)))),
                                No.Samples =  rowSums(!is.na(geno)),
                                marker_ID = rownames(geno)),
              x = ~marker_index, y = ~No.Samples, hoverinfo="text",
              text = ~paste('</br> Marker ID: ',marker_ID),
              name = "No. samples by marker",
              mode = "markers",
              type = "scatter")

      ply2 <- plotly::plot_ly(data = 
                                data.frame(sample_index = 
                                             seq(1:length(colSums(!is.na(geno)))),
                                No.Markers =  colSums(!is.na(geno)),
                                sample_ID = colnames(geno)),hoverinfo="text",
              x = ~sample_index, y = ~No.Markers,
              text = ~paste('</br> Sample ID: ',sample_ID), 
              name = "No. markers by sample",type = "scatter",
              mode = "markers")

      p <- subplot(ply1,ply2)
      return(list(ply = p,n_samples = rowSums(!is.na(geno)),
                  n_markers = colSums(!is.na(geno))))

    } else {

      p1 <- ggplot()+
        geom_point(mapping = aes(x = seq(1:length(rowSums(!is.na(geno)))),
                                 y = rowSums(!is.na(geno))))+
        theme_classic()+
        ylab("Number of samples")+xlab("markers index")+
        ggtitle("No. samples by marker")

      p2 <- ggplot()+
        geom_point(mapping = aes(x = seq(1:length(colSums(!is.na(geno)))),
                                 y = colSums(!is.na(geno))))+
        theme_classic()+
        ylab("Number of markers")+
        xlab("samples index")+
        ggtitle("No. markers by sample")
      p <- grid.arrange(p1, p2, nrow = 1)
    }
    return(list(plot = p,
                n_samples = rowSums(!is.na(geno)),
                n_markers = colSums(!is.na(geno))))
  }

  return(list(n_samples = rowSums(!is.na(geno)),
              n_markers = colSums(!is.na(geno))))
}

#' filterGT
#'
#' Filter markers or samples that have too many missing values
#'
#' @inheritParams countGT
#'
#' @param min_markers the minimum number of markers for a sample to be kept
#' @param min_samples the minimum number of samples for a marker to be kept
#'
#' @details
#' This function takes the \code{geno} data.frame and filter the data.frame by
#' the provided cut-offs.
#'
#' @author Ruqian Lyu
#'
#' @return
#' A filtered genotype matrix
#' @export
#'
filterGT <- function(geno, min_markers = 5, min_samples = 3){

  gt_counts <- countGT(geno,plot = FALSE)
  keep_markers <- gt_counts$n_samples >= min_samples
  keep_samples <- gt_counts$n_markers >= min_markers

  message(paste0( "filter out ",sum(keep_markers==FALSE)," marker(s)"))
  message(paste0( "filter out ",sum(keep_samples==FALSE)," sample(s)"))

  return(geno[keep_markers, keep_samples])

}

#' findDupSamples
#'
#' Find the duplicated samples by look at the number of matching genotypes
#' in all pair-wise samples
#'
#' @inherit countGT
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.text grid.rect gpar grid.circle
#' @importFrom circlize colorRamp2
#' @param threshold the frequency cut-off of number of matching genotypes out of
#' all geneotypes for determining whether the pair of samples
#' are duplicated, defaults to \code{0.99}. NAs are regarded as same genotypes for
#' two samples if they both have NA for a marker.
#'
#' @param plot whether a frequency plot should be generated for paired-samples
#' @param in_text whether text of frequencies should be displayed in the
#' heatmap cells
#' @export
#' @author Ruqian Lyu
#' @return
#' The paris of duplicated samples.
#'
#' @examples
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[,1] <- or_geno[,5]
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' dups <- findDupSamples(cr_geno,plot = TRUE)
findDupSamples <- function(geno, threshold = 0.99, plot =TRUE,
                           in_text =TRUE){

  # similarity.matrix<-apply(gt_matrix,2,function(x)colSums(identical(x,gt_matrix[,1])))
  # diag(similarity.matrix)<-0
  #
  n <- seq_len( ncol(geno) )

  #  Make combinations
  id <- expand.grid( n , n )

  #  Get result
  out <- matrix( colSums( (is.na(geno[ , id[,1] ]) & is.na(geno[ , id[,2] ])) | 
                            (geno[ , id[,1] ] == geno[ , id[,2] ]),
                          na.rm =TRUE) , ncol = length(n) )
  dimnames(out)[1] <- dimnames(geno)[2]
  dimnames(out)[2] <- dimnames(geno)[2]

  out_freq <- out / dim(geno)[1]



  if(plot){
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    if(in_text){
       htmap <- Heatmap(out_freq,
            column_title  = "Samples' pair-wise frequencies of having
            same genotype across all markers",
            col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                            gp = gpar(fontsize = 10))
                } else if(i > j) {
                  grid.rect(x = x, y = y,  width = width, height = height,
                              gp = gpar(col = "grey",
                                        fill = col_fun(out_freq[i, j])))
                } else {
                  grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                            gp = gpar(fontsize = 10))
                }
              }, cluster_rows = FALSE, cluster_columns = FALSE,
            name = "Sample pair-wise\ncorrelation")

    } else {
      htmap <- Heatmap(out_freq,
                       column_title  = "Samples' pair-wise frequencies of having
            same genotype across all markers",
                       col = col_fun,cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       name = "Sample pair-wise\ncorrelation",
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))

    }

  }

  l_out_freq <- out_freq
  l_out_freq[lower.tri(out_freq,diag=TRUE)] <- 0


  dups_index <- which(l_out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(l_out_freq)[x])
  if(plot) dups <- list(htmap,dups)
  return(dups)
}


#' findDupMarkers
#'
#' Find the duplicated markers by look at the number of matching genotypes across 
#' samples
#' 
#' @inheritParams countGT
#' @param threshold the frequency cut-off for determining whether the pair of 
#' markers are duplicated, defaults to \code{0.9}
#' @param in_text whether text of frequencies should be displayed in the
#' heatmap cells
#' @param plot whether a frequency heatmap plot should be generated. If the 
#' number of markers is large, do not plot.
#' @return
#' The paris of duplicated markers.
#'
#' @export
#' @author Ruqian Lyu
#'
#' @examples
#'   or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[,1] <- or_geno[,5]
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                      alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' dups <- findDupMarkers(cr_geno,plot = TRUE)


findDupMarkers <- function(geno, threshold = 0.99, plot =FALSE,
                           in_text = TRUE){

  n <- seq_len( nrow(geno) )

  #  Make combinations
  id <- expand.grid( n , n )

  #  Get result
  out <- matrix( rowSums( geno[  id[,1], ] == geno[ id[,2], ],
                          na.rm =TRUE) , nrow = length(n) )
  dimnames(out)[1] <- dimnames(geno)[1]
  dimnames(out)[2] <- dimnames(geno)[1]

  out_freq <- out / dim(geno)[2]

  if (plot) {

  #  heatmap(out_freq,scale = NULL,margins = c(7, 7),xlab ="Markers Pair-wise frequencies of having\nsame genotype across all samples")

    col_fun <- colorRamp2(c(0, 0.5, 0.8,1), c("blue", "white", "red","darkred"))

    if (in_text) {
      htmap <- Heatmap(out_freq,
                       column_title  = "Frequencies of pair-wise markers for having
                     same genotype across all samples",
                       col = col_fun,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.rect(x = x, y = y, width = width, height = height,
                                   gp = gpar(col = "grey", fill = NA))
                         if(i == j) {
                           grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                                     gp = gpar(fontsize = 10))
                         } else if(i > j) {
                           grid.rect(x = x, y = y,  width = width, height = height,
                                     gp = gpar(col = "grey",
                                               fill = col_fun(out_freq[i, j])))
                         } else {
                           grid.text(sprintf("%.1f", out_freq[i, j]), x, y,
                                     gp = gpar(fontsize = 10))
                         }
                       }, cluster_rows = FALSE, cluster_columns = FALSE,
                       name = "Marker pair-wise\ncorrelation",
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))
    } else {
      htmap <- Heatmap(out_freq,
                       column_title  = "Frequencies of pair-wise markers for having
                     same genotype across all samples",
                       col = col_fun,cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       name = "Marker pair-wise\ncorrelation",
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))

    }

  }

  l_out_freq <- out_freq
  l_out_freq[lower.tri(out_freq,diag=TRUE)] <- 0


  dups_index <- which(l_out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(l_out_freq)[x])
  if(plot) dups <- list(htmap,dups)
  return(dups)

}

#' getDistortedMarkers
#'
#' Marker disortation detection using chisq-test
#'
#' @details
#' We expect the genotypes to appear with the frequenceis of 1:1 homo_alt:hets.
#' We use chisq-test for finding markers that have genotypes among samples that
#' are significantly diffferent from the 1:1 ratio and report them
#'
#' @param p the expected geneotype ratio in a numeric vector, defauls to c(0.5,0.5)
#' @importFrom stats chisq.test p.adjust
#' @inheritParams countGT
#' @author Ruqian Lyu
#'
#' @export

getDistortedMarkers <- function(geno, p = c(0.5,0.5)){
 geno.table <- sapply(rownames(geno), function(marker){
   list(Het = ifelse(is.na(table(geno[marker,],useNA = "no")["Het"]),
                     0,
                     table(geno[marker,],useNA = "no")["Het"]),
        Homo_alt = ifelse(is.na(table(geno[marker,],useNA = "no")["Homo_alt"]),
                          0,
                          table(geno[marker,],useNA = "no")["Homo_alt"]))
 })
 geno.table <- data.frame(Markers = colnames(geno.table),
                          No.Hets =  as.character(unlist(geno.table["Het",])),
                          No.Homo_alt = as.character(unlist(geno.table["Homo_alt",])),
                          stringsAsFactors = FALSE)

 pvals <- sapply(as.character(geno.table$Markers), function(marker){
   ctest <- chisq.test(as.numeric(geno.table[geno.table$Markers==marker,2:3]),
                       p = p)
   ctest$p.value
 })

 #names(pvals) == geno.table$Markers
 geno.table$Pvals <- pvals
 geno.table$Adj.pvals <- p.adjust(pvals, method = "BH")

 return(geno.table[order(geno.table$Adj.pvals),])
}


#' plotGTFreq
#'
#' Function to plot the genotypes for all samples faceted by genotype

#' @importFrom plotly plot_ly subplot
#' @importFrom reshape2 melt
#' @inheritParams countGT
#' @author Ruqian Lyu
#' @param interactive, it determines whether an interactive
#' plot will be generated.
#' @param color_set, the RColorBrewer::brewer.pal color set names 
#' @export
#' @examples
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[1,] <- rep("Fail",dim(or_geno)[2])
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.,
#'                     chr = snp_geno$CHR)
#' ft_gt <- filterGT(cr_geno)
#' plotGTFreq(ft_gt)

plotGTFreq <- function(geno, interactive = FALSE,color_set="Set1"){

  geno.table <- sapply(colnames(geno), function(sample){
    list(Het = ifelse(is.na(table(geno[,sample],useNA = "no")["Het"]),
                      0,
                      table(geno[,sample],useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(geno[,sample],useNA = "no")["Homo_alt"])
                           ,0,
                           table(geno[,sample],useNA = "no")["Homo_alt"]))
  })
  geno.table <- data.frame(samples = colnames(geno.table),
                           Freq.Hets =  as.numeric(unlist(geno.table["Het",]))/nrow(geno),
                           Freq.Homo_alt = as.numeric(unlist(geno.table["Homo_alt",])/nrow(geno)),
                           stringsAsFactors = FALSE)
  pltdf <- melt(geno.table)

  if(interactive){
    ply1 <- plot_ly(pltdf, x=~.data$samples,y=~.data$value,type = "scatter",
                    color = ~.data$variable,mode = "markers",colors = color_set)
    return(ply1)
  } else {

    stplt1 <- ggplot(data = pltdf)+
      geom_point(mapping = aes_string(x = "samples", y = "value",
                               color = "variable"))+
      theme_classic()+
      ylab("Genotype Frequecies for each sample")+labs(color ="Genotype")+
      theme(axis.text.x = element_text(angle =-90))
    return(stplt1)

  }

}




