
#title: "YELN_2019_09_04_Genetic-mapsTest"
#author: "Ruqian Lyu"
#date: "9/4/2019"


#' @import methods
NULL

#' `label_gt` for changing genotypes in alleles to labels
#'
#' It turns a vector of Genotypes to a vector of Labels consist of
#' `Homo_ref`, `Homo_alt`, and `Het` given the genotypes for reference
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
#' alt, a vector of genotypes for alternative starin across markers
#'
#' @return
#' a vector of labels \code{Homo_ref, Homo_alt, Het} indicating the progeny's
#' genotypes across markers
#'
#'
#' @details
#' This function takes the a sample's genotype across each SNP marker in alleles
#' and compare with genotypes of in-bred reference and alternative strains to. If the
#' sample's genotype for a particular SNP marker is the same with the reference
#' strain, it is labelled as Homo_ref \code{homogeneous reference} for a particular
#' SNP marker; if the sample's genotype is the same with the alternative strain
#' it is labelled as Homo_alt \code{homogeneous alternative} for a particular SNP
#' marker; if the sample's genotype is heterozygous then it is labelled as
#' Het \code{heterozygous} for this particular genotypes. If it does not fall in any
#' of the three cases, it is labelled as `missing`.
#'
#' Note that the `Fail` genotype is still labelled as `Fail` after this function.
#' If there is a different label for failed genotype, provide the label using
#' the `fail` argument.
#'
#' @keywords internal
#' @author Ruqian Lyu

label_gt <- function(s_gt,ref,alt,fail ="Fail"){
  stopifnot(length(s_gt) == length(ref))
  stopifnot(length(ref) == length(alt))

  ## initialise the vector of GT with NAs
  tem <- rep("missing",length(s_gt))
  ## keep the Fail in
  tem[s_gt == fail] <- fail
  ## if the genotype is the same as referece:
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
  ## the wrong GT will be converted to Missing

  return(tem)
}





#' Function for correcting wrong genotypes
#'
#' If we see home_ref markers in the samples, then something is wrong.
#' We can change all home_ref genotypes to hets as a rough correction or
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
#' a vector of genotype labels with Home_ref corrected to \code{Het} or \code{Fail}
#' as specified by the argument change_to.

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
#' @param missing, the string used for encoding missing values default to \code{missing}
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


#' Count crossovers
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
  interval_df$gt <- s_gt

  interval_df$gt_before <- interval_df$gt
  interval_df$chr_before <-  interval_df$chrs

  interval_df$gt_before[2:nrow(interval_df)] <- interval_df$gt[1:(nrow(interval_df)-1)]
  interval_df$chr_before[2:nrow(interval_df)] <- interval_df$chrs[1:(nrow(interval_df)-1)]

  interval_df$is_rb <- (interval_df$gt != interval_df$gt_before &
                          interval_df$chrs == interval_df$chr_before)

  # interval_df$is_rb_na_asF <- interval_df$is_rb
  # interval_df$is_rb_na_asF[is.na(interval_df$is_rb_na_asF)] <- FALSE


  #  temp_df2$cum_rb <- cumsum(temp_df2$is_rb )
  interval_df <- data.frame(interval_df)


  if(type == "bool")
  {
    rownames(interval_df) <- paste0(interval_df$chrs,"_",interval_df$interval_s,"_",
                                    interval_df$interval_e)
    return (interval_df[,"is_rb",drop =FALSE])
  } else {
    interval_df <- interval_df %>% group_by(chrs) %>%
      mutate(cum_rb = cumsum(!is.na(is_rb) & is_rb))
    rownames(interval_df) <- paste0(interval_df$chrs,"_",interval_df$interval_s,"_",
                                    interval_df$interval_e)
    return (interval_df[,"cum_rb",drop =FALSE])
  }
}

#' correctGT, function for formatting and correction genotypes of markers
#'
#' This function changes genotype in alleles to genotype labels, change Homo_ref
#' to Hets/Fail, infer Failed genotype, and change "Failed" to NA for
#' counting crossover later
#'
#'
#' @param gt_matrix
#' the input genotype matrix of markers by samples with rownames as marker IDs
#' and column names as sample IDs
#'
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

correctGT <- function(gt_matrix, ref, alt, chr, ref_change_to = "Het",
                      infer_fail = FALSE, fail = 'Fail')
{
  stopifnot(ref_change_to %in% c('Het','Fail'))
  stopifnot(length(ref)==length(alt))
  stopifnot(length(ref)==dim(gt_matrix)[1])

  ### change  GT to labels
  row_names <- rownames(gt_matrix)
  gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                     ref = ref,
                     alt = alt)
  gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref, change_to = ref_change_to )

  ### infer and Fill in Fail or put NA in Fail
  if(infer_fail){
    gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail,
                       fail = "Fail",chr = as.factor(chr))
  }
  gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing,
                     missing = fail)

  rownames(gt_matrix) <- row_names

  gt_matrix

}


#' detectCO
#'
#' detect crossovers between every two markers by detecting if there is a change
#' of genotypes
#' @importFrom dplyr group_by
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#'
#' @param gt_matix
#' corrected genotype matrix returned by \code{correct_gt}
#' @param prefix
#' what prefix to add to sample IDs
#' @param type
#' whether return boolean value indicating whether crossover happens between this
#' marker and its preceding marker or returne (type=count) how many crossovers
#' has been detected from all preceding intevals for each chromosome.
#'
#' @param chrs the chromosomes for markers
#' @examples
#'  or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
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
detectCO <-function(gt_matrix, prefix = "Sample_",
                     chrs, chrPos, type = "bool"){
  row_names <- rownames(gt_matrix)
  colnames(gt_matrix) <- paste0(prefix,colnames(gt_matrix))
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

  gt_matrix_co <- apply(gt_matrix,2, count_cos,
                       chrs = chrs,
                       chrPos = chrPos,
                       interval_df = interval_df,type = type)

  cols_names <- names(gt_matrix_co)

  gt_matrix_co <-  Reduce(cbind,gt_matrix_co)
  colnames(gt_matrix_co) <- cols_names

  stopifnot(nrow(gt_matrix) == nrow(gt_matrix_co))

  # rownames(gt_matrix_co) <- row_names

  return(gt_matrix_co)
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
#'   or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
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
#' @return data.frame
#' data.frame for all markers with Haldane and Kosambi morgans calculated
#' @export

calGeneticMap <- function(co_geno){
  stopifnot(!is.null(rownames(co_geno)))
  stopifnot(!is.null(colnames(co_geno)))

  #gt_matrix_co
  rb_geno <- co_geno %>% mutate(t_counts = rowSums(.,na.rm = TRUE),
                     total_na = rowSums(is.na(.)),
                     f_counts = rowSums(.==FALSE,na.rm = TRUE),
                     total_calls = (f_counts + t_counts),
                     total_samples = length(.))
  rownames(rb_geno) <-rownames(co_geno)
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

  rb_geno_f <- rb_geno %>%  mutate(na_rate = total_na /total_samples,
           pointEst = t_counts / total_calls)
  rownames(rb_geno_f) <- rownames(rb_geno)

           # lower_ci = Hmisc::binconf(t_counts,total_calls,alpha = alpha)[,"Lower"],
           # upper_ci = Hmisc::binconf(t_counts,total_calls,alpha = alpha)[,"Upper"])

#co_rate = t_counts / total_calls,

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

  rb_geno_e <-
    rb_geno_f %>% mutate(haldane = -0.5 * log(1 - 2 * pointEst),
           kosambi = 0.25 * log((1 + 2 * pointEst) / (1 - 2 * pointEst)))
  rownames(rb_geno_e) <- rownames(rb_geno_f)
           # kosambi_lower = 0.25 * log((1 + 2 * lower_ci) / (1 - 2 * lower_ci)),
           # kosambi_upper = 0.25 * log((1 + 2 * upper_ci) / (1 - 2 * upper_ci)))


  return(rb_geno_e)
}

#' Plot markers with missing genotypes
#'
#' @author Ruqian Lyu
#'
#' @param gt_matrix the genotype matrix of markers by samples
#' @param plot_wg
#' whether to plot all markers across whole genome or just the missing markers
#'
#' @param missing
#' the label in the matrix that is used for encoding the missing or failed data
#' @import ggplot2
#' @return
#' a plot for markers with missing genotype across samples
#' @export
plotMissingGT <- function(gt_matrix, missing = "Fail", plot_wg = TRUE){

  mis_matrix <- apply(gt_matrix, 2, function(es){
    is.na(es) | es == missing
    })

  plot_df <- melt(mis_matrix)

  plot_df$CHR <- sapply(strsplit(as.character(plot_df$Var1),"_"),`[[`,1)
  plot_df$POS <- as.numeric(sapply(strsplit(as.character(plot_df$Var1),"_"),`[[`,2))
  plot_df$CHR <- factor(plot_df$CHR,levels = c(seq(1:19),"X","Y"))
  plot_df <- plot_df[order(plot_df$CHR, xtfrm(plot_df$POS)), ]
  plot_df$Var1 <- factor(plot_df$Var1,levels = unique(plot_df$Var1))
  #plot_df <- plot_df[plot_df$value,]
  if(plot_wg){
    ggplot(data = plot_df)+geom_point(mapping = aes(x = Var1, colour = value, y = Var2))+xlab("markers")+ylab("samples")+
    labs(colour = "is_missing")+scale_color_manual(values = c("TRUE" = "red",
                                                              "FALSE"= "lightgrey"))
  } else {

    plot_df <- plot_df[plot_df$value,]
    ggplot(data = plot_df)+geom_point(mapping = aes(x = Var1, colour = value, y = Var2))+xlab("markers")+ylab("samples")+
      labs(colour = "is_missing")+scale_color_manual(values = c("TRUE" = "red",
                                                                "FALSE"= "lightgrey"))
  }
}

#' countGT
#' count how many samples have genotypes calls across all markers
#' and count how many markers that each individual have genotypes
#'
#' @importFrom plotly plot_ly subplot
#' @importFrom gridExtra grid.arrange

#' @author Ruqian Lyu
#' @param gt_matrix, the genotype matrix of marker by samples after correction by
#' funcion \code{correctGT}
#' @param plot, it determines whether a plot will be generated.
#'
#' @export
#'
#' @return
#' A list of two elements including \code{n_markers} and \code{n_samples}

countGT <- function(gt_matrix, plot =TRUE,interactive=FALSE){
  if(plot){
    if(interactive){

      ply1 <- plotly::plot_ly(data = data.frame(marker_index = seq(1:length(rowSums(!is.na(gt_matrix)))),
                                No.Samples =  rowSums(!is.na(gt_matrix)),
                                marker_ID = rownames(gt_matrix)),
              x = ~marker_index, y = ~No.Samples,
              text = ~marker_ID,name = "No. samples by marker",mode = "markers",
              type = "scatter")

      ply2 <- plotly::plot_ly(data = data.frame(sample_index = seq(1:length(colSums(!is.na(gt_matrix)))),
                                No.Markers =  colSums(!is.na(gt_matrix)),
                                sample_ID = colnames(gt_matrix)),
              x = ~sample_index, y = ~No.Markers,
              text = ~sample_ID, name = "No. markers by sample",type = "scatter",
              mode = "markers")

      p <- plotly::subplot(ply1,ply2)
      return(list(ply = p,n_samples = rowSums(!is.na(gt_matrix)),
                  n_markers = colSums(!is.na(gt_matrix))))

    } else {
      # par(mfrow=c(1,2))
      # plot(rowSums(!is.na(gt_matrix)),ylab = "Number of samples",xlab="markers index",
      #      main = "No. samples by marker")
      # plot(colSums(!is.na(gt_matrix)),ylab = "Number of markers",xlab="samples index",
      #      main = "No. markers by sample")
      p1 <- ggplot()+geom_point(mapping = aes(x = seq(1:length(rowSums(!is.na(gt_matrix)))),
                                              y = rowSums(!is.na(gt_matrix))))+theme_classic()+
        ylab("Number of samples")+xlab("markers index")+ggtitle("No. samples by marker")

      p2 <- ggplot()+geom_point(mapping = aes(x = seq(1:length(colSums(!is.na(gt_matrix)))),
                                              y = colSums(!is.na(gt_matrix))))+theme_classic()+
        ylab("Number of markers")+xlab("samples index")+ggtitle("No. markers by sample")
      p <- gridExtra::grid.arrange(p1, p2, nrow = 1)
    }
    return(list(plot = p,
                n_samples = rowSums(!is.na(gt_matrix)),
                n_markers = colSums(!is.na(gt_matrix))))
  }

  return(list(n_samples = rowSums(!is.na(gt_matrix)),
              n_markers = colSums(!is.na(gt_matrix))))
}

#' filterGT
#'
#' Filter markers or samples that have too many missing values.
#' @param gt_matrx the genotype matrix of marker by samples after correction by
#' funcion \code{correctGT}
#' @param min_markers the minimum number of markers for a sample to be kept
#' @param min_samples the minimum number of samples for a marker to be kept
#'
#' @details
#' This function takes the \code{gt_matrix} and subset the matrix by
#' the provided cut-offs.
#' @return
#' A filtered genotype matrix
#' @export
filterGT <- function(gt_matrix, min_markers = 5, min_samples = 3){
  gt_counts <- countGT(gt_matrix,plot = FALSE)
  keep_markers <-gt_counts$n_samples >= min_samples
  keep_samples <- gt_counts$n_markers >= min_markers

  message(paste0( "filter out ",sum(keep_markers==FALSE)," marker(s)"))
  message(paste0( "filter out ",sum(keep_samples==FALSE)," sample(s)"))

  return(gt_matrix[keep_markers, keep_samples])

}

#' findDupSamples
#'
#' Find the duplicated samples by look at the number of matching genotypes
#' between all pairs of samples
#' @param gt_matrix the corrected and filtered genotype matrix of markers by samples
#' @param threshold the frequency cut-off of number of matching genotypes out of
#' all geneotypes for determining whether the pair of samples
#' are duplicated, defaults to \code{0.99}
#' @param plot whether a frequency plot should be generated
#' @export
findDupSamples <- function(gt_matrix, threshold = 0.99, plot =TRUE){

  # similarity.matrix<-apply(gt_matrix,2,function(x)colSums(identical(x,gt_matrix[,1])))
  # diag(similarity.matrix)<-0
  #
  n <- seq_len( ncol(gt_matrix) )

  #  Make combinations
  id <- expand.grid( n , n )

  #  Get result
  out <- matrix( colSums( gt_matrix[ , id[,1] ] == gt_matrix[ , id[,2] ],
                          na.rm =TRUE) , ncol = length(n) )
  dimnames(out)[1] <- dimnames(gt_matrix)[2]
  dimnames(out)[2] <- dimnames(gt_matrix)[2]

  out_freq <- out / dim(gt_matrix)[1]



  if(plot){
    heatmap(out_freq,margins = c(7, 7),xlab = "Samples' pair-wise frequencies of having\nsame genotype across all markers")
  }
  out_freq[lower.tri(out_freq,diag=TRUE)] <- 0


  dups_index <- which(out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(out_freq)[x])

  return(dups)
}


#' findDupMarkers
#'
#' Find the duplicated markers by look at the number of matching genotypes across samples
#' @param gt_matrix the corrected and filtered genotype matrix of markers by samples
#' @param threshold the frequency cut-off of number of matching genotypes out of
#' all geneotypes for determining whether the pair of samples
#' are duplicated, defaults to \code{0.9}
#' @param plot whether a frequency heatmap plot should be generated
#' @export
findDupMarkers <- function(gt_matrix, threshold = 0.99, plot =TRUE){

  # similarity.matrix<-apply(gt_matrix,2,function(x)colSums(identical(x,gt_matrix[,1])))
  # diag(similarity.matrix)<-0
  #
  n <- seq_len( nrow(gt_matrix) )

  #  Make combinations
  id <- expand.grid( n , n )

  #  Get result
  out <- matrix( rowSums( gt_matrix[  id[,1], ] == gt_matrix[ id[,2], ],
                          na.rm =TRUE) , nrow = length(n) )
  dimnames(out)[1] <- dimnames(gt_matrix)[1]
  dimnames(out)[2] <- dimnames(gt_matrix)[1]

  out_freq <- out / dim(gt_matrix)[2]

  if (plot) {
    # ComplexHeatmap::Heatmap(
    #   out_freq,
    #   column_title =
    #     "Markers frequencies of having same genotype across all samples",
    #   top_annotation = ComplexHeatmap::HeatmapAnnotation(chr = sapply(strsplit(
    #     colnames(out_freq), "_"
    #   ), `[[`, 1)),
    #   right_annotation = ComplexHeatmap::rowAnnotation(high_freq = ComplexHeatmap::anno_mark(
    #     at = as.numeric(which(out_freq > threshold, arr.ind = TRUE)[, 1]),
    #     labels = names(which(out_freq >
    #                            threshold, arr.ind = TRUE)[, 1])
    #   )),show_row_names = FALSE,show_column_names = FALSE)
    heatmap(out_freq,scale = NULL,margins = c(7, 7),xlab ="Markers Pair-wise frequencies of having\nsame genotype across all samples")


  }
  out_freq[upper.tri(out_freq,diag=TRUE)] <- 0

  dups_index <- which(out_freq > threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(out_freq)[x])

  return(dups)
}


### Look for distorted segregation patterns

### We expect the genotypes to appear with the frequenceis of 1:1 (homo_alt:hets)


#' getDistortedMarkers
#'
#' Marker disortation detection using chisq-test
#' We expect the genotypes to appear with the frequenceis of 1:1 homo_alt:hets, we
#' use chisq-test for finding markers that have genotypes significantly diffferent from
#' the 1:1 ratio and report them
#'
#' @param gt_matrix the corrected and filtered genotype matrix of markers by samples
#' @export

getDistortedMarkers <- function(gt_matrix){
 geno.table <- sapply(rownames(gt_matrix), function(marker){
   list(Het = ifelse(is.na(table(gt_matrix[marker,],useNA = "no")["Het"]),
                     0,
                     table(gt_matrix[marker,],useNA = "no")["Het"]),
        Homo_alt = ifelse(is.na(table(gt_matrix[marker,],useNA = "no")["Homo_alt"]),
                          0,
                          table(gt_matrix[marker,],useNA = "no")["Homo_alt"]))
 })
 geno.table <- data.frame(Markers = colnames(geno.table),
                          No.Hets =  as.character(unlist(geno.table["Het",])),
                          No.Homo_alt = as.character(unlist(geno.table["Homo_alt",])),
                          stringsAsFactors = FALSE)
 pvals <- sapply(as.character(geno.table$Markers), function(marker){
   ctest <- chisq.test(as.numeric(geno.table[geno.table$Markers==marker,2:3]),p = c(0.5,0.5))
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

#' @author Ruqian Lyu
#' @param gt_matrix, the genotype matrix of marker by samples after correction by
#' funcion \code{correctGT}
#' @param plot, it determines whether a plot will be generated.
#' @export

plotGTFreq <- function(gt_matrix, interactive = FALSE){
  geno.table <- sapply(colnames(gt_matrix), function(sample){
    list(Het = ifelse(is.na(table(gt_matrix[,sample],useNA = "no")["Het"]),
                      0,
                      table(gt_matrix[,sample],useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(gt_matrix[,sample],useNA = "no")["Homo_alt"])
                           ,0,
                           table(gt_matrix[,sample],useNA = "no")["Homo_alt"]))
  })
  geno.table <- data.frame(samples = colnames(geno.table),
                           Freq.Hets =  as.numeric(unlist(geno.table["Het",]))/nrow(gt_matrix),
                           Freq.Homo_alt = as.numeric(unlist(geno.table["Homo_alt",])/nrow(gt_matrix)),
                           stringsAsFactors = FALSE)
  pltdf <- melt(geno.table)

  if(interactive){
    ply1 <- plot_ly(pltdf, x=~samples,y=~value,type = "scatter",
                    color =  ~variable)
    return(ply1)
  } else {

    stplt1 <- ggplot(data = pltdf)+geom_point(mapping = aes(x = samples, y = value,
                                                  color = variable))+theme_classic()+
      ylab("Genotype Frequecies for each sample")
    return(stplt1)

  }

}




