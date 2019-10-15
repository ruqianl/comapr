
#title: "YELN_2019_09_04_Genetic-mapsTest"
#author: "Ruqian Lyu"
#date: "9/4/2019"


#' @import methods
NULL

#' Function label_gt for changing genotypes to labels
#'
#' It turns a vector of Genotypes to Labels: Homo_ref, Homo_alt, Het
#' given the known reference and alternative genotypes
#'
#'
#' @param s_gt
#' s_gt, a vector of genotypes for one sample across markers
#'
#' @param ref
#' ref, a vector of genotypes for reference across markers
#'
#' @param alt
#' alt, a vector of genotypes for alternative across markers
#'
#' @return
#' a vector of labels indicating the progeny's genotypes across markers
#' : Homo_ref, Homo_alt, Het
#'
#' @details
#'
#' This function takes the a sample's genotype vector and checks for each SNP marker
#' to see whether it matches with homogeneous reference, or homogeneous alternative
#' or heterozygous genotypes.
#' @keywords internal
#' @author Ruqian Lyu

label_gt <- function(s_gt,ref,alt,fail ="Fail"){
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

  ## the wrong GT will be NAs
  ## the wrong GTs will be neither home_ref, home_alt, Het or Fail

  return(tem)
}





#' Function for correct Markers
#'
#' If we see home ref markers in the samples, then something is wrong.
#' We can change all home_ref to hets
#'
#' @param s_gt
#' a vector of genotype labels as returned by \code{label_gt}
#'
#' @param change_to
#' by default, the Home_ref is changed to \code{Het}, but it can be changed to any labels
#' such as \code{Fail} fi we want to exclude these wrong labels.
#' @keywords internal
#' @return
#' a vector of genotype labels with Home_ref corrected to \code{Het} or \code{Fail} as
#' specified by the argument change_to.

correct_ref <- function(s_gt, change_to = 'Het'){

  if(! any(s_gt == "Homo_ref")){
    # if for this sample, there is no Homo_ref called across all marker,
    # nothing needs to be changed.
    return(s_gt)

  } else {
    # otherwise, the home_ref is changed to Het or NA.
    ref_index = grep("ref",s_gt)
    s_gt[ref_index] <- change_to

    return(s_gt)
  }
}




#' change 'missing' labels to \code{NA}
#'
#' @param s_gt, a column of labelled genotypes
#' @param missing, the string used for encoding missing values default to \code{missing}

change_missing <- function(s_gt, missing = "missing"){

  if(! any(s_gt == missing)){
    return(s_gt)
  } else {

    missing_index = grep(missing,s_gt)

    s_gt[missing_index] = NA
  }
  return(s_gt)
}



#' Function fill in FAIL

#' If we have a \code{FAIL} in the genotype data and the \code{FAIL} in a block of either Home_alt,
#' or Het, we fill in the \code{FAIL} using values of the ones adjacent to it,
#' otherwise we fill in "NA" as missing value.
#'
#' @param s_gt, a column of labelled genotypes
#' @param fail, the string that is used for encoding failed genotype results, default to \code{Fail}
#' @keywords internal

fill_fail <- function(s_gt,fail = "Fail"){
  if(! any(s_gt == fail)){
    return(s_gt)
  } else {

    fail_index = grep(fail,s_gt)
    for (i in fail_index){
      if (s_gt[i-1] == s_gt[i+1] & s_gt[i-1] != fail){
        s_gt[i] = s_gt[i-1]
      } else {
        s_gt[i] = "missing"
      }
    }
    return(s_gt)
  }
}


### Count crossovers

## count crossovers for each sample, per chromosome, cumulatively

## s_gt:: genotypes of a sample along markers on one chromosome
## chr:: a vector of characters indicating the chromosome locations of these markers
## type = "counts"
## type = "bool"

count_cos <- function(s_gt, chrs, type = "bool"){

  stopifnot(length(s_gt)==length(chrs))

  temp_df <- data.frame(marker_ID = paste0("ID_",seq(1:length(s_gt))),
                        gt = s_gt,
                        chrs = chrs,stringsAsFactors = FALSE)


  temp_df$gt_before <- temp_df$gt
  temp_df$chr_before <-  temp_df$chrs
  temp_df$gt_before[2:nrow(temp_df)] <- temp_df$gt[1:(nrow(temp_df)-1)]
  temp_df$chr_before[2:nrow(temp_df)] <- temp_df$chrs[1:(nrow(temp_df)-1)]

  temp_df$is_rb <- (temp_df$gt != temp_df$gt_before & temp_df$chrs == temp_df$chr_before)

  temp_df$is_rb_na_asF <- temp_df$is_rb
  temp_df$is_rb_na_asF[is.na(temp_df$is_rb_na_asF)] <- FALSE

  temp_df <- temp_df %>% group_by(chrs) %>% mutate(cum_rb = cumsum(is_rb_na_asF))
  #  temp_df2$cum_rb <- cumsum(temp_df2$is_rb )
  temp_df <- data.frame(temp_df)
  rownames(temp_df) <- temp_df$marker_ID
  co_counts <- temp_df[temp_df$marker_ID,]
  #
  # co_counts_perchr <- sapply(unique(as.character(temp_df2$chrs)), function(chr){
  #   s_gt_chr <- temp_df2$gt[chrs == chr]
  #   co_counts_chr <- rep(0,length(s_gt_chr))
  #
  #   for (i in c(2:length(s_gt_chr))){
  #     if(identical(s_gt_chr[i],s_gt_chr[i-1]) & !anyNA(c(s_gt_chr[i],s_gt_chr[i-1]))){
  #       co_counts_chr[i] <-  co_counts_chr[i-1]
  #     } else {
  #       if(is.na(s_gt_chr[i-1])){
  #         co_counts_chr[i]  <-  NA
  #       } else if(is.na(s_gt_chr[i-1])){
  #         co_counts_chr[i]  <- NA
  #       } else {
  #         co_counts_chr[i]  <-  co_counts_chr[i-1] +1
  #       }
  #
  #     }
  #
  #   }
  #   co_counts_chr
  # })
  # co_counts <- Reduce(append,co_counts_perchr)
  #
  if(type == "bool")
  {
    return (co_counts$is_rb)
  } else {
    return (co_counts$cum_rb)
  }
}

# count_cos <- function(s_gt, chrs, type = "status"){
#
#   stopifnot(length(s_gt)==length(chrs))
#   temp_df <- data.frame(marker_ID = paste0("ID_",seq(1:length(s_gt))),
#                         gt = s_gt,
#                         chrs = chrs,stringsAsFactors = FALSE)
#
#   temp_df2 <- temp_df[!is.na(temp_df$gt),]
#   temp_df2$gt_before <- temp_df2$gt
#   temp_df2$chr_before <-  temp_df2$chrs
#   temp_df2$gt_before[2:nrow(temp_df2)] <- temp_df2$gt[1:(nrow(temp_df2)-1)]
#   temp_df2$chr_before[2:nrow(temp_df2)] <- temp_df2$chrs[1:(nrow(temp_df2)-1)]
#   temp_df2$is_rb <- (temp_df2$gt !=temp_df2$gt_before & temp_df2$chrs == temp_df2$chr_before)
#   temp_df2 <- temp_df2 %>% group_by(chrs) %>% mutate(cum_rb = cumsum(is_rb))
#   #  temp_df2$cum_rb <- cumsum(temp_df2$is_rb )
#   temp_df2 <- data.frame(temp_df2)
#   rownames(temp_df2) <- temp_df2$marker_ID
#   co_counts <- temp_df2[temp_df$marker_ID,]
#   #
#   # co_counts_perchr <- sapply(unique(as.character(temp_df2$chrs)), function(chr){
#   #   s_gt_chr <- temp_df2$gt[chrs == chr]
#   #   co_counts_chr <- rep(0,length(s_gt_chr))
#   #
#   #   for (i in c(2:length(s_gt_chr))){
#   #     if(identical(s_gt_chr[i],s_gt_chr[i-1]) & !anyNA(c(s_gt_chr[i],s_gt_chr[i-1]))){
#   #       co_counts_chr[i] <-  co_counts_chr[i-1]
#   #     } else {
#   #       if(is.na(s_gt_chr[i-1])){
#   #         co_counts_chr[i]  <-  NA
#   #       } else if(is.na(s_gt_chr[i-1])){
#   #         co_counts_chr[i]  <- NA
#   #       } else {
#   #         co_counts_chr[i]  <-  co_counts_chr[i-1] +1
#   #       }
#   #
#   #     }
#   #
#   #   }
#   #   co_counts_chr
#   # })
#   # co_counts <- Reduce(append,co_counts_perchr)
#   #
#   if(type == "status")
#   {
#     return (co_counts$is_rb)
#   } else {
#     return (co_counts$cum_rb)
#   }
# }





### Test crossovers

# Test crossovers for each sample, per chromosome, whether a crossover happed for the marker


## s_gt:: genotypes of a sample along markers on one chromosome
## chr:: a vector of characters indicating the chromosome locations of these markers
#
# test_cos <- function(s_gt, chrs){
#   stopifnot(length(s_gt)==length(chrs))
#
#   co_counts_perchr <- sapply(unique(chrs), function(chr){
#     s_gt_chr <- s_gt[chrs == chr]
#     co_counts_chr <- s_gt_chr
#
#     for (i in c(2:length(s_gt_chr))){
#       if(identical(s_gt_chr[i],s_gt_chr[i-1]) & !anyNA(c(s_gt_chr[i],s_gt_chr[i-1]))){
#         co_counts_chr[i] <-  FALSE
#       } else {
#         if(is.na(s_gt_chr[i-1])){
#           co_counts_chr[i]  <-  NA
#         } else if(is.na(s_gt_chr[i-1])){
#           co_counts_chr[i]  <- NA
#         } else {
#           co_counts_chr[i]  <- TRUE
#         }
#
#       }
#
#     }
#     co_counts_chr
#   })
#   co_counts <- Reduce(append,co_counts_perchr)
#
#   return (co_counts)
# }
#



#' correctGT
#'
#' map GT to genotype labels, change Homo_ref to Hets, infer missing data,
#' and change "missing" to NA for counting crossover later
#'
#'
#' @param gt_matrix
#' the input from excel with rownames as marker IDs and column names as sample IDs
#'
#' @param ref
#'  is the reference genotype
#' @param alt
#' is the alternative genotype.
#'
#' @details
#' This function changes genotype code to labels by calling \code{lable_gt}, corrects
#' the Home reference geneotype calls to Hets by calling \code{correct_ref}, infer the
#' failed genotype calls by calling \code{fill_fail} and change missing values to NA by
#' calling \code{change_missing}.
#'
#' @export
#'
#' @author Ruqian Lyu

correctGT <- function(gt_matrix, ref, alt)
{

  ### change  GT to labels
  row_names <- rownames(gt_matrix)

  gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                     ref = ref,
                     alt = alt)
  gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref)



  ### Fill in Fail or put NA in Fail
  gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail)

  gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing)

  rownames(gt_matrix) <- row_names

  gt_matrix

}


#' detect crossovers
#'
#' detect crossovers between every two markers by examining change in genotype patterns
#'
#' @param gt_matix
#'
#' correct genotype matrix returned by `correct_gt`
#'
#' @param prefix
#'
#' what prefix to add to sample IDs
#'
#' @param  type
#'
#' whether return boolean value indicating whether crossover happens or integer value
#' indicating how many crossovers happen from all preceding intevals
#'
#' @param  chr
#'
#' the chromosomes for markers on the rows
#'
#' @return
#'
#' a matrix of marker by samples with values indicating crossovers
#' @export
#'
detect_co <-function(gt_matrix, prefix = "Sample_",
                     chrs, type = "bool"){
  row_names <- rownames(gt_matrix)
  colnames(gt_matrix) <- paste0(prefix,colnames(gt_matrix))
  #
  #   gt_matrix_co_counts <- apply(gt_matrix,2, count_cos,
  #                                chrs = chrs,
  #                                type = "counts")

  gt_matrix_co <-apply(gt_matrix,2, count_cos,
                       chrs = chrs, type = type)
  stopifnot(nrow(gt_matrix) == nrow(gt_matrix_co))

  rownames(gt_matrix_co) <- row_names

  return(gt_matrix_co)
}


#' Calculate genetic map length
#'
#' Given whether crossover happens in every two marker interval, calculate the
#' recombination fraction and then derive the Haldane or Kosambi morgans.
#'
#' @param gt_matrix_co_by_marker
#' data.frame, from `melt` the matrix returned by `detect_co`
#'
#' @return data.frame
#' data.frame for all markers with Haldane and Kosambi morgans calculated
#'

cal_geneticMap <- function(gt_matrix_co_by_marker){


  #gt_matrix_co

  gt_matrix_dst <-
    gt_matrix_co_by_marker %>%  group_by(CHR_POS) %>% summarise(
      t_counts = sum(Cross_over == TRUE, na.rm = TRUE),
      f_counts = sum(Cross_over == FALSE, na.rm = TRUE),
      total_calls = (f_counts +
                       t_counts),
      total_na = sum(is.na(Cross_over)),
      total_samples = length(Cross_over)
    )
  if (any(gt_matrix_dst$total_calls == 0)) {
    message(paste0(
      sum(gt_matrix_dst$total_calls == 0),
      " marker(s) do not have calls across all samples, they will be removed"
    ))
  }

  gt_matrix_dst <- gt_matrix_dst[gt_matrix_dst$total_calls != 0, ]
  gt_matrix_dst <-
    gt_matrix_dst %>%  group_by(CHR_POS) %>% mutate(co_rate = t_counts / total_calls,
                                                    na_rate = total_na /
                                                      total_samples)

  if (any(gt_matrix_dst$co_rate >= 0.5)) {
    warning(
      paste0(
        sum(gt_matrix_dst$co_rate >= 0.5),
        " markers have cross-over fraction larger or equal to 0.5,
        please check whether they are informative: ",
        paste0(gt_matrix_dst$CHR_POS[gt_matrix_dst$co_rate >= 0.5], collapse = ",")
      )
    )

  }

  gt_matrix_dst <- gt_matrix_dst[gt_matrix_dst$co_rate < 0.5, ]

  gt_matrix_dst <-
    gt_matrix_dst %>%  group_by(CHR_POS) %>%
    summarise(haldane = -0.5 * log(1 - 2 * co_rate),
              kosambi = 0.25 * log((1 + 2 * co_rate) / (1 - 2 * co_rate)))


  return(gt_matrix_dst)
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
#'
#' @return
#' a plot for markers with missing genotype across samples
#'
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
    theme(axis.text.x = element_blank())+
    labs(colour = "is_missing")+scale_color_manual(values = c("TRUE" = "red",
                                                              "FALSE"= "white"))
  } else {

    plot_df <- plot_df[plot_df$value,]
    ggplot(data = plot_df)+geom_point(mapping = aes(x = Var1, colour = value, y = Var2))+xlab("markers")+ylab("samples")+
      theme(axis.text.x = element_blank())+
      labs(colour = "is_missing")+scale_color_manual(values = c("TRUE" = "red",
                                                                "FALSE"= "white"))
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
#'
#' A filtered genotype matrix

filterGT <- function(gt_matrix, min_markers = 5, min_samples = 3){
  gt_counts <- countGT(gt_matrix)
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
#'
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
    heatmap(out_freq,main ="Samples' frequencies of having same genotype across all markers")
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
#'
findDupMarkers <- function(gt_matrix, threshold = 0.9, plot =TRUE){

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

  out_freq <- out / dim(gt_matrix)[1]



  if(plot){
    heatmap(out_freq,main ="Markers frequencies of having same genotype across all samples")
  }

  out_freq[lower.tri(out_freq,diag=TRUE)] <- 0


  dups_index <- which(out_freq == threshold, arr.ind = TRUE)
  dups <- apply(dups_index,1,function(x)rownames(out_freq)[x])

  return(dups)
}


### Look for distorted segregation patterns

### We expect the genotypes to appear with the frequenceis of 1:1 (homo_alt:hets)

###

#' getDisortationMarkers
#'
#' Marker disortation detection using chisq-test
#' We expect the genotypes to appear with the frequenceis of 1:1 homo_alt:hets, we
#' use chisq-test for finding markers that have genotypes significantly diffferent from
#' the 1:1 ratio and report them
#'

getDisortatedMarkers <- function(gt_matrix){
 geno.table <- sapply(rownames(gt_matrix), function(marker){
   list(Het = ifelse(is.na(table(gt_matrix[marker,],useNA = "no")["Het"]),0,table(gt_matrix[marker,],useNA = "no")["Het"]),
        Homo_alt = ifelse(is.na(table(gt_matrix[marker,],useNA = "no")["Homo_alt"]),0,table(gt_matrix[marker,],useNA = "no")["Homo_alt"]))
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
#' @importFrom plotly plot_ly subplot
#' @importFrom reshape2 melt
#' Function to plot the genotypes for all samples faceted by genotype

plotGTFreq <- function(gt_matrix, interactive = FALSE){
  geno.table <- sapply(colnames(gt_matrix), function(sample){
    list(Het = ifelse(is.na(table(gt_matrix[,sample],useNA = "no")["Het"]),0,table(gt_matrix[,sample],useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(gt_matrix[,sample],useNA = "no")["Homo_alt"]),0,table(gt_matrix[,sample],useNA = "no")["Homo_alt"]))
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

