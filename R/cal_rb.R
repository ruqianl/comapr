#title: "YELN_2019_09_04_Genetic-mapsTest"
#author: "Ruqian Lyu"
#date: "9/4/2019"




## Function label GT

## s_gt, a column of genotypes for one sample across all markers
## s_gt, a column of genotypes for reference across all markers
## s_gt, a column of genotypes for alternative across all markers

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





## Function correct Markers

#If we see home ref markers in the samples, then something is wrong. We change all home_ref to hets or marker as NA




correct_ref <- function(s_gt,change_to = 'Het'){

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

#(test_gt <- label_gt(mutate_inform$`56`,mutate_inform$`C57BL/6J`,mutate_inform$`FVB/NJ [i]`))

#(correct_ref(test_gt))





## Function fill in FAIL

# If we have a `FAIL` in the data and the `FAIL` in a block of either Home_alt, or Het, we fill in the FAIL using values of the ones adjacent to it, otherwise we fill in "NA" as missing value.


## change 'missing' to 'NA'
## s_gt, a column of labelled genotypes

change_missing <- function(s_gt ,missing = "missing"){
  if(! any(s_gt == missing)){
    return(s_gt)
  } else {

    missing_index = grep(missing,s_gt)

    s_gt[missing_index] = NA
  }
  return(s_gt)
}


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






# plot_df <- data.frame(full_matrix)
# plot_df$CHR = factor(plot_df$CHR,levels = c(1:19,"X"))
#
# ggplot(data = plot_df)+
#   geom_point(mapping = aes(x = POS,y = Sample_56))+facet_wrap(~CHR)
#
# ggplot(data = plot_df)+
#   geom_point(mapping = aes(x = POS,y = Sample_57))+facet_wrap(~CHR)+theme(axis.text.x = element_text(face="bold", color="#993333",size=8,
#                                                                                                      angle = 45))+ scale_x_continuous(labels = scales::comma)



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




## Count Crossovers

## gt_matrix is a matrix of marker by samples with values as genotypes such as
## AA,GC,GT etc

## ref is the reference genotype
## alt is the alternative genotype

## the function returns a matrix of marker by crossover happens

## type:: counts, counts of how many crossovers happened until this marker


detect_co <-function(gt_matrix, ref, alt,prefix = "Sample_",
                     chrs, type = "bool"){
  ### change  GT to labels
  row_names <- rownames(gt_matrix)
  gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                     ref = ref,
                     alt = alt)
  gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref)



  ### Fill in Fail or put NA in Fail
  gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail)

  gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing)

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



cal_geneticMap <- function(gt_matrix_co_by_marker){

  #gt_matrix_co

  gt_matrix_dst <- gt_matrix_co_by_marker %>%  group_by(CHR_POS) %>% summarise(t_counts = sum(Cross_over == TRUE,na.rm = TRUE),
                                                                               f_counts = sum(Cross_over == FALSE,na.rm = TRUE),
                                                                               total_calls = (f_counts+t_counts),
                                                                               total_na = sum(is.na(Cross_over)),
                                                                               total_samples = length(Cross_over))
  if(any(gt_matrix_dst$total_calls == 0)){
    message(paste0(sum(gt_matrix_dst$total_calls == 0)," marker(s) do not have calls across all samples, they will be removed"))
  }

  gt_matrix_dst <- gt_matrix_dst[gt_matrix_dst$total_calls != 0,]
  gt_matrix_dst <- gt_matrix_dst %>%  group_by(CHR_POS) %>% mutate(co_rate = t_counts/total_calls,
                                                                   na_rate = total_na/total_samples)

  if(any(gt_matrix_dst$co_rate >= 0.5)){
    message(paste0(sum(gt_matrix_dst$co_rate >= 0.5)," markers have cross-over fraction larger or equal to 0.5,they will be removed"))
  }

  gt_matrix_dst <- gt_matrix_dst[gt_matrix_dst$co_rate <0.5,]

  gt_matrix_dst <- gt_matrix_dst %>%  group_by(CHR_POS) %>% summarise(haldane = -0.5*log(1-2*co_rate),
                                                                      kosambi = 0.25*log( (1+2*co_rate)/(1-2*co_rate)))


  return(gt_matrix_dst)
}
