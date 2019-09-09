
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(chromoMap)
  library(UpSetR)
})

## Load data sheet


mutate_inform1 <- read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                            sheet = "mutant_25-6-19_informative")
mutate_inform2<- read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                            sheet = "mutant_30-08-2019_informative")

markers_1 <- paste0(mutate_inform1$CHR,"_",mutate_inform1$POS)

markers_2 <- paste0(mutate_inform2$CHR,"_",mutate_inform2$POS)


## same markers
listInput <- list(markers_1 = markers_1, markers_2 =markers_2)

upset(fromList(listInput), order.by = "freq")

mutate_inform_merge <- merge.data.frame(mutate_inform1,mutate_inform2,all = TRUE)

mutate_inform_merge$CHR_POS <- paste0(mutate_inform_merge$CHR,"_",mutate_inform_merge$POS)
## Example for mutate_inform_merge

gt_matrix <- mutate_inform_merge[,grep("[0-9].*[0-9]$",colnames(mutate_inform_merge))]

## sample IDs

### change  GT to labels

gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                   ref = mutate_inform_merge$`C57BL/6J`,
                   alt = mutate_inform_merge$`FVB/NJ [i]`)



### Correct ref genotypes

gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref)


### Fill in Fail or put NA in Fail
gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail)

gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing)

colnames(gt_matrix) <- paste0("Sample_",colnames(gt_matrix))
full_matrix <- cbind(mutate_inform_merge[,c(1:5,ncol(mutate_inform_merge))],
                     gt_matrix)

count.samples <- rowSums(is.na(full_matrix))

full_matrix$na_rate <- count.samples
### Add annotation back to the matrix


head(full_matrix)

#missing_rate <- apply(gt_matrix,1, is.na)

ggplot(data = full_matrix)+geom_point(mapping = aes(x = POS, y = na_rate))+
  facet_wrap(.~CHR)+geom_text(data=full_matrix[full_matrix$na_rate > (0.5*(ncol(full_matrix)-7)),],mapping =aes(x = POS, y = na_rate,label = rsID))

gt_matrix_co_counts <- apply(gt_matrix,2, count_cos,
                             chrs = mutate_inform_merge$CHR,
                             type = "counts")
gt_matrix_co <-apply(gt_matrix,2, count_cos,
                     chrs = mutate_inform_merge$CHR)
#gt_matrix_co



df  = cbind(mutate_inform_merge[,c(1:5,ncol(mutate_inform_merge))], gt_matrix_co_counts)
df = data.frame(df)

head(df)

# ggplot(data = df)+
#   geom_point(mapping = aes(x = POS,y = Sample_56))+facet_wrap(~CHR)




### Calculate recombination rate



gt_matrix_co <- cbind(mutate_inform_merge[,"CHR_POS"], gt_matrix_co)
gt_matrix_co_df = melt(gt_matrix_co,id.vars = "CHR_POS")
head(gt_matrix_co_df)

gt_matrix_dst <- gt_matrix_co_df %>%  group_by(CHR_POS) %>% summarise(t_counts = sum(value == TRUE,na.rm = TRUE),
                                                                      f_counts = sum(value == FALSE,na.rm = TRUE),
                                                                      co_rate = t_counts/(f_counts+t_counts),
                                                                      haldane = -0.5*log(1-2*co_rate),
                                                                      kosambi = 0.25*log( (1+2*co_rate)/(1-2*co_rate)))


head(gt_matrix_dst)

### remove NaN markers? how about markers with Inf
gt_matrix_dst = gt_matrix_dst[!is.nan(gt_matrix_dst$haldane),]

gt_matrix_dst$CHR =  sapply(strsplit(gt_matrix_dst$CHR_POS,"_"),`[[`,1)

gt_matrix_dst$POS =  as.numeric(sapply(strsplit(gt_matrix_dst$CHR_POS,"_"),`[[`,2))

gt_matrix_dst$CHR  <-  factor(gt_matrix_dst$CHR,levels = c(seq(1:19),"X"))
gt_matrix_dst <- gt_matrix_dst[order(gt_matrix_dst$CHR, xtfrm(gt_matrix_dst$POS)), ]

gt_matrix_dst <- gt_matrix_dst %>% group_by(CHR) %>% mutate(cum_haldane = cumsum(haldane),
                                                            cum_kosambi = cumsum(kosambi))

#gt_matrix_dst[sort(gt_matrix_dst$CHR,gt_matrix_dst$POS),]




## genetic length plot



ggplot(data = gt_matrix_dst)+geom_point(mapping = aes(x = POS, y = cum_haldane ),color = "blue")+
  facet_wrap(.~CHR)+geom_point(mapping = aes(x = POS, y = cum_kosambi ),color = "red")+ggtitle("cum_haldane (blue),cum_kosambi (red)")+
  theme(axis.text.x = element_blank())+ylab(label = "cumulative genetic length")


sum(gt_matrix_dst$haldane[!is.infinite(gt_matrix_dst$haldane)],na.rm = TRUE)
sum(gt_matrix_dst$kosambi[!is.infinite(gt_matrix_dst$kosambi)],na.rm = TRUE)

### Chrom Maps


anno_file <- data.frame(marker_ID =gt_matrix_dst$CHR_POS, CHR = gt_matrix_dst$CHR,
                        Start = gt_matrix_dst$POS, End = gt_matrix_dst$POS,
                        co_rate =  100*gt_matrix_dst$co_rate)
anno_file$End[1:(nrow(anno_file)-1)] <- anno_file$Start[2:nrow(anno_file)]

head(anno_file)
tail(anno_file)

write.table(anno_file,file = "data/annotataion_file.txt",sep = "\t",
            quote = FALSE,col.names = FALSE,row.names = FALSE)

cols <- colorRamp(c("blue", "red"))

chromoMap(ch.files = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/chromInfo.txt","data/annotataion_file.txt",
          data_based_color_map  = T,
          data_type = "numeric",
          data_colors = list(c("blue","red")),
          canvas_width = 2200,
          canvas_height = 1100,legend = T,
          chr_width = 10,
          chr_length = 6,
          chr_color = c("orange"), ch_gap = 6)

