
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(chromoMap)
})

## Load data sheet


mutate_inform <- read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                            sheet = "mutant_25-6-19_informative")

mutate_inform$CHR_POS <- paste0(mutate_inform$CHR,"_",mutate_inform$POS)
## Example for mutate_inform

gt_matrix <- mutate_inform[,grep("[0-9].*[0-9]$",colnames(mutate_inform))]

## sample IDs

### change  GT to labels

gt_matrix <- apply(as.matrix(gt_matrix),2, label_gt,
                   ref = mutate_inform$`C57BL/6J`,
                   alt = mutate_inform$`FVB/NJ [i]`)



### Correct ref genotypes

gt_matrix <- apply(as.matrix(gt_matrix),2, correct_ref)


### Fill in Fail or put NA in Fail
gt_matrix <- apply(as.matrix(gt_matrix),2, fill_fail)
gt_matrix <- apply(as.matrix(gt_matrix),2, change_missing)

colnames(gt_matrix) <- paste0("Sample_",colnames(gt_matrix))
full_matrix <- cbind(mutate_inform[,c(1:5,ncol(mutate_inform))], gt_matrix)


### Add annotation back to the matrix


head(full_matrix)

gt_matrix_co_counts <- apply(gt_matrix,2, count_cos,
                             chrs = mutate_inform$CHR,
                             type = "counts")
gt_matrix_co <-apply(gt_matrix,2, count_cos,
                     chrs = mutate_inform$CHR)
#gt_matrix_co



df  = cbind(mutate_inform[,c(1:5,ncol(mutate_inform))], gt_matrix_co_counts)
df = data.frame(df)

head(df)

# ggplot(data = df)+
#   geom_point(mapping = aes(x = POS,y = Sample_56))+facet_wrap(~CHR)




### Calculate recombination rate



gt_matrix_co <- cbind(mutate_inform[,"CHR_POS"], gt_matrix_co)
gt_matrix_co_df = melt(gt_matrix_co,id.vars = "CHR_POS")
head(gt_matrix_co_df)

gt_matrix_dst <- gt_matrix_co_df %>%  group_by(CHR_POS) %>% summarise(t_counts = sum(value == TRUE,na.rm = TRUE),
                                                                      f_counts = sum(value == FALSE,na.rm = TRUE),
                                                                      co_rate = t_counts/(f_counts+t_counts),
                                                                      haldane = -0.5*log10(1-2*co_rate),
                                                                      kosambi = 0.25*log10( (1+2*co_rate)/(1-2*co_rate)))

sum(gt_matrix_dst$co_rate*100,na.rm = TRUE)

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

sum(gt_matrix_dst$co_rate*100,na.rm = TRUE)



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

