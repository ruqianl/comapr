## Test functions

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
mutate_inform1 <- data.frame(mutate_inform1)
mutate_inform2<- read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                            sheet = "mutant_30-08-2019_informative")

mutate_inform2 <- data.frame(mutate_inform2)

markers_1 <- paste0(mutate_inform1$CHR,"_",mutate_inform1$POS)

markers_2 <- paste0(mutate_inform2$CHR,"_",mutate_inform2$POS)


## same markers
listInput <- list(markers_1 = markers_1, markers_2 =markers_2)

upset(fromList(listInput), order.by = "freq")

mutate_inform_merge <- merge.data.frame(mutate_inform1,mutate_inform2, all = TRUE)

wt_inform1 <-  read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                          sheet = "wildtype_2-5-19_informative")
wt_inform1 <- data.frame(wt_inform1)


wt_inform2 <-  read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                          sheet = "wt_30-08-2019_informative")
wt_inform2 <- data.frame(wt_inform2)

markers_1 <- paste0(wt_inform1$CHR,"_",wt_inform1$POS)

markers_2 <- paste0(wt_inform2$CHR,"_",wt_inform2$POS)

listInput <- list(markers_1 = markers_1, markers_2 =markers_2)

upset(fromList(listInput), order.by = "freq")

wt_inform_merge <- merge.data.frame(wt_inform1,wt_inform2,all = TRUE)

#############################







jobs <- list("aug_WT" = wt_inform2,
          "may_WT" = wt_inform1,
          "aug_MT" = mutate_inform2,
          "jun_MT" = mutate_inform1,
          "merge_WT" = wt_inform_merge,
          "merge_MT" = mutate_inform_merge)
dist_outs <- list()

names(jobs)
jobID <- 1
names(jobs)[jobID]

for(job in jobs){
#  print(head(job))
  sample_informative <- job


  sample_informative[duplicated(paste0(sample_informative$CHR,"_",sample_informative$POS)),]
  sample_informative$CHR_POS <- paste0(sample_informative$CHR,"_",sample_informative$POS)

  ## Example for sample_informative
  rownames(sample_informative) <- sample_informative$CHR_POS

  gt_matrix <- sample_informative[,grep("[0-9].*[0-9]$",colnames(sample_informative))]

  #rownames(gt_matrix)

  gt_matrix_co <- detect_co(gt_matrix = gt_matrix,
                            ref = sample_informative$C57BL.6J,
                            alt = sample_informative$FVB.NJ..i.,
                            chrs = sample_informative$CHR,
                            type = "bool")

  gt_matrix_co_by_marker <- melt(gt_matrix_co)

  colnames(gt_matrix_co_by_marker) <- c("CHR_POS","Sample","Cross_over")

  gt_matrix_co_by_marker$CHR_POS <- as.character(gt_matrix_co_by_marker$CHR_POS)

  #head(gt_matrix_co_by_marker)

  gt_matrix_dist <- cal_geneticMap(gt_matrix_co_by_marker = gt_matrix_co_by_marker )

  #head(gt_matrix_dist)

  ### Plotting


  gt_matrix_dist$CHR =  sapply(strsplit(gt_matrix_dist$CHR_POS,"_"),`[[`,1)

  gt_matrix_dist$POS =  as.numeric(sapply(strsplit(gt_matrix_dist$CHR_POS,"_"),`[[`,2))

  gt_matrix_dist$CHR  <-  factor(gt_matrix_dist$CHR,levels = c(seq(1:19),"X"))

  ## important to order the markers by chr then by position

  gt_matrix_dist <- gt_matrix_dist[order(gt_matrix_dist$CHR, xtfrm(gt_matrix_dist$POS)), ]

  gt_matrix_dist <- gt_matrix_dist %>% group_by(CHR) %>% mutate(cum_haldane = cumsum(haldane),
                                                                cum_kosambi = cumsum(kosambi))

  #gt_matrix_dst[sort(gt_matrix_dst$CHR,gt_matrix_dst$POS),]


  ## total centimorgans

  k_sum = sum(gt_matrix_dist$kosambi*100)

  ## genetic length plot


#
#   print(ggplot(data = gt_matrix_dist)+
#     facet_wrap(.~CHR)+geom_point(mapping = aes(x = POS, y = cum_kosambi ),color = "red")+ggtitle(paste0(names(jobs)[jobID],"_",k_sum," cum_kosambi "))+
#     theme(axis.text.x = element_blank())+ylab(label = "cumulative genetic length"))


  write.table(gt_matrix_dist,file = paste0("./data/",names(jobs)[jobID],".tsv"))
  gt_matrix_dist$from  <-  names(jobs)[jobID]
  dist_outs[[names(jobs)[jobID]]] <- gt_matrix_dist
  jobID  <-  jobID+1
}

##


## Plots

plot_df <- merge(dist_outs$merge_WT,dist_outs$merge_MT,all = TRUE)
########
head(plot_df)
plot_df <- plot_df[order(plot_df$CHR, xtfrm(plot_df$POS)), ]
head(plot_df)
ggplot(data = plot_df)+geom_point(mapping = aes(x = POS, y= cum_kosambi, colour = from))+facet_wrap(.~CHR)+
  ggtitle(label = paste0("total genetic length Mutant: \n",
                         round(sum(plot_df$kosambi[plot_df$from == "merge_MT"])*100,2)," vs WT",
                         round(sum(plot_df$kosambi[plot_df$from == "merge_WT"])*100 ,2)))
## Run for one job

ggsave("./data/merge.png")
plot_df <- merge(dist_outs$may_WT,dist_outs$jun_MT,all = TRUE)
########
head(plot_df)
plot_df <- plot_df[order(plot_df$CHR, xtfrm(plot_df$POS)), ]
head(plot_df)
ggplot(data = plot_df)+geom_point(mapping = aes(x = POS, y= cum_kosambi, colour = from))+facet_wrap(.~CHR)+
  ggtitle(label = paste0("total genetic length Mutant: \n",
                         round(sum(plot_df$kosambi[plot_df$from == "jun_MT"])*100,2)," vs WT",
                         round(sum(plot_df$kosambi[plot_df$from == "may_WT"])*100 ,2)))
## Run for one job

ggsave("./data/may_jun.png")

plot_df <- merge(dist_outs$aug_WT,dist_outs$aug_MT,all = TRUE)
########
head(plot_df)
plot_df <- plot_df[order(plot_df$CHR, xtfrm(plot_df$POS)), ]
head(plot_df)
ggplot(data = plot_df)+geom_point(mapping = aes(x = POS, y= cum_kosambi, colour = from))+facet_wrap(.~CHR)+
  ggtitle(label = paste0("total genetic length Mutant: \n",
                         round(sum(plot_df$kosambi[plot_df$from == "aug_MT"])*100,2)," vs WT",
                         round(sum(plot_df$kosambi[plot_df$from == "aug_WT"])*100 ,2)))

ggsave("./data/aug.png")
##########
# sample_informative <- wt_inform2
#
#
# sample_informative[duplicated(paste0(sample_informative$CHR,"_",sample_informative$POS)),]
# sample_informative$CHR_POS <- paste0(sample_informative$CHR,"_",sample_informative$POS)
#
# ## Example for sample_informative
# rownames(sample_informative) <- sample_informative$CHR_POS
#
# gt_matrix <- sample_informative[,grep("[0-9].*[0-9]$",colnames(sample_informative))]
#
# #rownames(gt_matrix)
#
# gt_matrix_co <- detect_co(gt_matrix = gt_matrix,
#                           ref = sample_informative$C57BL.6J,
#                           alt = sample_informative$FVB.NJ..i.,
#                           chrs = sample_informative$CHR,
#                           type = "bool")
#
# gt_matrix_co_by_marker <- melt(gt_matrix_co)
#
# colnames(gt_matrix_co_by_marker) <- c("CHR_POS","Sample","Cross_over")
#
# gt_matrix_co_by_marker$CHR_POS <- as.character(gt_matrix_co_by_marker$CHR_POS)
#
# #head(gt_matrix_co_by_marker)
#
# gt_matrix_dist <- cal_geneticMap(gt_matrix_co_by_marker = gt_matrix_co_by_marker )
#
# #head(gt_matrix_dist)
#
# ### Plotting
#
#
# gt_matrix_dist$CHR =  sapply(strsplit(gt_matrix_dist$CHR_POS,"_"),`[[`,1)
#
# gt_matrix_dist$POS =  as.numeric(sapply(strsplit(gt_matrix_dist$CHR_POS,"_"),`[[`,2))
#
# gt_matrix_dist$CHR  <-  factor(gt_matrix_dist$CHR,levels = c(seq(1:19),"X"))
#
# ## important to order the markers by chr then by position
#
# gt_matrix_dist <- gt_matrix_dist[order(gt_matrix_dist$CHR, xtfrm(gt_matrix_dist$POS)), ]
#
# gt_matrix_dist <- gt_matrix_dist %>% group_by(CHR) %>% mutate(cum_haldane = cumsum(haldane),
#                                                               cum_kosambi = cumsum(kosambi))
#
# #gt_matrix_dst[sort(gt_matrix_dst$CHR,gt_matrix_dst$POS),]
#
#
# ## total centimorgans
#
# sum(gt_matrix_dist$kosambi*100)
#
# sum(gt_matrix_dist$haldane*100)
#
# ## genetic length plot
#
#
#
# ggplot(data = gt_matrix_dist)+geom_point(mapping = aes(x = POS, y = cum_haldane ),color = "blue")+
#   facet_wrap(.~CHR)+geom_point(mapping = aes(x = POS, y = cum_kosambi ),color = "red")+ggtitle("cum_haldane (blue),cum_kosambi (red)")+
#   theme(axis.text.x = element_blank())+ylab(label = "cumulative genetic length")
#
#
# write.table(gt_matrix_dist,file = paste0("./data/",jobs,".tsv"))

