test_that("Calculate cie ", {
  data(snp_geno_gr)
  n_interval <- 5
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  #  co_geno$interval_ID <- rownames(co_geno)
  # melt_geno <- melt(co_geno,id= 'interval_ID')
  # colnames(melt_geno) <- c("interval_ID","SampleID","Cross_over")
  marker_gr_cos <- countCOs(snp_geno_gr)
  
  coc_df <- calCoC(marker_gr_cos,n_intervals = n_interval)
  expect_equal(length(unique(c(coc_df$X1,coc_df$X2))), n_interval)
  
  coc_df2 <- calCoC(marker_gr_cos,n_intervals = n_interval,
                    group_prefix = c("X1","X9"))
  expect_equal(length(unique(c(coc_df2$X1,coc_df2$X2))), n_interval)
  
  expect_equal(unique(coc_df2$group_prefix), c("X1","X9"))
})
# 
# coc_df %>% group_by(inter_dist,chr) %>% mutate(mean_coc = mean(na.omit(coc))) %>%
#   
#   ggplot()+geom_point(mapping = aes(x=inter_dist,
#                                     y= mean_coc))+
#   geom_smooth(mapping = aes(x=inter_dist,
#                             y= coc))+theme_classic()+facet_wrap(.~chr)
# 
# coc_df2 %>% group_by(inter_dist,chr) %>% mutate(mean_coc = mean(na.omit(coc))) %>%
#   
#   ggplot()+geom_point(mapping = aes(x=inter_dist,
#                                     y= mean_coc,color=group_prefix))+
#   geom_smooth(mapping = aes(x=inter_dist,
#                             y= coc,
#                             color=group_prefix))+theme_classic()+facet_wrap(.~chr)
