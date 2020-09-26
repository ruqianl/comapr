test_that("Calculate genetic distances", {
  data("snp_geno_gr")
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
#  co_geno$interval_ID <- rownames(co_geno)
  # melt_geno <- melt(co_geno,id= 'interval_ID')
  # colnames(melt_geno) <- c("interval_ID","SampleID","Cross_over")
  marker_gr_cos <- countCOs(snp_geno_gr)
  
  expect_equal(sum(marker_gr_cos[4:5,"X92"]$X92),1)
  dist <- calGeneticDist(marker_gr_cos)
  expect_true(abs(dist$kosambi_cm[4]-13.67603)<1e-3)

})


# test_that("Calculate genetic distances successfully for all failed marker", {
#   or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#   rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#   or_geno[1,] <- rep("Fail",dim(or_geno)[2])
#   cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#                        alt = snp_geno$FVB.NJ..i.,
#                        chr = snp_geno$CHR)
#   co_geno <- detectCO(cr_geno,
#                       chrs = snp_geno$CHR,
#                       chrPos = snp_geno$POS)
#   rb_rate <- calGeneticDist(co_geno)
#   expect_true(is.numeric(rb_rate["1_21638464_22665060",]$na_rate))
#   expect_equivalent(round(rb_rate["1_21638464_22665060",]$na_rate,3),0.182)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$t_counts,1)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$total_na,4)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$f_counts,17)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$total_calls,18)
# 
# })
# 
# 
# test_that("Calculate genetic distances successfully for rb>0.5 marker", {
#   or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
#   rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#   cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#                        alt = snp_geno$FVB.NJ..i.,
#                        chr = snp_geno$CHR)
#   co_geno <- detectCO(cr_geno,
#                       chrs = snp_geno$CHR,
#                       chrPos = snp_geno$POS)
#   co_geno[3,] <- rep(TRUE,dim(co_geno)[2])
# 
# 
#   rb_rate <- calGeneticDist(co_geno)
#   expect_true(is.numeric(rb_rate["1_21638464_22665060",]$na_rate))
#   expect_equivalent(round(rb_rate["1_21638464_22665060",]$na_rate,3),0.182)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$t_counts,1)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$total_na,4)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$f_counts,17)
#   expect_equivalent(rb_rate["1_21638464_22665060",]$total_calls,18)
#   expect_false("1_5595513_6057774" %in% rownames(rb_rate))
# })
