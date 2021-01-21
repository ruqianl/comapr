# test boostrap

test_that("bootstrapDist works", {
  data(snp_geno_gr)
  
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  #  co_geno$interval_ID <- rownames(co_geno)
  # melt_geno <- melt(co_geno,id= 'interval_ID')
  # colnames(melt_geno) <- c("interval_ID","SampleID","Cross_over")
  marker_gr_cos <- countCOs(snp_geno_gr)
  
  # expect_equal(sum(marker_gr_cos[4:5,"X92"]$X92),1)
  # dist <- calGeneticDist(marker_gr_cos)
  # expect_true(abs(dist$kosambi_cm[4]-13.67603)<1e-3)
  # 
  bootstrap_results <- bootstrapDist(marker_gr_cos,B=10,
                                     by_group  = c("X1","X9"))
  expect_true(is.matrix(bootstrap_results))
})
