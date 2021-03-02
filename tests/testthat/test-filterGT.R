test_that("filterGT works", {
  or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  or_geno[1,] <- rep("Fail",dim(or_geno)[2])
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.)
  ft_gt <- filterGT(cr_geno)
  expect_false("1_4526088" %in% rownames(ft_gt))
  ft_gt <- filterGT(cr_geno,min_samples = 0)
  expect_true("1_4526088" %in% rownames(ft_gt))
  ft_gt <- filterGT(cr_geno,min_markers = 31)
  expect_false( sum(c("X92","X99") %in% colnames(ft_gt))>0)

  })
