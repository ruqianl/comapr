
test_that("plotGTFreq works",{
  or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  or_geno[1,] <- rep("Fail",dim(or_geno)[2])
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.,
                       chr = snp_geno$CHR)
  ft_gt <- filterGT(cr_geno)
  p <- plotGTFreq(ft_gt)
  expect_true(class(p)[2]=="ggplot")
  
})
