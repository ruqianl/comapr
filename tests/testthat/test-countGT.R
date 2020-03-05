test_that("countGT works", {
  or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  or_geno[1,] <- rep("Fail",dim(or_geno)[2])
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.,
                       chr = snp_geno$CHR)
  p1 <- countGT(cr_geno,plot = FALSE)
  p2 <- countGT(cr_geno,plot = TRUE, interactive = FALSE)
  p3 <- countGT(cr_geno,plot = TRUE, interactive = TRUE)
  p4 <- countGT(cr_geno,plot_wg = "TRUE")

  expect_true(class(p1)[1] == 'gg')
  expect_true(class(p2)[2] == 'ggplot')
  expect_true(class(p3)[2] == 'ggplot')
  expect_true(class(p4)[1] == 'gg')
})
