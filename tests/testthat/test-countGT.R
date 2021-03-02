test_that("countGT works", {
  or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  or_geno[1,] <- rep("Fail",dim(or_geno)[2])
  cr_geno <- correctGT(or_geno,
                       ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.)
  p1 <- countGT(cr_geno,plot = FALSE)
  p2 <- countGT(cr_geno,plot = TRUE, interactive = FALSE)
  p3 <- countGT(cr_geno,plot = TRUE, interactive = TRUE)
  p4 <- countGT(cr_geno)

  expect_true(class(p1) == 'list')
  expect_true(class(p2) == 'list')
  expect_true(class(p3$ply)[1] == 'plotly')
  expect_equal(p1,p4[c(2,3)])
})
