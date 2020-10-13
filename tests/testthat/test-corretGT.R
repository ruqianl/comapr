test_that("correctGT works", {
  or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.,
                       chr = snp_geno$CHR)
  expect_false(any(cr_geno=="Fail",na.rm = TRUE))
})
