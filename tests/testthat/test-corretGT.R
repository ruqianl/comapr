test_that("correctGT works", {
  data(snp_geno)
  data("parents_geno")
  or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.)
  expect_false(any(cr_geno=="Fail",na.rm = TRUE))
  expect_true(is.na(cr_geno["1_21638464","X92"]))
})
