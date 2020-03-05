test_that("detect crossovers successfully", {
  or_geno <-snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.,
                       chr = snp_geno$CHR)
  co_geno <- detectCO(cr_geno,
                     chrs = snp_geno$CHR,
                     chrPos = snp_geno$POS)
  count_geno <- detectCO(cr_geno,
                         chrs = snp_geno$CHR,
                         chrPos = snp_geno$POS,
                         type="counts")

  expect_true(co_geno[16,1])
  expect_equal(count_geno[16,1],1)

  })
