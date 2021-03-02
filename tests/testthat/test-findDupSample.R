test_that("findDupSamples works", {
  or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
  rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
  or_geno[,1] <- or_geno[,5]
  cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
                       alt = snp_geno$FVB.NJ..i.)
  dups <- findDupSamples(cr_geno,)

  expect_equal(dups[,1], c("X92","X96"))


})

# 
# test_that("findDupMarkers works", {
#   or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#   rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#   or_geno[,1] <- or_geno[,5]
#   cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#                        alt = snp_geno$FVB.NJ..i.,
#                        chr = snp_geno$CHR)
#   dups <- findDupMarkers(cr_geno,plot = TRUE)
#   expect_true(length(dups)==2)
# 
# 
# })
