test_that("label gt error", {
  data(snp_geno)
  expect_error(.label_gt(s_gt = snp_geno[,grep("X",colnames(snp_geno))],
                        ref = snp_geno$C57BL.6J,
                        alt = snp_geno$FVB.NJ..i.))
})

test_that("label gt successfully coverts Fail", {
  data(snp_geno)
  converted_gt <- .label_gt(s_gt = snp_geno[,"X98"],
                        ref = snp_geno$C57BL.6J,
                        alt = snp_geno$FVB.NJ..i.)
  expect_equal(converted_gt[5],"Fail")
})

test_that("label gt successfully coverts unknow GT to missinng", {
  data(snp_geno)
  s_gt <- snp_geno[,"X98"]
  s_gt[2] <- "KK"
  converted_gt <- .label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.,
                           failed = "Fail")

  expect_equal(converted_gt[2],"Fail")
})


test_that("label gt successfully coverts homo_ref", {
  data(snp_geno)
  s_gt <- snp_geno$C57BL.6J
  converted_gt <- .label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)

  expect_equal(sum(converted_gt %in% "Homo_ref"),length(s_gt))
})


test_that("label gt successfully coverts homo_alt", {
  data(snp_geno)
  s_gt <- snp_geno$FVB.NJ..i.
  converted_gt <- .label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)

  expect_equal(sum(converted_gt %in% "Homo_alt"),length(s_gt))
})



test_that("correct_ref successfully corrects Homo_ref to NA", {
  data(snp_geno)
  s_gt <- snp_geno[,5]
  s_gt[4] <- "Homo_ref"
  converted_gt <- correctGT(gt_matrix = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)


  expect_true(is.na(converted_gt[4]))
})

