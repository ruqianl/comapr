test_that("label gt error", {
  expect_error(label_gt(s_gt = snp_geno[,grep("X",colnames(snp_geno))],
                        ref = snp_geno$C57BL.6J,
                        alt = snp_geno$FVB.NJ..i.))
})

test_that("label gt successfully coverts Fail", {
  converted_gt <- label_gt(s_gt = snp_geno[,"X98"],
                        ref = snp_geno$C57BL.6J,
                        alt = snp_geno$FVB.NJ..i.)
  expect_equal(converted_gt[5],"Fail")
})

test_that("label gt successfully coverts unknow GT to missinng", {
  s_gt <- snp_geno[,"X98"]
  s_gt[2] <- "KK"
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)

  expect_equal(converted_gt[2],"missing")
})


test_that("label gt successfully coverts homo_ref", {
  s_gt <- snp_geno$C57BL.6J
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)

  expect_equal(sum(converted_gt %in% "Homo_ref"),length(s_gt))
})


test_that("label gt successfully coverts homo_alt", {
  s_gt <- snp_geno$FVB.NJ..i.
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)

  expect_equal(sum(converted_gt %in% "Homo_alt"),length(s_gt))
})



test_that("correct_ref successfully corrects Homo_ref to Fail", {
  s_gt <- snp_geno[,5]
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  converted_gt[4] <- "Homo_ref"

  expect_equal(correct_ref(converted_gt,change_to = "Fail")[4],"Fail")
})



test_that("correct_ref successfully corrects Homo_ref to Het", {
  s_gt <- snp_geno[,5]
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  converted_gt[4] <- "Homo_ref"

  expect_equal(correct_ref(converted_gt,change_to = "Het")[4],"Het")
})



test_that("infer Failed genotype successfully", {
  s_gt <- snp_geno[,5]
  s_gt[10] <- "Fail"
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  infer_fail <- fill_fail(converted_gt)
  expect_equal(correct_ref(converted_gt,change_to = "Het")[4],"Het")
})


