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
  s_gt[9] <- "Fail"

  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  infer_fail <- fill_fail(converted_gt)
  expect_true(infer_fail[9]=="Homo_alt")
})



test_that("infer Failed genotype successfully to Fail", {
  s_gt <- snp_geno[,5]
  s_gt[length(s_gt)] <- "Fail"
## the 10th marker is the end of chromosome 1.
## thus we cannot infer the genotype of it

  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  infer_fail <- fill_fail(converted_gt)
  expect_false(infer_fail[length(s_gt)]=="Homo_alt")
  expect_true(infer_fail[length(s_gt)]=="Fail")
})



test_that("infer Failed genotype successfully", {
  s_gt <- snp_geno[,5]
  s_gt[10] <- "Fail"
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
  infer_fail <- fill_fail(converted_gt)
  expect_equal(infer_fail[4],"Homo_alt")
})



test_that("infer Failed genotype successfully by chrs", {
  s_gt <- snp_geno[,5]
  s_gt[c(10,11,23)] <- "Fail"
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.,)
  infer_fail <- fill_fail(converted_gt,chr = as.factor(snp_geno$CHR))
  expect_equal(infer_fail[10],"Fail")
  expect_equal(infer_fail[11],"Fail")
  expect_equal(infer_fail[23],"Fail")

})


test_that("change missing successfully change Fail to missing", {
  s_gt <- snp_geno[,5]
  s_gt[c(10,11,23)] <- "Fail"
  converted_gt <- label_gt(s_gt = s_gt,
                           ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.,)
  infer_fail <- fill_fail(converted_gt,chr = as.factor(snp_geno$CHR))
  change_miss <- change_missing(infer_fail)
  expect_true(is.na(change_miss[10]))
  expect_true(is.na(change_miss[11]))
  expect_true(is.na(change_miss[23]))

})

