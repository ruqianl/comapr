test_that("Plot counts by sample group with RSE", {

  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,
                               barcodeFile=NULL,
                               minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               minCellSNP = 1)

  s1_counts <- countCOs(s1_rse_state)
  s1_counts$sampleGroup <- "s1"
  expect_s3_class(plotCount(s1_counts,group_by = 'sampleGroup'),"ggplot")
})

test_that("Plot counts by chr and sample group with RSE", {

  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,
                               barcodeFile=NULL,
                               minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               minCellSNP = 1)

  s1_counts <- countCOs(s1_rse_state)
  s1_counts$sampleGroup <- "s1"

  expect_s3_class(plotCount(s1_counts,
                            group_by = 'sampleGroup',
                            by_chr = TRUE),
                  "ggplot")

})


test_that("Plot counts by sample group with GRanges input", {

  data(snp_geno_gr)
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  #  co_geno$interval_ID <- rownames(co_geno)
  # melt_geno <- melt(co_geno,id= 'interval_ID')
  # colnames(melt_geno) <- c("interval_ID","SampleID","Cross_over")
  marker_gr_cos <- countCOs(snp_geno_gr)

  expect_s3_class(plotCount(marker_gr_cos,group_by = 'sampleGroup'),"ggplot")
})



test_that("Plot counts by chr sample group with GRanges input", {

  data(snp_geno_gr)
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  #  co_geno$interval_ID <- rownames(co_geno)
  # melt_geno <- melt(co_geno,id= 'interval_ID')
  # colnames(melt_geno) <- c("interval_ID","SampleID","Cross_over")
  marker_gr_cos <- countCOs(snp_geno_gr)

  expect_s3_class(plotCount(marker_gr_cos,
                            group_by = 'sampleGroup',
                            by_chr = TRUE),
                  "ggplot")
})

