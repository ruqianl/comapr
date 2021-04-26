BiocParallel::register(BPPARAM =BiocParallel::SerialParam())

test_that("plotGeneticDist works", {

 # data(snp_geno_gr)
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)

  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  marker_gr_cos <- countCOs(snp_geno_gr)

  expect_equal(sum(marker_gr_cos[,"X92"]$X92),3)

  dist <- calGeneticDist(marker_gr_cos)
  expect_true(sum(dist$kosambi_cm)-169.53<0.01)
  #expect_true(abs(dist$kosambi_cm[1]-13.67603)<1e-3)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="1"]) - 47.71165 <0.001)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="2"]) - 51.27605 <0.001)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="4"]) - 42.82993 <0.001)
  error_plot <- plotGeneticDist(dist,bin = FALSE,chr ="chr1")
  expect_error(plot(error_plot))

  good_plot <- plotGeneticDist(dist,bin = FALSE,chr ="1")
  expect_s3_class(plot(good_plot),"ggplot")

  })

test_that("plotting genetic distances for bined intervals", {
  data(snp_geno_gr)
  corrected_geno <- correctGT(gt_matrix = mcols(snp_geno_gr),
                              ref = parents_geno$ref,
                              alt = parents_geno$alt)
  GenomicRanges::mcols(snp_geno_gr) <- corrected_geno
  marker_gr_cos <- countCOs(snp_geno_gr)

  dist <- calGeneticDist(marker_gr_cos,bin_size = 1e7)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="1"]) - 47.71165 <0.001)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="2"]) - 51.27605 <0.001)
  expect_true(sum(dist$kosambi_cm[as.character(seqnames(dist))=="4"]) - 42.82993 <0.001)

  good_plot <- plotGeneticDist(dist,bin = FALSE,chr ="1")
  expect_s3_class(plot(good_plot),"ggplot")
})



test_that("plot whole genome genetic distances for RangedSummarizedExperiment", {
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
            dist <- calGeneticDist(s1_counts,bin_size = 1e7)

            expect_true(sum(dist$X)-41.84747 <0.001)
            expect_s3_class(plotWholeGenome(dist),"ggplot")
    })

test_that("plot genetic distances works for RangedSummarizedExperiment ", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",
                               chroms=c("chr1"),
                               path=demo_path,
                               barcodeFile=NULL,
                               minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               minCellSNP = 1)

  s1_counts <- countCOs(s1_rse_state)
  dist <- calGeneticDist(s1_counts)

  expect_equivalent(colSums(as.matrix(assay(dist))), c(1,1,0,0,0 ))
  expect_true(sum(rowRanges(dist)$kosambi)-41.84747 <0.001)
  expect_s3_class(plotWholeGenome(rowRanges(dist)),"ggplot")

})

