test_that("get Cell AF datatracks works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,barcodeFile=NULL,minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,maxRawCO=10,
                               minCellSNP = 0)
  s1_counts <- countCOs(s1_rse_state)

  af_co_tracks <- getCellAFTrack(chrom ="chr1",
                                 path_loc = demo_path,
                                 sampleName = "s1",
                                 barcodeFile = paste0(demo_path,
                                                      "s1_barcodes.txt"),
                 cellBarcode = "BC1",
                 co_count = s1_counts)
  expect_true(length(af_co_tracks) == 2)
  expect_true(GenomicRanges::start(af_co_tracks$co_range)==3201)
  expect_true(GenomicRanges::end(af_co_tracks$co_range)==5499)

  })


test_that("get cell CO range works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,barcodeFile=NULL,minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,maxRawCO=10,
                               minCellSNP = 0)
  s1_counts <- countCOs(s1_rse_state)

  co_ranges <- getCellCORange(cellBarcode = "BC1",
                                 co_count = s1_counts)
  expect_true(is(co_ranges,"GRanges"))
  expect_true(GenomicRanges::start(co_ranges)==3201)

  expect_true(GenomicRanges::seqnames(co_ranges)=="chr1")
  expect_true(GenomicRanges::end(co_ranges)=="5499")


})
