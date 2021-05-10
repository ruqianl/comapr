test_that("get Mean DP track works works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  meanDP_track <- getMeanDPTrack(chrom ="chr1",
                                 path_loc = demo_path,
                                 barcodeFile = paste0(demo_path,"s1_barcodes.txt"),
                                 sampleName = "s1")
  expect_silent(Gviz::plotTracks(meanDP_track))
})


test_that("get celll DP track works works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  bc1_dp_track <- getCellDPTrack(chrom ="chr1",
                                 path_loc = demo_path,
                                 barcodeFile = paste0(demo_path,"s1_barcodes.txt"),
                                 sampleName = "s1",
                                 chunk = 2,
                                 cellBarcode = "BC1",
                                 log = TRUE)
  expect_silent(Gviz::plotTracks(bc1_dp_track))
})
