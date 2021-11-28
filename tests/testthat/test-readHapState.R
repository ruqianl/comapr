test_that("readHapState from raw files works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  sampleName <- "s1"
  s1_rse_state <- readHapState(sampleName,
                               chroms=c("chr1"),
                               path=demo_path,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               barcodeFile=NULL,
                               minSNP = 0,
                               minCellSNP = 1)

  expect_true(length(GenomeInfoDb::seqlevels(rowRanges(s1_rse_state))) ==1)

  s1_rse_state <- readHapState(sampleName,chroms=c("chr1"),
                               path=demo_path,
                               barcodeFile=paste0(demo_path,"s1_barcodes.txt"),
                               minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               minCellSNP = 1)
  expect_true(length(GenomeInfoDb::seqlevels(rowRanges(s1_rse_state))) ==1)

  s1_rse_state <- readHapState(sampleName,chroms=c("chr1","chr2"),
                               path=demo_path,
                               barcodeFile=paste0(demo_path,"s1_barcodes.txt"),
                               minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,
                               maxRawCO=10,
                               minCellSNP = 1)
  expect_true(length(GenomeInfoDb::seqlevels(rowRanges(s1_rse_state))) ==2)



})


test_that("readHapState from raw files works for multiple filters", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  sampleName <- "s1"

  s1_rse_state <- readHapState(sampleName,chroms=c("chr1","chr2"),
                               path=demo_path,
                               barcodeFile=paste0(demo_path,"s1_barcodes.txt"),
                               minSNP = c(0,0),
                               minlogllRatio = c(50,49),
                               bpDist = c(100,99.1),
                               maxRawCO=10,
                               minCellSNP = 1)
  expect_true(length(GenomeInfoDb::seqlevels(rowRanges(s1_rse_state))) ==2)



})
