test_that("countCOs works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,barcodeFile=NULL,minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,maxRawCO=10)
  count_rse <- countCOs(s1_rse_state)
  ## colsums for each cell is either 0 or 1
  expect_true(sum(colSums( as.matrix(assay(  count_rse))) !=0 & 
                    colSums( as.matrix(assay(  count_rse)))!=1)==0)
})
