test_that("combineHapState works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                               path=demo_path,barcodeFile=NULL,minSNP = 0,
                               minlogllRatio = 50,
                               bpDist = 100,maxRawCO=10,
                               minCellSNP = 0)
  
  s2_rse_state <- readHapState("s2",chroms=c("chr1","chr2"),
                               path=demo_path,
                               barcodeFile=paste0(demo_path,"s2_barcodes.txt"),minSNP = 0, 
                               minlogllRatio = 50,
                               bpDist = 100,maxRawCO=10,
                               minCellSNP = 0)
  sb <- combineHapState(s1_rse_state,s2_rse_state)
  expect_true(dim(sb)[2] == (dim(s1_rse_state) + dim(s2_rse_state))[2])
  
  expect_true(sum(colnames(sb) == c(colnames(s1_rse_state),
                                    colnames(s2_rse_state)))==dim(sb)[2] )
  
})
