test_that("diagnostic function for per cell chr works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  sampleName <- "s1"
  pcqc <- perCellChrQC(sampleName,
                               chroms=c("chr1"),
                               path=demo_path)
  expect_true(is.list(pcqc))
  expect_true(length(pcqc)==2)
  expect_s3_class(pcqc$plot,"ggplot")

})

test_that("diagnostic per seg works", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  sampleName <- "s1"
  psqc <- perSegChrQC(sampleName,
                       chroms=c("chr1"),
                       path=demo_path)

  expect_s3_class(psqc,"gtable")

})
