test_that("Read col MM works  for column 2 ", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_vi <-   readColMM(file = paste0(demo_path,"s1_chr1_vi.mtx"), which.col=2,
                       chunk=2)
  expect_true(Matrix::colSums(s1_vi)[1]==0)
})


test_that("Read col MM works for column 3 ", {
  demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
  s1_vi <-   readColMM(file = paste0(demo_path,"s1_chr1_vi.mtx"), which.col=3,
                       chunk=2)
  expect_true(Matrix::colSums(s1_vi)[1]==0 & Matrix::colSums(s1_vi)[2]==0 )

})
