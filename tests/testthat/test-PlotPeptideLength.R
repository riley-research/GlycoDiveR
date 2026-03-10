test_that("PlotPeptideLength works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotPeptideLength(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})
