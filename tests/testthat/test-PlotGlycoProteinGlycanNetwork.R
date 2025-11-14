test_that("PlotGlycoProteinGlycanNetwork", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoProteinGlycanNetwork(mydata, silent = TRUE)
  expect_s3_class(p, "recordedplot")
})
