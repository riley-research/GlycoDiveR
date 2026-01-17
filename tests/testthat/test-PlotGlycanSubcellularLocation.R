test_that("PlotGlycanSubcellularLocation count", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanSubcellularLocation(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycanSubcellularLocation intensity", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanSubcellularLocation(mydata, summaryFunction = "intensity", silent = TRUE)
  expect_s3_class(p, "gg")
})
