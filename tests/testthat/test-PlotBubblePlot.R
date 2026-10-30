test_that("multiplication works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotBubblePlot(mydata)
  expect_s3_class(p, "gg")
})
