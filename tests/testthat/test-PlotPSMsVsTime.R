test_that("PlotPSMsVsTime: defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotPSMsVsTime(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMsVsTime: gradientLength", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotPSMsVsTime(mydata, gradientLength = 20, silent = TRUE)
  p2 <- PlotPSMsVsTime(mydata, gradientLength = 2000, silent = TRUE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})

test_that("PlotPSMsVsTime: type", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotPSMsVsTime(mydata, type = "allGlyco", silent = TRUE)
  p2 <- PlotPSMsVsTime(mydata, type = "both", silent = TRUE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})

test_that("PlotPSMsVsTime: bindwidth", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotPSMsVsTime(mydata, binWidth = 12, silent = TRUE)
  expect_s3_class(p1, "gg")
})
