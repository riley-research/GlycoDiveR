test_that("PlotErrorVersusRt defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotErrorVersusRt(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotErrorVersusRt zoom", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotErrorVersusRt(mydata, lowerLimit = -10, upperLimit = 10, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotErrorVersusRt glyco", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotErrorVersusRt(mydata, type = "glyco", silent = TRUE)
  expect_s3_class(p, "gg")
})
