test_that("PlotCV glyco plot", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCV(mydata, type = "glyco", returnType = "plot", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotCV  all data", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCV(mydata, type = "all", returnType = "data", returnThreshold = 1000, silent = TRUE)
  expect_s3_class(p, "data.frame")
})

test_that("PlotCV  all data GD", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCV(mydata, type = "all", returnType = "GDdata", returnThreshold = 1000, silent = TRUE)
  expect_true(is.list(p))
})
