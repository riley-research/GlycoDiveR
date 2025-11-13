test_that("PlotPSMCount works", {
  p <- PlotPSMCount(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMCount works", {
  p <- PlotPSMCount(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})
