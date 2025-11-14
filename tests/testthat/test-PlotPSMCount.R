test_that("PlotPSMCount: technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotPSMCount(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMCount: biologicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotPSMCount(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMCount: condition biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotPSMCount(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
