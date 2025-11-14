test_that("PlotPSMCount technicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoPSMCount(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMCount biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoPSMCount(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSMCount condition biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoPSMCount(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
