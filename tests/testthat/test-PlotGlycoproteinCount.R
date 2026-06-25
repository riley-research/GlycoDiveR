test_that("PlotGlycoproteinCount technicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoproteinCount(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycoproteinCount biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoproteinCount(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycoproteinCount condition biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoproteinCount(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
