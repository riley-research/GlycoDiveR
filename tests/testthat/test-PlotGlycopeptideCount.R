test_that("PlotGlycopeptideCount technicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycopeptideCount(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycopeptideCount biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycopeptideCount(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycopeptideCount condition biologicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycopeptideCount(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
