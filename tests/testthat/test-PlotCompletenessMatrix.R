test_that("PlotCompletenessMatrix: peptideType = glyco", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCompletenessMatrix(mydata, peptideType = "glyco", silent = TRUE)
  expect_true(is(p, "Heatmap") || is(p, "HeatmapList"))
})

test_that("PlotCompletenessMatrix: peptideType = other", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCompletenessMatrix(mydata, peptideType = "other", silent = TRUE)
  expect_true(is(p, "Heatmap") || is(p, "HeatmapList"))
})
