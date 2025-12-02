test_that("PlotCompletenessHeatmap: peptideType = glyco", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCompletenessHeatmap(mydata, peptideType = "glyco", silent = TRUE)
  expect_true(is(p, "Heatmap") || is(p, "HeatmapList"))
})

test_that("PlotCompletenessHeatmap: peptideType = other", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotCompletenessHeatmap(mydata, peptideType = "other", silent = TRUE)
  expect_true(is(p, "Heatmap") || is(p, "HeatmapList"))
})
