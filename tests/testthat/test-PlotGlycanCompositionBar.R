test_that("PlotGlycanCompositionBar grouping condition", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionBar(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycanCompositionBar grouping technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionBar(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})
