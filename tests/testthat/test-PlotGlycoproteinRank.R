test_that("PlotGlycoproteinRank all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoproteinRank(mydata, bins = 20, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycoproteinRank all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycoproteinRank(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
