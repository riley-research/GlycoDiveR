test_that("PlotGlycanCompositionPie", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycansPerSite(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})
