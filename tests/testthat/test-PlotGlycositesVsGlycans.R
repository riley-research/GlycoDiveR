test_that("PlotGlycositesVsGlycans: default", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycositesVsGlycans(mydata, silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycositesVsGlycans: negative cutoffs", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycositesVsGlycans(mydata, labelGlycositeCutoff  = -10,
                               labelGlycanCutoff = -10, silent = TRUE)
  expect_s3_class(p, "gg")
})
