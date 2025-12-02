test_that("PlotGlycositeIdsPerRun", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]
  p <- PlotGlycositeIdsPerRun(mydata, whichProtein = testProtein, silent = TRUE)
  expect_s3_class(p, "gg")
})
