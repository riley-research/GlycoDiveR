test_that("PlotPCA defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) > 1, "Only one group")
  p <- PlotPCA(mydata, quantType = "normalized", dropNoQuant = TRUE, type = "glyco",
               label = TRUE, thresholdMode = "total", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPCA defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) > 1, "Only one group")
  p <- PlotPCA(mydata, quantType = "nonNormalized", dropNoQuant = FALSE, type = "all",
               label = FALSE, thresholdMode = "group", minPeptideCoverage = c(0.5, 1),
               silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPCA defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) > 1, "Only one group")
  p <- PlotPCA(mydata, returnType = "loadingsPlot", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPCA defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) > 1, "Only one group")
  p <- PlotPCA(mydata, returnType = "loadingsData", silent = TRUE)
  expect_s3_class(p, "data.frame")
})
