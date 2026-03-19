test_that("PlotVennDiagram: double", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) >= 2, "User data not loaded")
  groups <- unique(mydata$PSMTable$Condition)
  p1 <- PlotVennDiagram(mydata, groups = groups[1:2], silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotVennDiagram: double", {
  skip_if_not(exists("mydata"), "User data not loaded")
  skip_if_not(length(unique(mydata$PSMTable$Condition)) >= 3, "User data not loaded")
  groups <- unique(mydata$PSMTable$Condition)
  p1 <- PlotVennDiagram(mydata, groups = groups[1:3], level = "glycan", silent = TRUE)
  expect_s3_class(p1, "gg")
})
