test_that("PlotPTMQuantification: defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotPTMQuantification(mydata, protein = testProtein, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotPTMQuantification: linewidth", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotPTMQuantification(mydata, protein = testProtein, lineWidth = 10, silent = TRUE)
  p2 <- PlotPTMQuantification(mydata, protein = testProtein, lineWidth = NA, silent = TRUE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})

test_that("PlotPTMQuantification: rowFontSize", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotPTMQuantification(mydata, protein = testProtein, rowFontSize = 10, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotPTMQuantification: showRowNames", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotPTMQuantification(mydata, protein = testProtein, showRowNames = FALSE, silent = TRUE)
  expect_s3_class(p1, "gg")
})
