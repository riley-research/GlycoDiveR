test_that("PlotSiteGlycanComposition: glycoprotein quantification", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotSiteGlycanComposition(mydata, protein = testProtein, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotSiteGlycanComposition: horizontalPoints", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotSiteGlycanComposition(mydata, protein = testProtein, horizontalPoints = 10000, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotSiteGlycanComposition: yCorrection", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotSiteGlycanComposition(mydata, protein = testProtein, yCorrection = -1, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotSiteGlycanComposition: yNudge", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotSiteGlycanComposition(mydata, protein = testProtein, yNudge = 0.000001, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotSiteGlycanComposition: boxSpacing", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  p1 <- PlotSiteGlycanComposition(mydata, protein = testProtein, boxSpacing = 0.2, silent = TRUE)
  expect_s3_class(p1, "gg")
})
