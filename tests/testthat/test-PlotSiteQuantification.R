test_that("PlotSiteQuantification: default quantification", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]
  p1 <- PlotSiteQuantification(mydata, protein = testProtein,
                               site = testSite, silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotSiteQuantification: intensity thresholds", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]
  p1 <- PlotSiteQuantification(mydata, protein = testProtein, cutoff = 0.001,
                               site = testSite, silent = TRUE)
  p2 <- PlotSiteQuantification(mydata, protein = testProtein, cutoff = "0.001%",
                               site = testSite, silent = TRUE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})

test_that("PlotSiteQuantification: catch error", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotSiteQuantification(mydata, protein = "testProtein", cutoff = 0.001,
                               site = "testSite", silent = TRUE)
  expect_null(p1)
})
