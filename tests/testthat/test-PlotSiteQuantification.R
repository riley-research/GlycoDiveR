test_that("PlotSiteQuantification: default quantification", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]

  result <- PlotSiteQuantification(mydata, whichProtein = testProtein,
                           site = testSite, silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotSiteQuantification: intensity thresholds", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]

  result1 <- PlotSiteQuantification(mydata, whichProtein = testProtein, cutoff = 0.001,
                           site = testSite, silent = TRUE)

  if ((is.character(result1) && grepl("No data|No quantitative data", result1)) |
      is.null(result1)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result1, "gg")
  }

  result2 <- PlotSiteQuantification(mydata, whichProtein = testProtein, cutoff = "0.001%",
                           site = testSite, silent = TRUE)

  if ((is.character(result2) && grepl("No data|No quantitative data", result2)) |
      is.null(result2)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result2, "gg")
  }
})

test_that("PlotSiteQuantification: catch error", {
  skip_if_not(exists("mydata"), "User data not loaded")
  result <- PlotSiteQuantification(mydata, whichProtein = "testProtein", cutoff = 0.001,
                           site = "testSite", silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)){
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})
