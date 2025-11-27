test_that("PlotSiteQuantification: default quantification", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]

  result <- tryCatch(
    PlotSiteQuantification(mydata, protein = testProtein,
                           site = testSite, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(result$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotSiteQuantification: intensity thresholds", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]
  testSite <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & UniprotIDs == testProtein)$ModificationID[1]

  result1 <- tryCatch(
    PlotSiteQuantification(mydata, protein = testProtein, cutoff = 0.001,
                           site = testSite, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result1, "error")) {
    expect_match(result1$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result1, "gg")
  }

  result2 <- tryCatch(
    PlotSiteQuantification(mydata, protein = testProtein, cutoff = "0.001%",
                           site = testSite, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result2, "error")) {
    expect_match(result2$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result2, "gg")
  }
})

test_that("PlotSiteQuantification: catch error", {
  skip_if_not(exists("mydata"), "User data not loaded")
  result <- tryCatch(
    PlotSiteQuantification(mydata, protein = "testProtein", cutoff = 0.001,
                           site = "testSite", silent = TRUE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(result$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result, "gg")
  }
})
