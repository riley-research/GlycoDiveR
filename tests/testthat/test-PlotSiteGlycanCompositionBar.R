test_that("PlotSiteGlycanCompositionBar grouping = none", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotSiteGlycanCompositionBar(mydata, whichProtein = testProtein,
                                   grouping = "none", intensity = "raw", silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotSiteGlycanCompositionBar grouping = technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotSiteGlycanCompositionBar(mydata, whichProtein = testProtein,
                                         grouping = "technicalReps", intensity = "log2", silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotSiteGlycanCompositionBar grouping = biologicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotSiteGlycanCompositionBar(mydata, whichProtein = testProtein,
                                         grouping = "biologicalReps", intensity = "log10", silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotSiteGlycanCompositionBar grouping = condition", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotSiteGlycanCompositionBar(mydata, whichProtein = testProtein,
                                         grouping = "condition", intensity = "log10",
                                         scales = "stack", silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)
  } else {
    expect_s3_class(result, "gg")
  }
})
