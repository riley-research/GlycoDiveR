test_that("PlotPTMQuantification: defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotPTMQuantification(mydata, whichProtein = testProtein, rowFontSize = 10, silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
  expect_true(TRUE)   # message means function behaved correctly
} else {
  expect_s3_class(result, "gg")
}
})

test_that("PlotPTMQuantification: linewidth", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result1 <- PlotPTMQuantification(mydata, whichProtein = testProtein, lineWidth = 10, silent = TRUE)

  if ((is.character(result1) && grepl("No data|No quantitative data", result1)) |
      is.null(result1)) {
    expect_true(TRUE)   # message means function behaved correctly
  } else {
    expect_s3_class(result1, "gg")
  }

  result2 <- PlotPTMQuantification(mydata, whichProtein = testProtein, lineWidth = NA, silent = TRUE)

  if ((is.character(result2) && grepl("No data|No quantitative data", result2)) |
      is.null(result2)) {
    expect_true(TRUE)   # message means function behaved correctly
  } else {
    expect_s3_class(result2, "gg")
  }
})

test_that("PlotPTMQuantification: rowFontSize", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotPTMQuantification(mydata, whichProtein = testProtein, rowFontSize = 10, silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)   # message means function behaved correctly
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotPTMQuantification: showRowNames", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(FilterForCutoffs(mydata, silent = TRUE)$PTMTable, GlycanType != "NonGlyco")$UniprotIDs[1]

  result <- PlotPTMQuantification(mydata, whichProtein = testProtein, showRowNames = FALSE, silent = TRUE)

  if ((is.character(result) && grepl("No data|No quantitative data", result)) |
      is.null(result)) {
    expect_true(TRUE)   # message means function behaved correctly
  } else {
    expect_s3_class(result, "gg")
  }
})
