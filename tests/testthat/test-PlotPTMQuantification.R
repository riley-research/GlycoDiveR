test_that("PlotPTMQuantification: defaults", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]

  result <- tryCatch(
    PlotPTMQuantification(mydata, protein = testProtein, rowFontSize = 10, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(result$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotPTMQuantification: linewidth", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]

  result1 <- tryCatch(
    PlotPTMQuantification(mydata, protein = testProtein, lineWidth = 10, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result1, "error")) {
    expect_match(result1$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result1, "gg")
  }

  result2 <- tryCatch(
    PlotPTMQuantification(mydata, protein = testProtein, lineWidth = NA, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result2, "error")) {
    expect_match(result2$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result2, "gg")
  }
})

test_that("PlotPTMQuantification: rowFontSize", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]

  result <- tryCatch(
    PlotPTMQuantification(mydata, protein = testProtein, rowFontSize = 10, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(result$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result, "gg")
  }
})

test_that("PlotPTMQuantification: showRowNames", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testProtein <- subset(mydata$PTMTable, GlycanType != "NonGlyco" & GlycanQValue < 0.01)$UniprotIDs[1]

  result <- tryCatch(
    PlotPTMQuantification(mydata, protein = testProtein, showRowNames = FALSE, silent = TRUE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    expect_match(result$message, "No quantitative data found. Aborting.")
  } else {
    expect_s3_class(result, "gg")
  }
})
