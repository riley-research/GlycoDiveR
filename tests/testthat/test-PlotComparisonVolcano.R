test_that("PlotComparisonVolcano: default", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison)
  expect_s3_class(p1, "gg")
})

test_that("PlotComparisonVolcano: statistic pvalue", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison,
                              statistic = "pvalue")
  expect_s3_class(p1, "gg")
})

test_that("PlotComparisonVolcano: statisticalCutoff 0.1", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison,
                              statisticalCutoff = 0.1)
  expect_s3_class(p1, "gg")
})

test_that("PlotComparisonVolcano: whichLabel none", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison,
                              whichLabel = "none")
  expect_s3_class(p1, "gg")
})

test_that("PlotComparisonVolcano: testLabel", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  testLabel <- paste(comparison$Proteins, comparison$ModificationID, sep ="-")[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison,
                              whichLabel = testLabel)
  expect_s3_class(p1, "gg")
})

test_that("PlotComparisonVolcano: maxOverlaps NA", {
  skip_if_not(exists("comparison"), "User data not loaded")
  testComparison <- unique(comparison$Label)[1]
  p1 <- PlotComparisonVolcano(comparison, whichComparison = testComparison,
                              maxOverlaps = NA)
  expect_s3_class(p1, "gg")
})
