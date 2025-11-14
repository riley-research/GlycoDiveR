test_that("PlotQuantificationQC: whichQuantification both", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotQuantificationQC(mydata, whichQuantification = "both", silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotQuantificationQC: whichQuantification corrected", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotQuantificationQC(mydata, whichQuantification = "corrected", silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotQuantificationQC: whichQuantification raw", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotQuantificationQC(mydata, whichQuantification = "raw", silent = TRUE)
  expect_s3_class(p1, "gg")
})

test_that("PlotQuantificationQC: whichPeptide raw", {
  skip_if_not(exists("mydata"), "User data not loaded")
  testPeptides <- c(unique(mydata$PTMTable$ModifiedPeptide)[1:30], NA)
  p1 <- PlotQuantificationQC(mydata, whichPeptide = testPeptides, silent = TRUE)
  expect_s3_class(p1, "gg")
})
