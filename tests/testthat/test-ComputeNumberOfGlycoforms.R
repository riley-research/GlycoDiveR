test_that("ComputeNumberOfGlycoforms", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- ComputeNumberOfGlycoforms(mydata, whichProtein = unique(mydata$PSMTable$UniprotIDs)[1], silent = TRUE)
  expect_true(is.character(p) || is.numeric(p))
})
