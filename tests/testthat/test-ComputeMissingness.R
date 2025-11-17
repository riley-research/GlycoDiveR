test_that("PlotPSMCount technicalReps works", {
  skip_if_not(exists("mydata"), "User data not loaded")
  rslt <- ComputeMissingness(mydata, silent = TRUE)
  expect_type(rslt, "character")
})
