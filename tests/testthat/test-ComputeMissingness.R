test_that("ComputeMissingness all data", {
  skip_if_not(exists("mydata"), "User data not loaded")
  rslt <- ComputeMissingness(mydata, type = "all", silent = TRUE)
  expect_type(rslt, "character")
})

test_that("ComputeMissingness glyco data", {
  skip_if_not(exists("mydata"), "User data not loaded")
  rslt <- ComputeMissingness(mydata, type = "glyco", silent = TRUE)
  expect_type(rslt, "character")
})
