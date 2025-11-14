test_that("PlotUpSet: default", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, silent = TRUE)
  expect_s3_class(p1, "upset")
})

test_that("PlotUpSet: grouping biologicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p1, "upset")
})

test_that("PlotUpSet: grouping technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p1, "upset")
})

test_that("PlotUpSet: type all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, type = "all", silent = TRUE)
  expect_s3_class(p1, "upset")
})

test_that("PlotUpSet: level all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, level = "protein", silent = TRUE)
  expect_s3_class(p1, "upset")
})

test_that("PlotUpSet: nintersect 3 and 3000", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, nintersects = 3, silent = TRUE)
  p2 <- PlotUpSet(mydata, nintersects = 3000, silent = TRUE)
  expect_s3_class(p1, "upset")
  expect_s3_class(p2, "upset")
})
