test_that("PlotUpSet: default", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
  } else {
    expect_message(PlotUpSet(mydata, silent = FALSE))
  }
})

test_that("PlotUpSet: grouping biologicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, grouping = "biologicalReps", silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
  } else {
    expect_message(PlotUpSet(mydata, grouping = "biologicalReps", silent = FALSE))
  }
})

test_that("PlotUpSet: grouping technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, grouping = "technicalReps", silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
  } else {
    expect_message(PlotUpSet(mydata, grouping = "technicalReps", silent = FALSE))
  }
})

test_that("PlotUpSet: type all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, type = "all", silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
  } else {
    expect_message(PlotUpSet(mydata, type = "all", silent = FALSE))
  }
})

test_that("PlotUpSet: level all", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, level = "protein", silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
  } else {
    expect_message(PlotUpSet(mydata, level = "protein", silent = FALSE))
  }
})

test_that("PlotUpSet: nintersect 3 and 3000", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p1 <- PlotUpSet(mydata, nintersects = 3, silent = TRUE)
  p2 <- PlotUpSet(mydata, nintersects = 3000, silent = TRUE)
  if (class(p1) == "upset") {
    expect_s3_class(p1, "upset")
    expect_s3_class(p2, "upset")
  } else {
    expect_message(PlotUpSet(mydata, nintersects = 3, silent = FALSE))
    expect_message(PlotUpSet(mydata, nintersects = 3000, silent = FALSE))
  }
})
