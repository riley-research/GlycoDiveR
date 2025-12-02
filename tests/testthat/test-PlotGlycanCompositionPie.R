test_that("PlotGlycanCompositionPie: grouping = technicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionPie(mydata, grouping = "technicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycanCompositionPie: grouping = biologicalReps", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionPie(mydata, grouping = "biologicalReps", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycanCompositionPie: grouping = condition", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionPie(mydata, grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotGlycanCompositionPie: protein = unique(mydata$PSMTable$UniprotIDs)[1:100]", {
  skip_if_not(exists("mydata"), "User data not loaded")
  p <- PlotGlycanCompositionPie(mydata, whichProtein = unique(mydata$PSMTable$UniprotIDs)[1:200],
                                grouping = "condition", silent = TRUE)
  expect_s3_class(p, "gg")
})
