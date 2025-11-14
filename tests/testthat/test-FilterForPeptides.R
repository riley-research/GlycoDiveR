test_that("FilterForPeptides PSMTable vector test", {
  skip_if_not(exists("mydata"), "User data not loaded")
  peptideSubsetLength <- ceiling(0.5 * length(unique(mydata$PSMTable$ModifiedPeptide)))
  peptideSubset <- unique(mydata$PSMTable$ModifiedPeptide)[1:peptideSubsetLength]
  returnedPeptides <- FilterForPeptides(mydata$PSMTable,peptideSubset)$ModifiedPeptide

  expect_setequal(returnedPeptides, intersect(returnedPeptides, peptideSubset))
  expect_true(all(returnedPeptides %in% peptideSubset))
})

test_that("FilterForPeptides PTMTable vector test", {
  skip_if_not(exists("mydata"), "User data not loaded")
  peptideSubsetLength <- ceiling(0.5 * length(unique(mydata$PTMTable$ModifiedPeptide)))
  peptideSubset <- unique(mydata$PTMTable$ModifiedPeptide)[1:peptideSubsetLength]
  returnedPeptides <- FilterForPeptides(mydata$PTMTable,peptideSubset)$ModifiedPeptide

  expect_setequal(returnedPeptides, intersect(returnedPeptides, peptideSubset))
  expect_true(all(returnedPeptides %in% peptideSubset))
})

test_that("FilterForPeptides PSMTable dataframe test", {
  skip_if_not(exists("mydata"), "User data not loaded")
  peptideSubsetLength <- ceiling(0.5 * length(unique(mydata$PSMTable$ModifiedPeptide)))
  peptideSubset <- data.frame(ModifiedPeptide = unique(mydata$PSMTable$ModifiedPeptide)[1:peptideSubsetLength],
                              extraColumn = NA)
  returnedPeptides <- FilterForPeptides(mydata$PSMTable,peptideSubset)$ModifiedPeptide

  expect_setequal(returnedPeptides, intersect(returnedPeptides, peptideSubset$ModifiedPeptide))
  expect_true(all(returnedPeptides %in% peptideSubset$ModifiedPeptide))
})

test_that("FilterForPeptides PTMTable dataframe test", {
  skip_if_not(exists("mydata"), "User data not loaded")
  peptideSubsetLength <- ceiling(0.5 * length(unique(mydata$PTMTable$ModifiedPeptide)))
  peptideSubset <- data.frame(ModifiedPeptide = unique(mydata$PTMTable$ModifiedPeptide)[1:peptideSubsetLength],
                              extraColumn = NA)
  returnedPeptides <- FilterForPeptides(mydata$PTMTable,peptideSubset)$ModifiedPeptide

  expect_setequal(returnedPeptides, intersect(returnedPeptides, peptideSubset$ModifiedPeptide))
  expect_true(all(returnedPeptides %in% peptideSubset$ModifiedPeptide))
})
