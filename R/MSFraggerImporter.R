source(here::here("R/ImportConverters.R"))

MSFraggerImporter <- function(path, annotation, peptideScoreCutoff, glycanScoreCutoff){
  #path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlyco"
  #annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlyco/annotation.csv"
  unfiltereddf <- data.frame()
  annotationdf <- utils::read.csv(annotation)

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("psm", fileList)]

  for(file in fileList){
    temptable <- utils::read.table(paste0(path, "/", file), sep = "\t", header = T)
    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
  }

  lengthList <- length(strsplit(unfiltereddf$Spectrum.File[1], "\\", fixed = T)[[1]])
  unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][lengthList-1])

  filtereddf <- MSFraggerConverter(unfiltereddf, annotationdf)

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               annotation = annotationdf,
               glycoType = "N",
               searchEngine = "MSFragger",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff)

  class(data)

  return(data)
}

#MSFraggerImporter(path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlyco",
#                  annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlyco/annotation.csv",
#                  peptideScoreCutoff = 40,
#                  glycanScoreCutoff = -100)
