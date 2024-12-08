source(here::here("R/ImportConverters.R"))

MSFraggerImporter <- function(path, annotation, fastaPath, peptideScoreCutoff, glycanScoreCutoff){
  unfiltereddf <- data.frame()
  annotationdf <- utils::read.csv(annotation)

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("psm.tsv", fileList)]

  for(file in fileList){
    temptable <- utils::read.table(paste0(path, "/", file), sep = "\t", header = T)
    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
  }

  lengthList <- length(strsplit(unfiltereddf$Spectrum.File[1], "\\", fixed = T)[[1]])
  unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][lengthList-1])

  filtereddf <- MSFraggerConverter(unfiltereddf, annotationdf, fastaPath)

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

#C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/NGlycosylation
#"C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/annotation.csv

#MSFraggerImporter(path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/NGlycosylation",
                  #                  annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/annotation.csv",
                  #                  fastaPath = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/2024-02-25-decoys-2024-02-25-SeerDB-SI4.fasta.fas"
                  #                  peptideScoreCutoff = 40,
                  #                  glycanScoreCutoff = -100)
