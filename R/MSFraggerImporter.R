#source(here::here("R/ImportConverters.R"))
#source(here::here("R/GlycoDiveRUtils.R"))

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

  if(sum(!is.na(filtereddf$Intensity)) > 0){
    quantAvailable <- TRUE
  }else{
    quantAvailable <- FALSE
    }

  PTMdf <- PSMToPTMTable(filtereddf)

  filtereddf$TotalGlycanComposition <- sapply(filtereddf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])
  PTMdf$TotalGlycanComposition <- sapply(PTMdf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               glycoType = "N",
               searchEngine = "MSFragger",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               quantAvailable = quantAvailable)

  class(data)

  return(data)
}

#test <- MSFraggerImporter(path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/NGlycosylation",
                  #                  annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/annotation.csv",
                  #                  fastaPath = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/2024-02-25-decoys-2024-02-25-SeerDB-SI4.fasta.fas",
                  #                  peptideScoreCutoff = -1,
                  #                  glycanScoreCutoff = 10)

#MSFraggerImporter(path = "C:/Users/tveth/OneDrive - UW/Documents/Research/GlycoDiveR/Testdata/NGlycosylation",
#                  annotation = "C:/Users/tveth/OneDrive - UW/Documents/Research/GlycoDiveR/Testdata/NGlycosylation/annotation.csv",
#                  fastaPath = "C:/Users/tveth/OneDrive - UW/Documents/Research/GlycoDiveR/Testdata/NGlycosylation/2024-02-25-decoys-2024-02-25-SeerDB-SI4.fasta.fas",
#                  peptideScoreCutoff = 40,
#                  glycanScoreCutoff = -100)

testing = FALSE
if(testing){
  testcase <- testdf$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", AssignedModifications))
  extractt <- unique(testcase$Genes)[15]
  testcase <- testcase %>% dplyr::filter(Genes == extractt)

  ggplot2::ggplot(testcase) +
    geom_line(data = data.frame(x = seq(1,850), y = 1), aes(x = x, y = y), size = 2) +
    geom_point(data = testcase, aes(x= ProteinPTMLocalization, y = 1), color = "white", fill = "red") +
    ggrepel::geom_label_repel(data = distinct(testcase[c("ProteinPTMLocalization", "ModificationID")]), aes(x =ProteinPTMLocalization, y = 1, label = ModificationID)) +
    theme_void()

  PlotPTMQuantification(test, unique(test$PTMTable$UniprotIDs)[18])

  testing <- subset(testdf$PTMTable, ModifiedPeptide == "YVTSAPM[147]PEPQAPGR")
  testing <- subset(testing, Run == "240419_ES_R00006_MG_c00001_SA_NPB_2_i1")

  GetMeanTechReps(testing)

  colorScheme <- c("#BAA5CC", "#9ADCEE", "#BAD97C", "#EEAED0", "#FAD821",
                   "#94D8C3", "#F7B8D2", "#A7C7E7", "#FFE87C", "#C0E4D0",
                   "#A1A9F2", "#C1D87F", "#E3B7E2", "#B1D3C2", "#F9A9B6",
                   "#D1D2E3", "#A4EFA1", "#D9D07A", "#98C9C7", "#F4D1A1")

  usethis::use_data(colorScheme, internal = TRUE, overwrite = TRUE)
  }
