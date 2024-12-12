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

#MSFraggerImporter(path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/NGlycosylation",
                  #                  annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/NGlycoV1/annotation.csv",
                  #                  fastaPath = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/2024-02-25-decoys-2024-02-25-SeerDB-SI4.fasta.fas"
                  #                  peptideScoreCutoff = 40,
                  #                  glycanScoreCutoff = -100)

testcase <- subset(testdf$PSMTable, UniprotIDs == "P09871,F8WCZ6" & Run == "240419_ES_R00006_MG_c00001_SA_NPB_1_i1")
testcase <- testcase %>%
  dplyr::filter(!is.na(AssignedModifications) & AssignedModifications != "") %>%
  tidyr::separate_rows(AssignedModifications, sep = ",") %>%
  dplyr::filter(!grepl("C\\(57.0215\\)", AssignedModifications)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(PeptidePTMLocalization = as.numeric(regmatches(AssignedModifications, gregexpr("^[0-9]+", AssignedModifications))[[1]]),
                ProteinPTMLocalization = PeptidePTMLocalization + ProteinStart,
                ModificationSite = sub(".*([A-Za-z])\\(.*", "\\1", AssignedModifications),
                ModificationID = paste0(ModificationSite,ProteinPTMLocalization))

ggplot2::ggplot(testcase) +
  geom_line(data = data.frame(x = seq(1,850), y = 1), aes(x = x, y = y), size = 2) +
  geom_point(data = testcase, aes(x= ProteinPTMLocalization, y = 1), color = "white", fill = "red") +
  ggrepel::geom_label_repel(data = distinct(testcase[c("ProteinPTMLocalization", "ModificationID")]), aes(x =ProteinPTMLocalization, y = 1, label = ModificationID)) +
  theme_void()
