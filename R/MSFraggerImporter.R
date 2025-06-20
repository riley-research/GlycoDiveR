#' Title
#'
#' @param path Path to search engine output folder
#' @param annotation Path to annotation file
#' @param fastaPath Path to FASTA file
#' @param peptideScoreCutoff Peptide score cutoff
#' @param glycanScoreCutoff Glycan score cutoff
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data
#'
#' @returns Formatted dataframes
#' @export
#'
#' @examples MSFraggerImporter(path = "Z:/Folder",
#' annotation = "Z:/Folder/annotation.csv",
#' fasta = "Z:/fasta.fasta",
#' peptideScoreCutoff = 0,
#' glycanScoreCutoff = 0.05,
#' scrape = FALSE)
MSFraggerImporter <- function(path, annotation, fastaPath, peptideScoreCutoff, glycanScoreCutoff,
                              scrape = TRUE){
  unfiltereddf <- data.frame()
  annotationdf <- utils::read.csv(annotation)

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("psm.tsv", fileList)]

  if(length(fileList) == 0){
    stop("No files found")
  }
  for(file in fileList){
    fmessage(paste0("Now importing: ", file))
    temptable <- data.table::fread(paste0(path, "/", file), sep = "\t", check.names = TRUE, fill = TRUE)
    #temptable <- utils::read.table(paste0(path, "/", file), sep = "\t", header = T, quote = "")
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

  if(scrape){
    fmessage("Now scraping Uniprot. Set 'scrape = FALSE' to the importer function to skip this step.\n
             This might take a while for large datasets.")

    #pb <- utils::txtProgressBar(min = 0, max = nrow(PTMdf), style = 3)

    #sapply(seq(1, nrow(PTMdf)), function(x) {
    #  GetUniprotGlycoInfo(as.character(PTMdf[,c("UniprotIDs", "ProteinPTMLocalization", "GlycanType")][x,][1]),
    #                      as.character(PTMdf[,c("UniprotIDs", "ProteinPTMLocalization", "GlycanType")][x,][2]),
    #                      as.character(PTMdf[,c("UniprotIDs", "ProteinPTMLocalization", "GlycanType")][x,][3]))
    #  setTxtProgressBar(pb, x)})
    #close(pb)

    PTMdf <- PTMdf %>%
      dplyr::group_by(UniprotIDs, ProteinPTMLocalization, GlycanType) %>%
      dplyr::mutate(UniprotGlycoEvidence = GetUniprotGlycoInfo(accVec = UniprotIDs, PTMLocalization = ProteinPTMLocalization, type = GlycanType))

    fmessage("Finished scrape.")
  }

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

#test <- MSFraggerImporter(path = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/CD59_N",
#                          annotation = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/CD59_N/annotation.csv",
#                          fastaPath = "C:/Users/tim_v/Documents/PostDoc/R_Glycopeptide/CD59_N/2025-01-12-decoys-contam-CD59_Leu26-Asn102_linkHis10.fasta.fas",
#                          peptideScoreCutoff = 0,
#                          glycanScoreCutoff = 0.01)

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
