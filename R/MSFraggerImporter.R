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
#' @examples \dontrun{MSFraggerImporter(path = "Z:/Folder",
#' annotation = "Z:/Folder/annotation.csv",
#' fasta = "Z:/fasta.fasta",
#' peptideScoreCutoff = 0,
#' glycanScoreCutoff = 0.05,
#' scrape = FALSE)}
MSFraggerImporter <- function(path, annotation, fastaPath, peptideScoreCutoff, glycanScoreCutoff,
                              scrape = FALSE){
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
