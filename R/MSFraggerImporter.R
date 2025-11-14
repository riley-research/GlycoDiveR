#' Title
#'
#' @param path Path to search engine output folder
#' @param annotation Path to annotation file
#' @param fastaPath Path to FASTA file
#' @param peptideScoreCutoff Peptide score cutoff
#' @param glycanScoreCutoff Glycan score cutoff
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data
#' @param normalization The (glyco)peptide modification.
#' Choose between "none" or "median"
#' @param convertFPModCodeToMass MSFragger uses modification code in peptide
#' modified sequences. This replaces the code with the mass of the modification.
#' Keep this enabled when importing MSstats comparison results as MSstats uses
#' modification masses instead of modification codes. Default = TRUE
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
                              scrape = FALSE, normalization = "median", convertFPModCodeToMass = TRUE){
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
    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
  }

  unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][length(strsplit(x, "\\", fixed = T)[[1]])-1])

  filtereddf <- MSFraggerConverter(unfiltereddf, annotationdf, fastaPath,
                                   scrape, normalization, convertFPModCodeToMass)

  PTMdf <- PSMToPTMTable(filtereddf)

  filtereddf$TotalGlycanComposition <- sapply(filtereddf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])
  PTMdf$TotalGlycanComposition <- sapply(PTMdf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "MSFragger",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff)

  return(data)
}
