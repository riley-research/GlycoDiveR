#' ImportMSFragger
#'
#' @param path Path to search engine output folder
#' @param annotation Path to annotation file
#' @param fastaPath Path to FASTA file
#' @param peptideScoreCutoff Peptide score cutoff
#' @param glycanScoreCutoff Glycan score cutoff
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param normalization The (glyco)peptide normalization used.
#' Choose between "none" or "median" (default = "median").
#' @param convertFPModCodeToMass MSFragger uses modification code in peptide
#' modified sequences. This replaces the code with the mass of the modification.
#' Keep this enabled when importing MSstats comparison results as MSstats uses
#' modification masses instead of modification codes. Default = TRUE
#' @param filterForNoNSequon Filter for peptides without an N-sequon. Only works
#' for OPair data (default = FALSE).
#' @param OPairLevelConversion Convert the OPair site probability levels to glycan
#' q-scores. Provide a vector that has values for each level. The default is
#' c(0,0,0.05,0.1), which means level1 = 0, level1b = 0, level2 = 0.05, and
#' level3 = 0.1). This qscore filtering is used with the Importers
#' glycanScoreCutoff argument.
#'
#' @returns Formatted GlycoDiveR data file.
#' @export
#'
#' @examples \dontrun{ImportMSFragger(path = "Z:/Folder",
#' annotation = "Z:/Folder/annotation.csv",
#' fasta = "Z:/fasta.fasta",
#' peptideScoreCutoff = 0,
#' glycanScoreCutoff = 0.05,
#' scrape = FALSE)}
ImportMSFragger <- function(path, annotation, fastaPath, peptideScoreCutoff, glycanScoreCutoff,
                            scrape = FALSE, normalization = "median",
                            convertFPModCodeToMass = TRUE, filterForNoNSequon = FALSE,
                            OPairLevelConversion = c(0,0,0.05,0.1)){
  unfiltereddf <- data.frame()
  quantdf <- data.frame()
  annotationdf <- utils::read.csv(annotation)
  CheckAnnotation(annotationdf)

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

  if(normalization %in% c("FP_Normalized", "FP_MaxLFQ")){
    quantPath <- list.files(path, recursive = TRUE)
    quantPath <- quantPath[grepl("combined_modified_peptide.tsv", quantPath)]

    if(length(quantPath) == 0){
      stop("No combined_modified_peptide.tsv files found. Select the right location or using different quantification.")
    }
    for(quant in quantPath){
      fmessage(paste0("Now importing: ", quant))
      temptable <- data.table::fread(paste0(path, "/", quant), sep = "\t", check.names = TRUE, fill = TRUE)
      quantdf <- plyr::rbind.fill(quantdf, temptable)
    }
  }

  filtereddf <- MSFraggerConverter(unfiltereddf, annotationdf, fastaPath, quantdf,
                                   scrape, normalization, convertFPModCodeToMass,
                                   OPairLevelConversion)

  PTMdf <- PSMToPTMTable(filtereddf)

  filtereddf$TotalGlycanComposition <- sapply(filtereddf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])
  PTMdf$TotalGlycanComposition <- sapply(PTMdf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "MSFragger",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = filterForNoNSequon)

  return(data)
}
