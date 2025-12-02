#' ImportMSFragger
#'
#' @param path Path to search engine output folder
#' @param annotation Path to annotation file
#' @param fastaPath Path to FASTA file
#' @param peptideScoreCutoff Peptide score cutoff
#' @param glycanScoreCutoff Glycan score cutoff
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param normalization The (glyco)peptide normalization used.
#' Choose between "median" (default), FP_Normalized, FP_MaxLFQ, or none.
#' median: performs median normalization.
#' FP_Normalized: extracts the intensity values in the Intensity columns of the
#' combined_modified_peptide.tsv files.
#' FP_MaxLFQ: extracts the intensity values in the MaxLFQ.Intensity columns of the
#' combined_modified_peptide.tsv files.
#' none: uses the raw intensity values from the PSM.tsv files.
#' @param convertFPModCodeToMass MSFragger uses modification code in peptide
#' modified sequences. This replaces the code with the mass of the modification.
#' Keep this enabled when importing MSstats comparison results as MSstats uses
#' modification masses instead of modification codes. Default = TRUE
#' @param filterForNoNSequon Filter for peptides without an N-sequon. Only works
#' for OPair data (default = FALSE).
#' @param confidenceLevel What OPair confidence levels to accept, options are
#' Level1, Level1b, Level2, Level3 (default = FALSE). Provide like this:
#' confidenceLevel = c("Level1", "Level1b")
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
ImportMSFragger <- function(path, annotation, fastaPath, peptideScoreCutoff = 0,
                            glycanScoreCutoff = 0.01, scrape = TRUE, normalization = "median",
                            convertFPModCodeToMass = TRUE, filterForNoNSequon = FALSE,
                            confidenceLevel = FALSE){
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
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = filterForNoNSequon,
               confidenceLevels = confidenceLevel,
               deltaModCutoff = FALSE)

  return(data)
}
