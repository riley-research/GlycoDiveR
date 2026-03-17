#' ImportpGlyco
#'
#' This is the importer for pGlyco data. It needs to untouched pGlyco output, which
#' is a .txt file named 'pGlycoDB-GP-FDR-Pro-Quant-Site.txt'. This script will look
#' for all files in the folder and all subfolders.
#'
#' @param path The folder path, which will be used to find all
#' 'pGlycoDB-GP-FDR-Pro-Quant-Site.txt' files in this folder and all subfolders.
#' @param annotation The annotation dataframe. Please generate a template using
#' GetAnnotationTemplate("path", tool = "pGlyco")
#' @param fastaPath The path to the fasta file that was used.
#' @param peptideScoreCutoff The score cutoff for the "PeptideFDR" column. The cutoff
#' is the upper limit.
#' @param glycanScoreCutoff The score cutoff for the "GlycanFDR" column. The score is the
#' upper limit.
#' @param normalization The (glyco)peptide normalization used.
#' Choose between "none" or "median" (default = "median").
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param cutoffFilter Filter for the glycan score, peptide score, N Sequon, and
#' confidence levels.
#' @param dropNoQuant Remove or keep PSMs without a quantitative value or have a quant
#' of 0.
#' @param minPeptideCoverage Works together with the thresholdMode argument.
#' If thresholdMode = "total", provide a single numeric value. Values between 0
#' and 1 are interpreted as proportions (e.g., 0.1 = 10%), while values >= 1 are
#' interpreted as absolute counts (e.g., 2 = at least two peptides identified).
#' If thresholdMode = "group", provide a numeric vector of length two (e.g., c(0.1, 3)).
#' The first value is interpreted in the same way as for "total" (proportion or
#' absolute count). The second value specifies the minimum number of groups in
#' which this threshold must be met.
#' @param thresholdMode Character string specifying the filtering mode: "total"
#' applies a global threshold across all data, while "group" identifies the threshold
#' within groups.
#' @param GlyToucan set to TRUE will connect to the GlyCosmos API to retrieve
#' GlyToucan identifiers.
#'
#' @returns Formatted GlycoDiveR data file.
#' @export
#'
#' @examples \dontrun{ImportpGlyco(path = "Z:/Folder",
#' annotation = "Z:/Folder/annotation.csv",
#' fasta = "Z:/fasta.fasta",
#' peptideScoreCutoff = 0.01,
#' glycanScoreCutoff = 0.01,
#' scrape = FALSE)}
ImportpGlyco <- function(path, annotation, fastaPath, peptideScoreCutoff = 0.01,
                         glycanScoreCutoff = 0.01, normalization = "median", scrape = TRUE,
                         cutoffFilter = TRUE, dropNoQuant = FALSE, minPeptideCoverage = FALSE,
                         thresholdMode = c("group", "total"), GlyToucan = TRUE){
  unfiltereddf <- data.frame()
  annotationdf <- utils::read.csv(annotation)
  CheckAnnotation(annotationdf)

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("pGlycoDB-GP-FDR-Pro-Quant-Site.txt", fileList)]
  unfiltereddf <- data.frame()

  if(length(fileList) == 0){stop("No files found.")}
  for(file in fileList){
    print(paste0(path, "/", file))
    temptable <- data.table::fread(paste0(path, "/", file), sep = "\t", check.names = TRUE, fill = TRUE)
    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
  }

  unfiltereddf$ID <-   seq(1:nrow(unfiltereddf))

  modification_df <- unfiltereddf %>%
    dplyr::select("Mod") %>%
    tidyr::separate_longer_delim(cols = "Mod", delim=";") %>%
    dplyr::filter(.data$Mod != "" & !is.na(.data$Mod)) %>%
    dplyr::mutate(Mod = stringr::str_extract(.data$Mod, "(?<=,).*?(?=\\[)")) %>%
    dplyr::distinct() %>%
    dplyr::rename("FullName" = "Mod") %>%
    dplyr::left_join(.modEnv$ModificationDatabase, by = "FullName")

  missing_masses <- modification_df %>%
    dplyr::filter(is.na(.data$ModificationMass) | .data$ModificationMass == "")

  if(nrow(missing_masses) > 0){
    stop(
      "Missing masses for the following modifications. Please refer to the XX function how to add those:\n",
      paste(missing_masses$FullName, collapse = ", ")
    )
  }

  filtereddf <- pGlycoConverter(unfiltereddf, annotationdf, fastaPath, modification_df,
                  normalization, scrape, GlyToucan)

  if(cutoffFilter){
    filtereddf <- FilterPSMTable(filtereddf,
                                 peptideScoreCutoff,
                                 glycanScoreCutoff,
                                 filterForNoNSequon = FALSE,
                                 confidenceLevels = FALSE,
                                 deltaModCutoff = FALSE,
                                 searchEngine = "pGlyco")}

  if(dropNoQuant){filtereddf <- filtereddf %>% dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0)}

  if(!identical(minPeptideCoverage, FALSE)){filtereddf <- FilterForMinPeptides(filtereddf, minPeptideCoverage, thresholdMode)}

  PTMdf <- PSMToPTMTable(filtereddf)

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "pGlyco",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = FALSE,
               confidenceLevels = FALSE,
               deltaModCutoff = FALSE)

  return(data)
  }
