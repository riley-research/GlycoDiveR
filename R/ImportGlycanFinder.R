#' ImportGlycanFinder
#'
#' This imports search results output by GlycanFinder. It will work for
#' N-glycopeptide searches and O-glycopeptides searches. During import, this
#' function will look for all the glycan.glycopsms.csv files in the
#' specified folder and all the subfolders.
#'
#' @param path Path to search engine output folder
#' @param annotation Path to annotation file
#' @param fastaPath Path to FASTA file
#' @param peptideScoreCutoff Peptide score cutoff, which is the Glycopeptide Score
#' column in the search engine output file.
#' @param glycanScoreCutoff Glycan score cutoff, which is the Glycan Score
#' column in the search engine output file.
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param normalization The (glyco)peptide normalization used.
#' Choose between "median" (default) or "none".
#' @param cutoffFilter Filter for the glycan score, peptide score, and
#' confidence levels during data import
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
#' @param useExtendedOGlycanCategories set to FALSE will classify all O-glycans as OGlycan
#' Set to TRUE will classify them in "Sialofucosylated", "Sialylated", "Fucosylated",
#' and "OGlycan".
#' @param GlyToucan set to TRUE will connect to the GlyCosmos API to retrieve
#' GlyToucan identifiers.
#' @param SScoreCutoff The lower limit of the S score used to filter the data, e.g.
#' SScoreCutoff = 10. This will retain S Scores that are 10 or higher.
#' @param structureConfidenceCutoff The Structure Cutoffs to be retained. For example,
#' structureConfidenceCutoff = c("LOW", "HIGH")
#' @param AScoreCutoff The lower limit of the A score used to filter the data, e.g.
#' AScoreCutoff = 10. This will retain A Scores 10 or higher. The highest score is
#' retained for For multiply glycosylated peptides.
#'
#' @returns A formatted GlycoDiveR file.
#' @export
#'
#' @examples
#' \dontrun{
#' ImportGlycanFinder(
#'   path = "Z:/Analysis 2",
#'   annotation = "Z:/annotation.csv",
#'   fastaPath = "Z:/glycan.proteins.fasta",
#'   normalization = "none",
#'   GlyToucan = FALSE,
#'   structureConfidenceCutoff = "HIGH",
#'   AScoreCutoff = 10
#' )
#' }
ImportGlycanFinder <- function(path, annotation, fastaPath, peptideScoreCutoff = 0,
                            glycanScoreCutoff = 0, scrape = TRUE, normalization = "median",
                            cutoffFilter = TRUE, dropNoQuant = FALSE, minPeptideCoverage = FALSE,
                            thresholdMode = c("group", "total"), useExtendedOGlycanCategories = FALSE,
                            GlyToucan = TRUE, SScoreCutoff = FALSE, structureConfidenceCutoff = FALSE,
                            AScoreCutoff = FALSE){
  unfiltereddf <- data.frame()
  quantdf <- data.frame()
  annotationdf <- utils::read.csv(annotation,
                                  colClasses = c(
                                    Run = "character",
                                    Condition = "character",
                                    Alias = "character",
                                    BioReplicate = "integer",
                                    TechReplicate = "integer"
                                  ))
  CheckAnnotation(annotationdf)

  if(identical(useExtendedOGlycanCategories, TRUE)){
    .modEnv$useExtendedOGlycanCategories <- useExtendedOGlycanCategories
  }else{.modEnv$useExtendedOGlycanCategories <- FALSE}

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("glycan.glycopsms.csv", fileList)]

  if(length(fileList) == 0){
    stop("No files found")
  }


  for(file in fileList){
    fmessage(paste0("Now importing: ", file))
    temptable <- data.table::fread(paste0(path, "/", file), sep = ",", check.names = TRUE, fill = TRUE)

    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
  }

  unfiltereddf$ID <-   seq(1:nrow(unfiltereddf))

  filtereddf <- GlycanFinderConverter(unfiltereddf, annotationdf, fastaPath, scrape,
                        normalization, GlyToucan)

  if(cutoffFilter){
    filtereddf <- FilterPSMTable(filtereddf,
                                 peptideScoreCutoff,
                                 glycanScoreCutoff,
                                 filterForNoNSequon = FALSE,
                                 confidenceLevels = FALSE,
                                 deltaModCutoff = FALSE,
                                 searchEngine = "GlycanFinder",
                                 SScoreCutoff,
                                 structureConfidenceCutoff,
                                 AScoreCutoff)

    if(!all(is.na(filtereddf$Intensity)) &
       sum(filtereddf$Intensity, na.rm = TRUE) != 0 &
       normalization == "median") {
      globalMedian = stats::median(filtereddf$Intensity[filtereddf$Intensity != 0], na.rm = TRUE)
      filtereddf <- filtereddf %>%
        dplyr::mutate(
          .by = .data$Run,
          Intensity = medianNormalization(
            intensityVec = .data$Intensity,
            globalMedian = globalMedian
          )
        )
      filtereddf <- filtereddf %>%
        dplyr::mutate(Intensity = dplyr::coalesce(.data$Intensity, 0))

      fmessage("Successfully median normalized after quality filtering.")
    }
  }

  if(dropNoQuant){filtereddf <- filtereddf %>% dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0)}

  if(!identical(minPeptideCoverage, FALSE)){
    filtereddf <- FilterForMinPeptides(filtereddf, minPeptideCoverage, thresholdMode)

    if(!all(is.na(filtereddf$Intensity)) &
       sum(filtereddf$Intensity != 0) &
       normalization == "median") {
      globalMedian = stats::median(filtereddf$Intensity[filtereddf$Intensity != 0], na.rm = TRUE)
      filtereddf <- filtereddf %>%
        dplyr::mutate(
          .by = .data$Run,
          Intensity = medianNormalization(
            intensityVec = .data$Intensity,
            globalMedian = globalMedian
          )
        )
      filtereddf <- filtereddf %>%
        dplyr::mutate(Intensity = dplyr::coalesce(.data$Intensity, 0))

      fmessage("Successfully median normalized after quality filtering.")
    }
  }

  PTMdf <- PSMToPTMTable(filtereddf)

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "GlycanFinder",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = FALSE,
               confidenceLevels = FALSE,
               deltaModCutoff = FALSE,
               SScoreCutoff = SScoreCutoff,
               structureConfidenceCutoff = structureConfidenceCutoff,
               AScoreCutoff = AScoreCutoff)

  return(data)
}
