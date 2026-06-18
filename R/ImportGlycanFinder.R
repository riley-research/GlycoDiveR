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
                                 filterForNoNSequon,
                                 confidenceLevel,
                                 deltaModCutoff = FALSE,
                                 searchEngine = "GlycanFinder",
                                 SScoreCutoff,
                                 structureConfidenceCutoff,
                                 AScoreCutoff)

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
