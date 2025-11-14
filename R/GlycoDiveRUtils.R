CombineTwoCols <- function(row){
  if (is.na(row[2]) | row[2] == ""){
    returnVal <- as.character(row[1])
    return(returnVal)
  }else{
    returnVal <- as.character(paste(row[1], row[2], sep = ","))
    if(substring(returnVal, 1, 1) == ","){
      returnVal <- substring(returnVal, 2, nchar(returnVal))
    }
    return(returnVal)
    }}

CombineProtCols <- function(row){
  if (is.na(row[2]) | row[2] == ""){
    return(row[1])
  }else{
    addcol <- unlist(regmatches(row[2], gregexpr("(?<=\\|)[^|]+(?=\\|)", row[2], perl = TRUE)))
    addcol <- addcol[!grepl(",|_", addcol)]
    return(as.character(paste(row[1], paste(addcol, collapse=","), sep = ",")))
  }
}

GetProteinLength <- function(IDVec, fastaFile){
  inputVector <- strsplit(IDVec[1], ",")[[1]]
  uniprotID <- inputVector[1]
  IDhit <- fastaFile[grepl(uniprotID, names(fastaFile)) & !grepl("rev", names(fastaFile))]
  protLength <- nrow(as.data.frame(IDhit[[1]]))
  return(as.integer(protLength))
}

PSMToPTMTable <- function(PSMTable){
  tempdf <- PSMTable %>%
    dplyr::filter(!is.na(.data$AssignedModifications) & .data$AssignedModifications != "") %>%
    tidyr::separate_rows(.data$AssignedModifications, sep = ",")

  tempdf$AssignedModifications <- gsub("N-term", "1", tempdf$AssignedModifications)

  tempdf <- tempdf %>%
    dplyr::mutate(
      PeptidePTMLocalization = as.numeric(stringr::str_extract(.data$AssignedModifications, "\\d+")),
      ProteinPTMLocalization = .data$PeptidePTMLocalization + .data$ProteinStart -1,
      ModificationSite = stringr::str_extract(.data$AssignedModifications, "[A-Za-z](?=\\()"),
      ModificationID = paste0(.data$ModificationSite, .data$ProteinPTMLocalization))

  tempdf$GlycanType <- apply(tempdf[,c("AssignedModifications", "TotalGlycanComposition")],
                             1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))

  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Generated PTM table.\033[0m")

  return(tempdf)
}

GlycanComptToGlycanType <- function(mod, glycanComp){
  modType <- c()
  allModsVec <- strsplit(mod, ",")[[1]]

  if(is.na(mod) | mod == ""){
    return("Unmodified")}

  for(i in allModsVec){
    modifiedResidue <- strsplit(i, "\\(")[[1]][1]
    modifiedResidue <- substr(modifiedResidue, nchar(modifiedResidue) , nchar(modifiedResidue))

    if(is.na(mod) | mod == ""){
      modType <- append(modType, "NonGlyco")
    }

    glycanMass <- sub(".*\\((.*)\\).*", "\\1", i)

    if(glycanMass %in% c("15.9949", "57.0214", "57.0215", "42.0106")){
      modType <- append(modType, "NonGlyco")
    }else if((glycanComp != "" & modifiedResidue == "N") | (!is.na(glycanComp) & modifiedResidue == "N")){
      if(TRUE %in% grepl(glycanMass, mod)){
        hexNAc_count <- suppressWarnings(as.numeric(sub(".*N\\(([0-9]+)\\).*", "\\1", glycanComp)))
        hex_count <- suppressWarnings(as.numeric(sub(".*H\\(([0-9]+)\\).*", "\\1", glycanComp)))

        glycanCat <- dplyr::case_when(
          grepl("A|G", glycanComp) & grepl("F", glycanComp) ~ "Sialyl+Fucose",
          grepl("A|G", glycanComp) ~ "Sialyl",
          grepl("F", glycanComp) ~ "Fucose",
          !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count == 2 & hex_count == 3 ~ "Paucimannose",
          (!is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 2) |
            ( !is.na(hexNAc_count) & !is.na(hex_count) & hex_count < 3) ~ "Truncated",
          !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 3 & hex_count > 3 ~ "High Mannose",
          TRUE ~ "Complex/Hybrid"
        )
        modType <- append(modType, glycanCat)
      }else{
        modType <- append(modType, "UndefinedGlyco")
      }}else if((glycanComp != "" & modifiedResidue %in% c("S", "T")) | (!is.na(glycanComp) & modifiedResidue %in% c("S", "T"))){
        modType <- append(modType, "OGlycan")
      }else if(glycanComp != "" & !is.na(glycanComp)){
        modType <- append(modType, "NonCanonicalGlyco")
      }else{
        modType <- append(modType, "NonGlyco")
      }
  }

  return(as.vector(modType))
}

GetPeptide <- function(pep, modpep){
  if(!is.na(modpep) & modpep != ""){
    return(modpep)}
  else {
    return(pep)
  }
}

GetMeanTechReps <- function(df){
  #Keep highest intensity per technical rep
  df <- df %>%
    dplyr::arrange(dplyr::desc(.data$Intensity)) %>%
    dplyr::distinct(.data$Run, .data$AssignedModifications, .data$ModifiedPeptide,
                    .data$Condition, .data$BioReplicate, .data$TechReplicate,
                    .keep_all = TRUE)

  if("ModificationID" %in% names(df)){
    #Take median of technical reps together
    df <- df %>%
      dplyr::mutate(.by = c("ModifiedPeptide", "AssignedModifications",
                            "Condition", "BioReplicate"),
                    Intensity = stats::median(.data$Intensity, na.rm = TRUE)) %>%
      dplyr::distinct(.data$ModifiedPeptide, .data$Condition, .data$BioReplicate,
                      .data$ModificationID, .keep_all = TRUE)

    return(df)
  }else{
    df <- df %>%
      dplyr::mutate(.by = c("ModifiedPeptide", "AssignedModifications",
                            "Condition", "BioReplicate"),
                    Intensity = stats::median(.data$Intensity, na.rm = TRUE)) %>%
      dplyr::distinct(.data$ModifiedPeptide, .data$Condition, .data$BioReplicate,
                      .keep_all = TRUE)
    return(df)
  }

}

CleanGlycanNames <- function(glycan){
  glycan <- gsub("HexNAc", "N", glycan)
  glycan <- gsub("Hex", "H", glycan)
  glycan <- gsub("NeuAc", "A", glycan)
  glycan <- gsub("Fuc", "F", glycan)
  glycan <- gsub("NeuGc", "G", glycan)
  return(glycan)
}

FilterForCutoffs <- function(input, silent = FALSE){
  if(input$searchEngine %in% c("MSFragger")){
    if(!silent){
      fmessage(paste0("Filtering for PSMScore >= ", input$peptideScoreCutoff, " and glycan score <= ", input$glycanScoreCutoff))
      fmessage(paste0("Prefilter number of rows PSM table: ", nrow(input$PSMTable), ". Prefilter number of rows PTM table: ", nrow(input$PTMTable)))
    }
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter((.data$PSMScore >= input$peptideScoreCutoff & .data$GlycanQValue <= input$glycanScoreCutoff) |
          (.data$PSMScore >= input$peptideScoreCutoff & is.na(.data$GlycanQValue)))

    input$PTMTable <- input$PTMTable %>%
      dplyr::filter((.data$PSMScore >= input$peptideScoreCutoff & .data$GlycanQValue <= input$glycanScoreCutoff) |
          (.data$PSMScore >= input$peptideScoreCutoff & is.na(.data$GlycanQValue)))

    if(!silent){
      fmessage(paste0("Postfilter number of rows PSM table: ", nrow(input$PSMTable), ". Postfilter number of rows PTM table: ", nrow(input$PTMTable)))
    }
    return(input)
  }else{
    warning("No search engine recognized, returning unfiltered dataframe.")
    return(input)
  }
}

GetGlycoSitesPerProtein <- function(IDVec, fastaFile){
  inputVector <- strsplit(IDVec[1], ",")[[1]]
  uniprotID <- inputVector[1]
  IDhit <- fastaFile[grepl(uniprotID, names(fastaFile)) & !grepl("rev", names(fastaFile))]
  seq <- paste(IDhit[[1]], collapse = "")
  seq <- toupper(seq)

  pattern_NXS <- "N[^P]S"
  pattern_NXT <- "N[^P]T"
  pattern_S <- "S"
  pattern_T <- "T"

  count_NXS <- length(regmatches(seq, gregexpr(pattern_NXS, seq))[[1]])
  count_NXT <- length(regmatches(seq, gregexpr(pattern_NXT, seq))[[1]])
  count_S <- length(regmatches(seq, gregexpr(pattern_S, seq))[[1]])
  count_T <- length(regmatches(seq, gregexpr(pattern_T, seq))[[1]])

  rslt <- paste0(sum(count_NXS, count_NXT), ";", sum(count_S, count_T))

  return(rslt)
}

fmessage <- function(m){
  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: ", m, "\033[0m")
}

GetUniprotGlycoInfo <- function(accVec, PTMLocalization, type){
  # if(type[1] == "" | is.na(type[1]) | type[1] %in% c("NonGlyco", "Unmodified")){
  #   return(NA)
  # }
  #
  # acc <- strsplit(accVec[1], ",")[[1]][1]
  # geturl <- paste0("https://rest.uniprot.org/uniprotkb/search?query=accession:", acc, "&format=tsv&fields=ft_carbohyd")
  # scrape <- read.csv(URLencode(geturl), sep = "\t")
  # scrape <- scrape %>%
  #   tidyr::separate_longer_delim(cols = "Glycosylation", delim = " CARBOHYD ") %>%
  #   dplyr::mutate(Glycosylation = gsub("CARBOHYD ", "", Glycosylation)) %>%
  #   tidyr::separate_wider_delim(cols = "Glycosylation", delim = "; /note=", names = c("Site", "Info"))
  #
  # scrape <- subset(scrape, Site == as.character(PTMLocalization[1]))
  #
  # if(is.null(scrape) | nrow(scrape) == 0 | is.null(nrow(scrape))){
  #   return("No evidence")
  # }

  return()
}

GetUniprotSubcellularInfo <- function(UniprotIDs){
  UniprotIDs_df <- data.frame(UniprotIDs = UniprotIDs, stringsAsFactors = FALSE) %>%
    dplyr::mutate(UID = purrr::map_chr(.data$UniprotIDs, ~ stringr::str_split(.x, ",")[[1]][1]))

  UniprotIDsClean <-  UniprotIDs_df %>%
    dplyr::distinct(.data$UID)

  tempdf <- UniprotR::GetSubcellular_location(UniprotIDsClean$UID)

  tempdf$SubcellularLocalization <- sapply(tempdf$Subcellular.location..CC., function(x) cleanSubcellularLocation(x))
  tempdf$TopoDomain <- sapply(tempdf$Topological.domain, function(x) getTopoDomain(x))
  tempdf$TransmembraneDomain <- sapply(tempdf$Transmembrane, function(x) getTransmembraneDomain(x))

  tempdf$UID <- rownames(tempdf)

  tempdf <- tempdf %>%
    dplyr::mutate(Domains = dplyr::case_when(is.na(.data$TopoDomain) & is.na(.data$TransmembraneDomain) ~ NA_character_,
                                             is.na(.data$TopoDomain) ~ .data$TransmembraneDomain,
                                             is.na(.data$TransmembraneDomain) ~ .data$TopoDomain,
                                             TRUE ~ paste(TopoDomain, TransmembraneDomain, sep = ";"))) %>%
    dplyr::select(.data$UID, .data$SubcellularLocalization, .data$Domains)

  UniprotIDs_df <- UniprotIDs_df %>%
    dplyr::left_join(tempdf, by = c("UID"))

  return(UniprotIDs_df[c("SubcellularLocalization", "Domains")])
}

cleanSubcellularLocation <- function(rawString){
  if(!grepl("SUBCELLULAR LOCATION" , rawString)){
    return(NA)
  }else if(grepl("\\{ECO:", rawString)){
    formattedString <- stringr::str_extract_all(rawString, "(?<=: |\\. )([^\\{]+?)(?= \\{ECO:)") %>%
      unlist() %>%
      stringr::str_trim() %>%
      unique()

    formattedString <- as.data.frame(formattedString) %>%
      dplyr::filter(!grepl("Note", .data$formattedString)) %>%
      tidyr::separate_longer_delim(cols = .data$formattedString, delim = ", ") %>%
      dplyr::mutate(formattedString = dplyr::if_else(grepl("\\:", .data$formattedString),
                                              stringr::str_extract(.data$formattedString, "(?<=\\]: ).*"),
                                              formattedString),
        formattedString = stringr::str_to_title(.data$formattedString)) %>%
      dplyr::distinct() %>%
      dplyr::pull(.data$formattedString) %>%
      paste(collapse = ";")

    return(formattedString)
  }else if(stringr::str_count(rawString, "LOCATION") > 1){
    formattedString <- stringr::str_extract_all(rawString, "(?<=\\]: ).*?(?=\\.|;|$)") %>%
      unlist()

    formattedString <- as.data.frame(formattedString) %>%
      tidyr::separate_longer_delim(cols = .data$formattedString, delim = ", ") %>%
      dplyr::mutate(formattedString = dplyr::if_else(grepl("\\:", .data$formattedString),
                                                     stringr::str_extract(.data$formattedString, "(?<=\\]: ).*"),
                                                     formattedString),
                    formattedString = stringr::str_to_title(.data$formattedString)) %>%
      dplyr::distinct() %>%
      dplyr::pull(.data$formattedString) %>%
      paste(collapse = ";")

    return(formattedString)
  }else{
    formattedString <- strsplit(rawString, "LOCATION: ")[[1]][2]

    if(any(grepl(",", formattedString))){
      formattedString <- as.data.frame(formattedString) %>%
        tidyr::separate_longer_delim(cols = .data$formattedString, delim = ", ") %>%
        dplyr::mutate(formattedString = dplyr::if_else(grepl("\\:", .data$formattedString),
                                                       stringr::str_extract(.data$formattedString, "(?<=\\]: ).*"),
                                                       formattedString),
                      formattedString = stringr::str_to_title(.data$formattedString)) %>%
        dplyr::distinct() %>%
        dplyr::pull(.data$formattedString) %>%
        paste(collapse = ";")
    }
    formattedString <- gsub("\\.", "", formattedString)

    return(stringr::str_to_title(formattedString))
    }
}

getTopoDomain <- function(rawString){
  formattedString <- ""

  if(is.na(rawString) | rawString == ""){
    return(NA)
  }

  rawStringVec <- strsplit(rawString, "TOPO_DOM")[[1]]

  for(i in 2:length(rawStringVec)){
    firstAA <- 0
    lastAA <- 0
    domainType <- ""
    tempString <- trimws(rawStringVec[i])

    firstAA <- sub("\\..*", "", tempString)
    lastAA <- sub(".*\\.\\.(\\d+).*", "\\1", tempString)
    domainType <- sub(".*note=([^;]+);.*", "\\1", tempString)

    if(formattedString == ""){
      formattedString <- paste0(domainType, "(", firstAA, "-", lastAA, ")")
    }else{
      formattedString <- paste(formattedString, paste0(domainType, "(", firstAA, "-", lastAA, ")"), sep = ";")
    }
  }
return(formattedString)
}

getTransmembraneDomain <- function(rawString){
  formattedString <- ""

  if(is.na(rawString) | rawString == ""){
    return(NA)
  }

  rawStringVec <- strsplit(rawString, "TRANSMEM")[[1]]

  for(i in 2:length(rawStringVec)){
    firstAA <- 0
    lastAA <- 0
    tempString <- trimws(rawStringVec[i])

    firstAA <- sub("\\..*", "", tempString)
    lastAA <- sub(".*\\.\\.(\\d+).*", "\\1", tempString)
    domainType <- sub(".*note=([^;]+);.*", "\\1", tempString)

    if(formattedString == ""){
      formattedString <- paste0("Transmembrane domain", "(", firstAA, "-", lastAA, ")")
    }else{
      formattedString <- paste(formattedString, paste0("Transmembrane domain", "(", firstAA, "-", lastAA, ")"), sep = ";")
    }
  }
  return(formattedString)
}

calculateElbowCoords <- function(xVec, yVec, return = "x"){
  coord1_x <- as.numeric(strsplit(xVec, ";")[[1]][1])
  coord1_y <- as.numeric(strsplit(yVec, ";")[[1]][1])

  coord4_x <- as.numeric(strsplit(xVec, ";")[[1]][2])
  coord4_y <- as.numeric(strsplit(yVec, ";")[[1]][2])

  yMean <- mean(c(coord1_y, coord4_y), na.rm = TRUE)

  coord2_x <- coord1_x
  coord2_y <- yMean

  coord3_x <- coord4_x
  coord3_y <- yMean

  if(return == "x"){
    return(paste(coord1_x, coord2_x, coord3_x, coord4_x, sep = ";"))
  }else{
    return(paste(coord1_y, coord2_y, coord3_y, coord4_y, sep = ";"))
  }
}

medianNormalization <- function(intensityVec, globalMedian){
  nonzero <- intensityVec != 0 & !is.na(intensityVec)
  intensityVec_log2 <- rep(NA_real_, length(intensityVec))
  intensityVec_log2[nonzero] <- log2(intensityVec[nonzero])
  localMedian <- stats::median(intensityVec_log2[nonzero], na.rm = TRUE)
  deltaMedian <- log2(globalMedian) - localMedian
  intensityVec_log2[nonzero] <- intensityVec_log2[nonzero] + deltaMedian
  return(2^intensityVec_log2)
}

FPModCodeToModMass <- function(modifiedPep, assignedMods){
  #1. Generate a dataframe with the modification codes and modification masses
  #2. Look over to rows to replace the codes with the masses
  # C-term mods should be included too (but are not yet)
  # FPModCodeToModMass(modifiedPep = mydata$PSMTable$ModifiedPeptide, assignedMods = mydata$PSMTable$AssignedModifications)
  modLookupTable <- data.frame(modCode = as.character(),
                               modMass = as.character())

  tempdf <- data.frame(modifiedPep = modifiedPep, assignedMods = assignedMods)

  modTable <- tempdf %>%
    dplyr::filter(!is.na(.data$assignedMods) & .data$assignedMods != "") %>%
    tidyr::separate_longer_delim(cols = assignedMods, delim = ",") %>%
    dplyr::filter(!grepl("C\\(57\\.02", .data$assignedMods))

  for(i in 1:nrow(modTable)){
    AssignedMod <- modTable[i,"assignedMods"]
    ModifiedPep <- modTable[i,"modifiedPep"]

    modMassi <- stringr::str_extract(AssignedMod, "(?<=\\()[^)]*(?=\\))")

    if(grepl("N-term", AssignedMod)){
      splitAA <- 3
    }else{
      splitAA <- as.numeric(stringr::str_extract(AssignedMod, "^[0-9]+"))
      splitAA <- gregexpr("[A-Z]", ModifiedPep)[[1]][splitAA] + 2
    }

    modCodei <- substr(ModifiedPep, splitAA, nchar(ModifiedPep))
    modCodei <- stringr::str_extract(modCodei, "^[0-9]+")

    if(!any(modLookupTable$modCode == modCodei & modLookupTable$modMass == modMassi)){
      modLookupTable <- rbind(modLookupTable, data.frame(modCode = modCodei,
                                             modMass = modMassi))
    }else{
      next
    }
  }

  modLookupTable$modCode <- paste("[", modLookupTable$modCode, "]", sep = "")
  modLookupTable$modMass <- paste("[", modLookupTable$modMass, "]", sep = "")

  tempdf$correctedPep <- tempdf$modifiedPep

  for(i in 1:nrow(modLookupTable)){
    modMassE <- gsub("\\[|\\]", "", modLookupTable$modMass[i])
    tempdf <- tempdf %>%
      dplyr::mutate(correctedPep = ifelse(grepl(modMassE, .data$assignedMods),
                                          gsub(modLookupTable$modCode[i],
                                               modLookupTable$modMass[i],
                                               .data$correctedPep, fixed = TRUE),
                                          .data$correctedPep))
  }

  return(tempdf$correctedPep)
}

GetGlycanMasses <- function(){
  tempdf <- data.frame(monosaccharide = c("N", "H", "A", "F", "G"),
             averageMass = c(203.07937, 162.05282, 291.09542, 146.05791, 307.0903))
  return(tempdf)
}

TTest_log2FC <- function(val1, val2){
  log2Fc <- NA
  pval <- NA

  if(all(is.na(val1)) && all(is.na(val2))){
    return("NA;NA")
  }

  if(sum(is.finite(val1)) > 1 & sum(is.finite(val2)) > 1){
    val1 <- val1[is.finite(val1)]
    val2 <- val2[is.finite(val2)]

    pval <- stats::t.test(val1, val2)$p.value
  }else{
    pval <- NA
  }

  meanval1 <- mean(val1, na.rm=TRUE)
  meanval2 <- mean(val2, na.rm=TRUE)

  if(is.na(meanval1)){
    log2Fc <- -Inf
  }else if(is.na(meanval2)){
    log2Fc <- Inf
  }else{
    log2Fc <- meanval1 -meanval2
  }

  return(paste(log2Fc, pval, sep = ";"))
}

FilterForPeptides <- function(rawdf, whichPeptides){
  if(identical(whichPeptides, NA)){
    return(rawdf)
  }else if(is.data.frame(whichPeptides) && "ModifiedPeptide" %in% names(rawdf) && "ModifiedPeptide" %in% names(whichPeptides)){
    returnVec <- rawdf %>%
      dplyr::filter(.data$ModifiedPeptide %in% whichPeptides$ModifiedPeptide)
  }else if(is.vector(whichPeptides) && length(whichPeptides) > 0){
    returnVec <- rawdf %>%
      dplyr::filter(.data$ModifiedPeptide %in% whichPeptides)
  }else{
    stop("whichPeptides input is invalid")
  }

  if(nrow(returnVec) == 0){
    stop("No data is left after filtering for peptides")
  }else{
    return(returnVec)
  }
}

CheckForQuantitativeValues <- function(intensityValues){
  valid <- !is.na(intensityValues) & intensityValues != 0 & is.finite(intensityValues)

  if (!any(valid)) {
    stop("No quantitative data found. Aborting.")
  }

  invisible(NULL)
}
