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

  tempdf$TotalGlycanComposition <- CleanGlycanComp(tempdf$AssignedModifications,
                                                   tempdf$TotalGlycanComposition)

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

CleanGlycanComp <- function(AssModVec, TotGlycoVec){
  tdf <- data.frame(AssMod = AssModVec,
                       TotGlyco = TotGlycoVec,
                       CorrectGlyco = character(length(TotGlycoVec)))

  for(i in seq_len(nrow(tdf))){
    TotGlycoi <- tdf$TotGlyco[i]
    if(is.na(TotGlycoi) || TotGlycoi == ""){
      tdf$CorrectGlyco[i] <- ""
      next
    }else if(grepl(",", TotGlycoi)){
    TotGlycoSplit <- strsplit(TotGlycoi, ",")[[1]]

    for(j in TotGlycoSplit){
      glycoMass <- ComputeGlycanMass(j)
      AssModMass <- as.double(stringr::str_extract(tdf$AssMod[i], "(?<=\\()[0-9.]+(?=\\))"))

      if(abs(glycoMass - AssModMass) < 0.02){
        tdf$CorrectGlyco[i] <- j
        break}
      }
    }else{
      tdf$CorrectGlyco[i] <- TotGlycoi
    }
    }

  return(tdf$CorrectGlyco)
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

    if(glycanMass %in% .modEnv$ModificationDatabase$ModificationMass){
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
  GlycanDatabaseClean <- GlycanDatabase %>%
    dplyr::mutate(len = nchar(.data$FullName)) %>%
    dplyr::arrange(dplyr::desc("len"))

  GlycanReplacements <- stats::setNames(GlycanDatabaseClean$ShortName, GlycanDatabaseClean$FullName)

  glycan <- stringr::str_replace_all(glycan, GlycanReplacements)
  return(glycan)
}

FilterForCutoffs <- function(input, silent = FALSE){
  if(input$searchEngine %in% c("MSFragger", "Byonic")){
    existingColsPSMTable <- names(input$PSMTable)
    existingColsPTMTable <- names(input$PTMTable)
    if(!silent){
      fmessage(paste0("Prefilter number of rows PSM table: ", nrow(input$PSMTable), ". Prefilter number of rows PTM table: ", nrow(input$PTMTable)))
      if("PSMScore" %in% existingColsPSMTable & !identical(input$peptideScoreCutoff,FALSE)){
        fmessage(paste0("Filtering for PSMScore >= ", input$peptideScoreCutoff))
      }
      if("GlycanQValue" %in% existingColsPSMTable & !identical(input$glycanScoreCutoff,FALSE)){
        fmessage(paste0("Filtering for GlycanQValue <= ", input$glycanScoreCutoff))
      }
      if(input$filterForNoNSequon){
        fmessage(paste0("Filtering for peptides without an N-sequon (OPair peptides only)"))
      }
      if(!identical(input$confidenceLevels,FALSE)){
        fmessage(paste0("Filtering for peptides confidence levels (O-glycopeptides only): ", paste(input$confidenceLevels, collapse= ";")))
      }
      if(!identical(input$deltaModCutoff,FALSE)){
        fmessage(paste0("Filtering for Delta Mod >= ", input$deltaModCutoff))
      }
    }
    #PSMTable####
    if("PSMScore" %in% existingColsPSMTable){
      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$PSMScore >= input$peptideScoreCutoff | is.na(.data$PSMScore))
    }
    if("GlycanQValue" %in% existingColsPSMTable){
      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$GlycanQValue <= input$glycanScoreCutoff | is.na(.data$GlycanQValue))
    }
    if("HasNSequon" %in% existingColsPSMTable && input$filterForNoNSequon){
      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(!.data$HasNSequon | is.na(.data$HasNSequon))
    }

    if("ConfidenceLevel" %in% existingColsPSMTable && !identical(input$confidenceLevels,FALSE)){
      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$ConfidenceLevel %in% input$confidenceLevels | is.na(.data$HasNSequon))
    }

    if("DeltaMod" %in% existingColsPSMTable && !identical(input$deltaModCutoff,FALSE)){
      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$DeltaMod >= input$deltaModCutoff | is.na(.data$DeltaMod))
    }

    #PTMTable####
    if("PSMScore" %in% existingColsPTMTable){
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$PSMScore >= input$peptideScoreCutoff | is.na(.data$PSMScore))
    }
    if("GlycanQValue" %in% existingColsPTMTable){
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$GlycanQValue <= input$glycanScoreCutoff | is.na(.data$GlycanQValue))
    }
    if("HasNSequon" %in% existingColsPTMTable && input$filterForNoNSequon){
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(!.data$HasNSequon | is.na(.data$HasNSequon))
    }

    if("ConfidenceLevel" %in% existingColsPTMTable && !identical(input$confidenceLevels,FALSE)){
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$ConfidenceLevel %in% input$confidenceLevels | is.na(.data$HasNSequon))
    }

    if("DeltaMod" %in% existingColsPTMTable && !identical(input$deltaModCutoff,FALSE)){
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$DeltaMod >= input$deltaModCutoff | is.na(.data$DeltaMod))
    }

    if(!silent){
      fmessage(paste0("Postfilter number of rows PSM table: ", nrow(input$PSMTable), ". Postfilter number of rows PTM table: ", nrow(input$PTMTable)))
    }
    return(input)
  }else if(input$searchEngine %in% c("pGlyco")){
    if(!silent){
      fmessage(paste0("Filtering for PSMScore <= ", input$peptideScoreCutoff, " and glycan score <= ", input$glycanScoreCutoff))
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
  pattern_NX <- "N[^P][ST]"
  pattern_ST <- "[ST]"

  uniprotID <- sub(",.*", "", IDVec[1])

  IDhit <- fastaFile[grepl(uniprotID, names(fastaFile)) & !grepl("rev", names(fastaFile))]
  seq <- toupper(paste(IDhit[[1]], collapse = ""))

  count_NX <- lengths(regmatches(seq, gregexpr(pattern_NX, seq)))
  count_ST <- lengths(regmatches(seq, gregexpr(pattern_ST, seq)))

  rslt <- paste0(count_NX, ";", count_ST)
  return(rslt)
}

fmessage <- function(m){
  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: ", m, "\033[0m")
}

GetRawUniprotInfo <- function(accVec, silent = FALSE){
  #https://www.uniprot.org/help/return_fields
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  fieldsUrl <- "&format=tsv&fields=accession,cc_subcellular_location,ft_intramem,ft_topo_dom,ft_transmem"
  errorIDs <- c()
  scrape_df <- data.frame()

  tempdf <- data.frame(UniprotIDs = accVec) %>%
    dplyr::distinct(.data$UniprotIDs) %>%
    dplyr::mutate(Entry = gsub(",.*", "", .data$UniprotIDs))

  if(!silent){
    fmessage(paste0("Now connecting to Uniprot for ", nrow(tempdf), " proteins...
                    Use scrape = FALSE to the importer to skip this step."))
  }

  totalNum <- nrow(tempdf)
  for (i in seq(1, nrow(tempdf), by = 25)) {
    cat("\rGetting information from Uniprot. Now at protein ", i, "of ", totalNum, "...")
    rows <- tempdf[i:min(i + 25 - 1, nrow(tempdf)), "Entry"]

    fullUrl <- paste0(baseUrl, paste(rows, collapse = "%20OR%20accession:"), fieldsUrl)

    result <- tryCatch({
      result <- suppressWarnings(utils::read.csv(utils::URLencode(fullUrl), header = TRUE, sep = "\t"))
      },
      error = function(e){
        Sys.sleep(0.5)
        res <- httr::GET(fullUrl)
        res <- httr::content(res, as="text", encoding = "UTF-8")

        matches <- stringr::str_extract_all(gsub("'accession'", "", res), "'(.*?)'")[[1]]
        matches <- gsub("'", "", matches)

        errorIDs <- c(errorIDs, matches)

        rows <- rows[!(rows %in% errorIDs)]
        if(length(rows) == 0) next

        fullUrl <- paste0(baseUrl, paste(rows, collapse = "%20OR%20accession:"), fieldsUrl)

        result <- suppressWarnings(utils::read.csv(utils::URLencode(fullUrl), header = TRUE, sep = "\t"))
      })

    scrape_df <- dplyr::bind_rows(scrape_df, result)
    Sys.sleep(0.5)
  }

  scrape_df <- scrape_df %>%
    dplyr::distinct(.data$Entry, .keep_all=TRUE) %>%
    dplyr::right_join(tempdf, by = "Entry") %>%
    dplyr::select(-"Entry")

  return(scrape_df)
}

GetUniprotInfo <- function(UniprotIDs){
  UniprotIDs_df <- data.frame(UniprotIDs = UniprotIDs, stringsAsFactors = FALSE)

  scrape_df <- GetRawUniprotInfo(UniprotIDs_df$UniprotIDs)

  if(nrow(scrape_df) == 0){return(data.frame(
    SubcellularLocalization = rep(NA, nrow(UniprotIDs_df)),
    Domains = rep(NA, nrow(UniprotIDs_df))))}

  scrape_df$SubcellularLocalization <- sapply(scrape_df$Subcellular.location..CC., function(x) cleanSubcellularLocation(x))
  scrape_df$TopoDomain <- sapply(scrape_df$Topological.domain, function(x) getTopoDomain(x))
  scrape_df$TransmembraneDomain <- sapply(scrape_df$Transmembrane, function(x) getTransmembraneDomain(x))

  scrape_df <- scrape_df %>%
    dplyr::mutate(Domains = dplyr::case_when(is.na(.data$TopoDomain) & is.na(.data$TransmembraneDomain) ~ NA_character_,
                                             is.na(.data$TopoDomain) ~ .data$TransmembraneDomain,
                                             is.na(.data$TransmembraneDomain) ~ .data$TopoDomain,
                                             TRUE ~ paste(TopoDomain, TransmembraneDomain, sep = ";"))) %>%
    dplyr::select("UniprotIDs", "SubcellularLocalization", "Domains")

  UniprotIDs_df <- UniprotIDs_df %>%
    dplyr::left_join(scrape_df, by = "UniprotIDs")

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

FilterForProteins <- function(rawdf, whichProtein, exactProteinMatch = TRUE){
  if(identical(whichProtein, NA)){
    return(rawdf)
  }else if(is.data.frame(whichProtein) && "UniprotIDs" %in% names(rawdf) && "UniprotIDs" %in% names(whichProtein)){
    if(exactProteinMatch){
      rawdf <- rawdf %>%
        dplyr::filter(.data$UniprotIDs %in% whichProtein$UniprotIDs)
    }else{
      testList <- sapply(unique(rawdf$UniprotIDs), function(x) strsplit(x, ","))
      check <- sapply(testList, function(x) any(whichProtein$UniprotIDs %in% x))
      toInclude <- names(testList[check])
      rawdf <- rawdf %>%
        dplyr::filter(.data$UniprotIDs %in% toInclude)
    }
  }else if(is.vector(whichProtein) && length(whichProtein) > 0){
    if(exactProteinMatch){
      rawdf <- rawdf %>%
        dplyr::filter(.data$UniprotIDs %in% whichProtein)
    }else{
      testList <- sapply(unique(rawdf$UniprotIDs), function(x) strsplit(x, ","))
      check <- sapply(testList, function(x) any(whichProtein %in% x))
      toInclude <- names(testList[check])
      rawdf <- rawdf %>%
        dplyr::filter(.data$UniprotIDs %in% toInclude)
    }}else{
      stop("whichPeptides input is invalid")
  }

  if(nrow(rawdf) == 0){
    return(rawdf)
  }else{
    return(rawdf)
  }
}

CheckForQuantitativeValues <- function(intensityValues){
  valid <- !is.na(intensityValues) & intensityValues != 0 & is.finite(intensityValues)

  if (!any(valid)) {
    stop("No quantitative data found. Aborting.")
  }

  invisible(NULL)
}

CheckAnnotation <- function(annotation){
  clNames <- names(annotation)
  if(!("Condition") %in% clNames){
    stop("The column 'Condition' is missing")
  }
  if(!("Run") %in% clNames){
    stop("The column 'Run' is missing")
  }
  if(!("BioReplicate") %in% clNames){
    stop("The column 'BioReplicate' is missing")
  }
  if(!("TechReplicate") %in% clNames){
    stop("The column 'TechReplicate' is missing")
  }
  if(anyNA(annotation)){
    warning("Spotted NA values in your annotation file. Did you forget to edit the file?")
  }
  if(anyDuplicated(annotation$Alias) > 0){
    warning("The annotation file contains duplicated values in the Alias column.
            This may result in unexpected behaviour.")
  }
}

ConvertByonicAssignedModifications <- function(ModPepVec, AssignedModsVec){
  # ModPepVec <- "SGLAPNPT[2158.803]N[2158.803]ATTKAAGGGGSGGGSHHHHHHHHHH"
  # AssignedModsVec <- "N[2158.803],T[2158.803]"
  tempdf <- data.frame(ModPep = ModPepVec, AssignedMods = AssignedModsVec)
  tempdf$AssignedMods <- gsub(" ", "", tempdf$AssignedMods)
  tempdf$AssignedMods <- gsub("\\]\\[" , "\\],\\[", tempdf$AssignedMods)
  tempdf$AssignedMods <- gsub("\\+" , "", tempdf$AssignedMods)

  #Convert to exact masses and adds the AA if not listed
  for(i in seq_len(nrow(tempdf))){
    if(is.na(tempdf[i,"AssignedMods"]) | tempdf[i,"AssignedMods"] == ""){
      next
    }
    mods <- strsplit(tempdf[i,"AssignedMods"], ",")[[1]]
    for(j in seq_len(length(mods))){
      modsj <- stringr::str_extract_all(tempdf[i,"ModPep"], "\\[.*?\\]")[[1]]
      tempList <- as.list(as.double(gsub("\\[|\\+|\\]", "", modsj)))
      names(tempList) <- modsj

      modsSub <- mods[j]

      modAssMod <- as.double(stringr::str_extract(modsSub, "(?<=\\[)\\d+\\.?\\d*(?=\\])"))
      indexOfSub <- which.min(abs(unlist(tempList) - modAssMod))

      # Prepend previous AA if modsSub starts with "["
      if(substring(modsSub,1,1) == "[") {
        fullMod <- paste0(substring(mods[j-1], 1, 1), modsSub)
        tempdf[i,"AssignedMods"] <- gsub(paste0(",", modsSub),
                                         paste0(",", fullMod),
                                         tempdf[i,"AssignedMods"],
                                         fixed = TRUE)
      }

      # Replace numeric mass with 3 decimals
      formattedValue <- sprintf("%.3f", tempList[[indexOfSub]][1])
      tempdf[i,"AssignedMods"] <- sub(
        paste0("\\[", modAssMod, "\\]"),
        paste0("[", formattedValue, "]"),
        tempdf[i,"AssignedMods"]
      )
    }
}

  #Add the AA number
  tempdf$formattedAssignedMods <- NA
  for(i in seq_len(nrow(tempdf))){
    if(is.na(tempdf[i,"AssignedMods"]) | tempdf[i,"AssignedMods"] == ""){
      next
    }
    if(grepl("\\*",tempdf[i,"AssignedMods"])){
      rawString <- strsplit(tempdf[i,"AssignedMods"], "\\*")[[1]][1]
      timesA <- as.integer(substring(tempdf[i,"AssignedMods"], nchar(tempdf[i,"AssignedMods"]), nchar(tempdf[i,"AssignedMods"])))
      tempdf[i,"AssignedMods"] <- paste(rep(rawString,timesA), collapse=",")
    }
    mods <- strsplit(tempdf[i,"AssignedMods"], ",")[[1]]

    Pep <- tempdf[i,"ModPep"]
    cleanMod <- c()

    for(j in mods){
      substr <- strsplit(Pep, j, fixed=TRUE)[[1]][1]
      substr <- gsub("[^A-Z]", "", substr)
      loc <- nchar(substr) + 1

      Pep <- sub(j, substring(j,1,1), Pep, fixed=TRUE)

      cleanMod <- c(cleanMod, paste0(loc,
                                     substring(j,1,1),
                                     "[",
                                     substring(j,3,nchar(j))))
    }
    tempdf$formattedAssignedMods[i] <- paste(cleanMod, collapse=",")
  }

  #Clean and return
  tempdf$formattedAssignedMods <- gsub("\\[", "\\(", tempdf$formattedAssignedMods)
  tempdf$formattedAssignedMods <- gsub("\\]", "\\)", tempdf$formattedAssignedMods)

  returnVal <- as.vector(tempdf$formattedAssignedMods)

  return(returnVal)
}

ComputeGlycanMass <- function(glycanComposition){
  if(is.na(glycanComposition) || glycanComposition == ""){
    return(NA)
  }
  glycanMass_df <- GlycanDatabase %>%
    dplyr::mutate(count = 0, mass = 0)

  for(i in seq_len(nrow(glycanMass_df))){
    mono <- glycanMass_df$ShortName[i]

    pattern <- paste0("(?<![A-Za-z])", mono, "\\(([0-9]+)\\)")
    match <- stringr::str_match(glycanComposition, pattern)[,2]

    count <- as.numeric(match)

    if(!is.na(count)){
      glycanMass_df$count[i] <- count
      glycanMass_df$mass[i] <- count * glycanMass_df$GMass[i]
    }
  }

  totalMass <- sum(glycanMass_df$mass, na.rm = TRUE)

  return(totalMass)
}

AssignedModsToGlycanComp_Byonic <- function(AssModVec, GlycanTable){
  AssMod <- data.frame(AssMod = AssModVec,
                          TotGlyComp = character(length(AssModVec)))
  TotalGlycanComp <- c()
  for(i in seq_len(nrow(AssMod))){
    AssModVeci <- AssMod$AssMod[i]

    TotGlycoi <- c()

    if(AssModVeci == "" || is.na(AssModVeci)){
      AssMod$TotGlyComp[i] <- NA
      next
    }

    AssModSplit <- strsplit(AssModVeci, ",")[[1]]

    for(j in AssModSplit){
      modMass <- as.double(stringr::str_extract(j, "(?<=\\()[0-9.]+(?=\\))"))

      GlycCompj <- GlycanTable %>%
        dplyr::filter(abs(.data$Mass - modMass) < 0.02)

      if(nrow(GlycCompj) > 0){
        GlycCompj <- GlycCompj %>%
          dplyr::slice_min(order_by = abs(.data$Mass - modMass), n = 1) %>%
          dplyr::pull("Rule")

        TotGlycoi <- c(TotGlycoi, GlycCompj)
      }else{
        TotGlycoi <- c(TotGlycoi, "")
      }
    }
    AssMod$TotGlyComp[i] <- gsub(" ", "", paste(TotGlycoi[TotGlycoi != ""], collapse = ","))
  }
  return(AssMod$TotGlyComp)
}

GetModifiedPeppGlyco <- function(PepVec, ModVec, ModDB){
  tempdf <- data.frame(Pep = PepVec,
                       Mod = ModVec,
                       cleanName = character(length(PepVec)))

  for(i in seq_len(nrow(tempdf))){
    modi <- tempdf[i,"Mod"]
    if(modi == "" | is.na(modi)){
      tempdf$cleanName[i] <- tempdf[i,"Pep"]
    }else{
      modi <- strsplit(modi, ";")[[1]]
      pepi <- tempdf[i,"Pep"]
      addAA <- 0

      num <- as.numeric(stringr::str_extract(modi, "^[0-9]+"))

      modi <- modi[order(num)]

      for(j in seq_len(length(modi))){
        residueNumber <- as.numeric(stringr::str_extract(modi[j], "^[0-9]+"))
        ModName <- stringr::str_extract(modi[j], "(?<=,).*?(?=\\[)")
        ModMass <- .modEnv$ModificationDatabase %>%
          dplyr::filter(.data$FullName == ModName) %>%
          dplyr::slice(1) %>%
          dplyr::pull(.data$ModificationMass)
        massLength <- nchar(ModMass) + 2

        #print(paste0(modi[j], "---", residueNumber, "---", ModName, "---", ModMass, "---", massLength))

        substr1 <- substr(pepi, 1, residueNumber + addAA)
        substr2 <- substr(pepi, residueNumber + addAA + 1, nchar(pepi))

        pepi <- paste0(substr1, "[", ModMass, "]", substr2)
        addAA <- addAA + massLength
      }
      tempdf$cleanName[i] <- pepi
    }
  }

  tempdf$cleanName <- gsub("\\[57\\.0215\\]", "", tempdf$cleanName)
  tempdf$cleanName <- gsub("J", "N", tempdf$cleanName)
  return(tempdf$cleanName)
}

AssignmedModspGlyco <- function(modVec, glymassVec, glysiteVec,
                                peptideVec, modDB){
  tempdf <- data.frame(mod = modVec, glymass = glymassVec,
                       glysite = glysiteVec, formattedMod = character(length(modVec)),
                       peptide = gsub("J", "N", peptideVec),
                       formattedGlycan = character(length(modVec)),
                       formattedAssignedMod = character(length(modVec)))

  #First clean mod column
  for(i in seq_len(nrow(tempdf))){
    modi <- tempdf[i,"mod"]
    if(modi == "" | is.na(modi)){
      tempdf$formattedMod[i] <- ""
    }else{
      modi <- strsplit(modi, ";")[[1]]
      cleanedModi <- c()

      for(j in seq_len(length(modi))){
        ResidueNumber <- as.numeric(stringr::str_extract(modi[j], "^[0-9]+"))
        ModName <- stringr::str_extract(modi[j], "(?<=,).*?(?=\\[)")
        ModMass <- .modEnv$ModificationDatabase %>%
          dplyr::filter(.data$FullName == ModName) %>%
          dplyr::slice(1) %>%
          dplyr::pull(.data$ModificationMass)
        ModResidue <- stringr::str_extract(modi[j], "(?<=\\[)[^]]+(?=\\])")

        cleanedModi <- c(cleanedModi, paste0(ResidueNumber, ModResidue, "(", ModMass, ")"))
      }

      tempdf$formattedMod[i] <- paste(cleanedModi, collapse=",")
  }
  }

  #Add the glycan part
  for(i in seq_len(nrow(tempdf))){
    glymassi <- tempdf[i,"glymass"]
    if(glymassi == "" | is.na(glymassi)){
      tempdf$formattedMod[i] <- ""
    }else{
      glymassi <- strsplit(as.character(glymassi), ";")[[1]]
      glysitei <- strsplit(as.character(tempdf[i,"glysite"]), ";")[[1]]
      pepi <- tempdf[i,"peptide"]

      cleanGlycan <- c()

      for(j in seq_len(length(glymassi))){
        ResidueNumber <- glysitei[j]
        ModMass <- glymassi[j]
        ModResidue <- substring(pepi, as.numeric(ResidueNumber), as.numeric(ResidueNumber))

        cleanGlycan <- c(cleanGlycan, paste0(ResidueNumber, ModResidue, "(", ModMass, ")"))
      }
      tempdf$formattedGlycan[i] <- paste(cleanGlycan, collapse=",")
    }
    }

    #Stitch them together and get them in the correct order
  tempdf$formattedAssignedMod <- paste(tempdf$formattedMod, tempdf$formattedGlycan, sep = ",")
  tempdf$formattedAssignedMod <- sub("^,+", "", tempdf$formattedAssignedMod)

  for(i in seq_len(nrow(tempdf))){
    modi <- tempdf[i,"formattedAssignedMod"]
    if(modi == "" | is.na(modi)){
      tempdf$formattedAssignedMod[i] <- tempdf[i,"formattedAssignedMod"]
    }else{
      modi <- strsplit(modi, ",")[[1]]
      num <- as.numeric(stringr::str_extract(modi, "^[0-9]+"))
      modi <- modi[order(num)]
      tempdf$formattedAssignedMod[i] <- paste(modi, collapse = ",")
    }
    }

  return(tempdf$formattedAssignedMod)
}

GetPeptideLocInProtein <- function(uniprotID, pep, fastaFile){
  uniprotID <- strsplit(uniprotID[1], ",")[[1]][1]
  IDhit <- fastaFile[grepl(uniprotID, names(fastaFile)) & !grepl("rev", names(fastaFile))]
  seq <- paste(IDhit[[1]], collapse = "")
  seq <- toupper(seq)

  uppercase_only <- gsub("[^A-Z]", "", pep[1])

  AANumber <- regexpr(uppercase_only, seq)

  return(AANumber)
}

UpdateFPIntensities <- function(rawdata, quantdata, normalization){
  uniqueRundf <- data.frame(Run = unique(rawdata$Run))

  if(normalization == "FP_Normalized"){
    uniqueRundf$colName <- paste("X", uniqueRundf$Run, ".Intensity", sep = "")
  }else if(normalization == "FP_MaxLFQ"){
    uniqueRundf$colName <- paste("X", uniqueRundf$Run, ".MaxLFQ.Intensity", sep = "")
  }else{
    stop("The normalization was not recognized: ", normalization)
  }

  if(any(!uniqueRundf$colName %in% names(quantdata))){
    notfound <- uniqueRundf$Run[!uniqueRundf$colName %in% names(quantdata)]

    warning(
      "Did not find the following runs in combined_peptide.tsv:\n",
      paste(" -", notfound, collapse = "\n")
    )
  }

  quantdata <- quantdata %>%
    dplyr::select("ModifiedPeptide" = "Modified.Sequence",
                  dplyr::any_of(uniqueRundf$colName)) %>%
    dplyr::mutate(ModifiedPeptide = gsub("\\[57\\.0214\\]|\\[57\\.0215\\]", "", .data$ModifiedPeptide)) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(uniqueRundf$colName), names_to = "colName", values_to = "Intensity") %>%
    dplyr::filter(!is.na(.data$Intensity))%>%
    dplyr::left_join(uniqueRundf, by = "colName") %>%
    dplyr::select(-"colName")

  rawdata <- rawdata %>%
    dplyr::select(-"Intensity") %>%
    dplyr::full_join(quantdata, by = c("Run", "ModifiedPeptide"))

  #Fill the rawdata data with the rest of the data
  rawdata <- rawdata %>%
    dplyr::group_by(.data$ModifiedPeptide) %>%
    tidyr::fill(dplyr::any_of(c("ModifiedPeptide", "AssignedModifications", "TotalGlycanComposition",
                  "IsUnique", "UniprotIDs", "Genes", "ProteinLength", "ConfidenceLevel", "NumberOfNSites",
                  "NumberOfOSites", "ProteinStart", "GlycanType", "SubcellularLocalization",
                  "Domains", "RetentionTime", "ID", "GlycanQValue", "PSMScore")),
                .direction = "downup") %>%
    dplyr::ungroup()

  if(any(is.na(rawdata$Intensity))){
    warning("NA intensity values detected after importing combined_peptide.tsv")
  }

  rawdata <- rawdata %>%
    dplyr::mutate(Intensity = dplyr::case_when(is.na(.data$RawIntensity) & .data$Intensity == 0 ~ NA_real_,
                                               TRUE ~ .data$Intensity)) %>%
    dplyr::filter(!(is.na(.data$RawIntensity) & !is.na(.data$Intensity)))

  return(rawdata)
}

Databases <- function(){
  GlycanDatabase <- data.frame(
    FullName = c("HexNAc", "Hex", "NeuAc", "Fuc", "NeuGc", "Pent",
                 "KDN", "HexA", "pseudaminic"),
    ShortName = c("N", "H", "A", "F", "G", "P", "Kdn", "HexA", "p"),
    GMass = c(203.07937, 162.05282, 291.09542, 146.05791, 307.09033, 132.0423,
              250.06889, 176.03209, 232.10592),
    stringsAsFactors = FALSE
  )

  ModificationDatabase <- data.frame(
    FullName = c("Oxidation", "CCarbamidomethylation1", "CCarbamidomethylation2",
                 "NAcetylation", "Carbamidomethyl"),
    ModificationMass = c("15.9949", "57.0214", "57.0215", "42.0106", "57.0215")
  )

  GlycanColors = data.frame(GlycanType = c("Complex/Hybrid", "Sialyl+Fucose", "Sialyl",
                                           "Fucose", "High Mannose", "Truncated",
                                           "Paucimannose", "OGlycan", "NonCanonicalGlyco"),
                            color = c("#00394a", "#ff7f2a", "#2475b5", "#aaaaaa", "#28b36d",
                                      "#D0A5C0", "#8B1E3F", "#f2d46f", "#6a4c8b"))

  colorScheme <- c(
    "#BAA5CC", "#9ADCEE", "#BAD97C", "#EEAED0", "#FAD821", "#94D8C3", "#F7B8D2", "#A7C7E7",
    "#FFE87C", "#C0E4D0", "#A1A9F2", "#C1D87F", "#E3B7E2", "#B1D3C2", "#F9A9B6", "#D1D2E3",
    "#A4EFA1", "#D9D07A", "#98C9C7", "#F4D1A1"
  )

  #usethis::use_data(GlycanDatabase, colorScheme, ModificationDatabase, GlycanColors, internal = TRUE, overwrite = TRUE)
}
