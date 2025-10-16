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

  # tempdf <- tempdf %>%
  #     dplyr::rowwise() %>%
  #     dplyr::mutate(PeptidePTMLocalization = as.numeric(regmatches(AssignedModifications, gregexpr("[0-9]+", AssignedModifications))[[1]][1]),
  #                   ProteinPTMLocalization = PeptidePTMLocalization + ProteinStart,
  #                   ModificationSite = sub(".*([A-Za-z])\\(.*", "\\1", AssignedModifications),
  #                   ModificationID = paste0(ModificationSite,ProteinPTMLocalization))

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
      }else if(!(modifiedResidue %in% c("S", "T", "N"))){
        modType <- append(modType, "NonGlyco")
        }else if((glycanComp != "" & modifiedResidue == "N") | (!is.na(glycanComp) & modifiedResidue == "N")){
      glycanMass = strsplit(glycanComp, "%")[[1]][2]
      glycanMass = substring(glycanMass, 2, nchar(glycanMass) - 1 )

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
      }else if((glycanComp != "" & modifiedResidue %in% c("S", "T")) | (!is.na(glycanComp) & modifiedResidue %in% c("S", "T"))){
        modType <- append(modType, "OGlycan")
      }else{
        modType <- append(modType, "NonGlyco")
      }
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
      dplyr::mutate(.by = c(.data$ModifiedPeptide, .data$AssignedModifications,
                            .data$Condition, .data$BioReplicate),
                    Intensity = stats::median(.data$Intensity, na.rm = TRUE)) %>%
      dplyr::distinct(.data$ModifiedPeptide, .data$Condition, .data$BioReplicate,
                      .data$ModificationID, .keep_all = TRUE)

    return(df)
  }else{
    df <- df %>%
      dplyr::mutate(.by = c(.data$ModifiedPeptide, .data$AssignedModifications,
                            .data$Condition, .data$BioReplicate),
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
  localMedian <- median(intensityVec[intensityVec != 0], na.rm = TRUE)
  deltaMedian <- globalMedian - localMedian
  intensityVec[intensityVec != 0] <- intensityVec[intensityVec != 0] + deltaMedian
  return(intensityVec)
}
