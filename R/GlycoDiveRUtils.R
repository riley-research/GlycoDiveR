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
  ids <- sub(",.*", "", IDVec)

  seqs <- fastaFile[ids]

  return(as.vector(nchar(seqs)))
}

PSMToPTMTable <- function(PSMTable, silent = FALSE){
  tempdf <- PSMTable %>%
    dplyr::filter(!is.na(.data$AssignedModifications) & .data$AssignedModifications != "") %>%
    tidyr::separate_rows(.data$AssignedModifications, sep = ",")

  tempdf <- tempdf %>%
    dplyr::mutate(AssignedModifications =
                    dplyr::if_else(grepl("N-term", .data$AssignedModifications),
                                   paste0("1", stringr::str_extract(.data$ModifiedPeptide, "[A-Z]"), gsub("N-term", "", .data$AssignedModifications)),
                                   .data$AssignedModifications))

  unique_mods <- tempdf %>%
    dplyr::distinct(.data$AssignedModifications, .data$TotalGlycanComposition)

  unique_mods <- unique_mods %>%
    dplyr::mutate(
      TotalGlycanCompositionCorrect = mapply(
        CleanGlycanComp,
        AssModVec = .data$AssignedModifications,
        TotGlycoVec = .data$TotalGlycanComposition
      ) %>% as.character()
    )

  tempdf <- tempdf %>%
    dplyr::left_join(unique_mods, by = c("AssignedModifications", "TotalGlycanComposition")) %>%
    dplyr::select(-"TotalGlycanComposition", "TotalGlycanComposition" = "TotalGlycanCompositionCorrect")

  # tempdf$TotalGlycanComposition <- CleanGlycanComp(tempdf$AssignedModifications,
  #                                                  tempdf$TotalGlycanComposition)

  tempdf <- tempdf %>%
    dplyr::mutate(
      PeptidePTMLocalization = as.numeric(stringr::str_extract(.data$AssignedModifications, "\\d+")),
      ProteinPTMLocalization = .data$PeptidePTMLocalization + .data$ProteinStart -1,
      ModificationSite = stringr::str_extract(.data$AssignedModifications, "[A-Za-z](?=\\()"),
      ModificationID = paste0(.data$ModificationSite, .data$ProteinPTMLocalization))

  unique_mods <- tempdf %>%
    dplyr::distinct(.data$AssignedModifications, .data$TotalGlycanComposition)

  unique_mods <- unique_mods %>%
    dplyr::mutate(
      GlycanTypeCorrect = mapply(
        GlycanComptToGlycanType,
        mod = .data$AssignedModifications,
        glycanComp = .data$TotalGlycanComposition
      ) %>% as.character()
    )

  tempdf <- tempdf %>%
    dplyr::left_join(unique_mods, by = c("AssignedModifications", "TotalGlycanComposition"))  %>%
    dplyr::select(-"GlycanType", "GlycanType" = "GlycanTypeCorrect") %>%
    dplyr::mutate(TotalGlycanComposition = ifelse(.data$TotalGlycanComposition == "", NA, .data$TotalGlycanComposition))

  # tempdf$GlycanType <- apply(tempdf[,c("AssignedModifications", "TotalGlycanComposition")],
  #                            1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))

  if(!silent){
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Generated PTM table.\033[0m")
  }

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

GlycanComptToGlycanType_legacy <- function(mod, glycanComp){
  modType <- c()
  allModsVec <- if (grepl(",", mod, fixed = TRUE)) strsplit(mod, ",")[[1]] else mod

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

        hexNAc_count[is.na(hexNAc_count)] <- 0
        hex_count[is.na(hex_count)] <- 0

        glycanCat <- dplyr::case_when(
          grepl("A|G", glycanComp) & grepl("F", glycanComp) ~ "Sialofucosylated",
          grepl("A|G", glycanComp) ~ "Sialylated",
          grepl("F", glycanComp) ~ "Fucosylated",
          grepl("Phospho", glycanComp) ~ "Phosphomannose",
          !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count == 2 & hex_count == 3 ~ "Truncated",
          (!is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 2) |
            ( !is.na(hexNAc_count) & !is.na(hex_count) & hex_count < 3) ~ "Truncated",
          !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 3 & hex_count > 3 ~ "Oligomannose",
          TRUE ~ "Complex/Hybrid"
        )
        modType <- append(modType, glycanCat)
      }else{
        modType <- append(modType, "UndefinedGlyco")
      }}else if((glycanComp != "" & modifiedResidue %in% c("S", "T")) | (!is.na(glycanComp) & modifiedResidue %in% c("S", "T"))){
        if(.modEnv$useExtendedOGlycanCategories){
          glycanCat <- dplyr::case_when(
            grepl("A|G", glycanComp) & grepl("F", glycanComp) ~ "Sialofucosylated",
            grepl("A|G", glycanComp) ~ "Sialylated",
            grepl("F", glycanComp) ~ "Fucosylated",
            TRUE ~ "OGlycan"
          )
          modType <- append(modType, glycanCat)
        }else{
          modType <- append(modType, "OGlycan")
        }
      }else if(glycanComp != "" & !is.na(glycanComp)){
        modType <- append(modType, "NonCanonicalGlyco")
      }else{
        modType <- append(modType, "NonGlyco")
      }
  }

  return(paste(modType, collapse = ", "))
}

GlycanComptToGlycanType <- function(mod, glycanComp) {
  rule_df <- .modEnv$GlycanCategories

  if (is.na(mod) || mod == "") return("Unmodified")

  allModsVec <- if (grepl(",", mod, fixed = TRUE)) strsplit(mod, ",")[[1]] else mod

  h_count <- as.numeric(stringr::str_match(glycanComp, "H\\((\\d+)\\)")[,2]) %>% tidyr::replace_na(0)
  n_count <- as.numeric(stringr::str_match(glycanComp, "N\\((\\d+)\\)")[,2]) %>% tidyr::replace_na(0)

  modType <- sapply(allModsVec, function(i) {
    res_part <- sub("\\(.*", "", i)
    modifiedResidue <- substr(res_part, nchar(res_part), nchar(res_part))

    glycanMass <- sub(".*\\((.*)\\).*", "\\1", i)

    if (glycanMass %in% .modEnv$ModificationDatabase$ModificationMass) {
      return("NonGlyco")
    }

    if (!is.na(glycanComp) && glycanComp != "") {

      match <- rule_df %>%
        dplyr::filter(
          stringr::str_detect(modifiedResidue, .data$Residue),
          stringr::str_detect(glycanComp, .data$Pattern),
          h_count >= .data$Min_H & h_count <= .data$Max_H,
          n_count >= .data$Min_N & n_count <= .data$Max_N
        ) %>%
        dplyr::arrange(.data$Priority) %>%
        dplyr::slice(1)

      if (nrow(match) > 0) {
        if (modifiedResidue %in% c("S", "T") && !.modEnv$useExtendedOGlycanCategories) {
          return("OGlycan")
        }
        return(match$Category)
      }

      if (modifiedResidue %in% c("S", "T")) return("OGlycan")
      return("NonCanonicalGlyco")

    } else {
      return("NonGlyco")
    }
  })

  return(paste(modType, collapse = ", "))
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

  if(nrow(df) < 1){return(df)}

  if("ModificationID" %in% names(df)){
    #Take median of technical reps together
    alias_levels <- levels(df$Alias)

    df <- df %>%
      dplyr::mutate(.by = c("Condition", "BioReplicate"),
                    Alias = min(as.character(.data$Alias), na.rm = TRUE) |> factor(),
                    TechReplicate = min(.data$TechReplicate, na.rm = TRUE)) %>%
      dplyr::mutate(Alias = factor(.data$Alias, levels = alias_levels)) %>%
      dplyr::mutate(.by = c("ModifiedPeptide", "AssignedModifications",
                            "Condition", "BioReplicate"),
                    Intensity = stats::median(.data$Intensity, na.rm = TRUE)) %>%
      dplyr::distinct(.data$ModifiedPeptide, .data$Condition, .data$BioReplicate,
                      .data$ModificationID, .keep_all = TRUE)
    df$Alias <- droplevels(df$Alias)

    return(df)
  }else{
     df <- df %>%
       dplyr::mutate(.by = c("Condition", "BioReplicate"),
                     Alias = min(as.character(.data$Alias), na.rm = TRUE) |> factor(),
                     TechReplicate = min(.data$TechReplicate, na.rm = TRUE)) %>%
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
  existingColsPSMTable <- names(input$PSMTable)
  existingColsPTMTable <- names(input$PTMTable)
  preFilterPSMRows <- nrow(input$PSMTable)
  preFilterPTMRows <- nrow(input$PTMTable)
  if(input$searchEngine %in% c("MSFragger", "Byonic")){
    if(!silent){
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
  }else if(input$searchEngine %in% c("pGlyco")){
    if(!silent){
      fmessage(paste0("Filtering for PSMScore <= ", input$peptideScoreCutoff, " and glycan score <= ", input$glycanScoreCutoff))
    }
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$PSMScore <= input$peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue <= input$glycanScoreCutoff | is.na(.data$GlycanQValue))

    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$PSMScore <= input$peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue <= input$glycanScoreCutoff | is.na(.data$GlycanQValue))

  }else if(input$searchEngine %in% c("GlycanFinder")){
    if(!silent){
      fmessage(paste0("Filtering for PSMScore >= ", input$peptideScoreCutoff, " and glycan score >= ", input$glycanScoreCutoff))
    }

    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$PSMScore >= input$peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue >= input$glycanScoreCutoff | is.na(.data$GlycanQValue))

    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$PSMScore >= input$peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue >= input$glycanScoreCutoff | is.na(.data$GlycanQValue))

    if("SScore" %in% existingColsPSMTable & !identical(input$SScoreCutoff,FALSE)){
      if(!silent){fmessage(paste0("Filtering for S Score >= ", input$SScoreCutoff))}

      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$SScore >= input$SScoreCutoff | is.na(.data$SScore))
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$SScore >= input$SScoreCutoff | is.na(.data$SScore))
    }
    if("StructureConfidence" %in% existingColsPSMTable & !identical(input$structureConfidenceCutoff,FALSE)){
      if(!silent){fmessage(paste0("Filtering for structure Confidence == ", input$structureConfidenceCutoff))}

      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$StructureConfidence %in% input$structureConfidenceCutoff | is.na(.data$StructureConfidence))
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$StructureConfidence %in% input$structureConfidenceCutoff | is.na(.data$StructureConfidence))
    }
    if("AScore" %in% existingColsPSMTable & !identical(input$AScoreCutoff,FALSE)){
      if(!silent){fmessage(paste0("Filtering for AScore >= ", input$AScoreCutoff))}

      input$PSMTable <- input$PSMTable %>%
        dplyr::filter(.data$AScore >= input$AScoreCutoff | is.na(.data$AScore))
      input$PTMTable <- input$PTMTable %>%
        dplyr::filter(.data$AScore >= input$AScoreCutoff | is.na(.data$AScore))
    }

  }else{
    warning("No search engine recognized, returning unfiltered dataframe.")
    return(input)
  }

  if(!silent){
    fmessage(paste0("Filtered PSM table: ", preFilterPSMRows, " to ", nrow(input$PSMTable), " rows and ",
                    "PTM table: ", preFilterPTMRows, " to ", nrow(input$PTMTable), " rows."))
  }
  return(input)
}

FilterPSMTable <- function(df,
                           peptideScoreCutoff,
                           glycanScoreCutoff,
                           filterForNoNSequon,
                           confidenceLevels,
                           deltaModCutoff,
                           searchEngine,
                           SScoreCutoff,
                           structureConfidenceCutoff,
                           AScoreCutoff){
  existingColsPSMTable <- names(df)
  preFilterPSMRows <- nrow(df)

  if(searchEngine %in% c("MSFragger", "Byonic")){
      if("PSMScore" %in% existingColsPSMTable & !identical(peptideScoreCutoff,FALSE)){
        fmessage(paste0("Filtering for PSMScore >= ", peptideScoreCutoff))
      }
      if("GlycanQValue" %in% existingColsPSMTable & !identical(glycanScoreCutoff,FALSE)){
        fmessage(paste0("Filtering for GlycanQValue <= ", glycanScoreCutoff))
      }
      if(filterForNoNSequon){
        fmessage(paste0("Filtering for peptides without an N-sequon (OPair peptides only)"))
      }
      if(!identical(confidenceLevels,FALSE)){
        fmessage(paste0("Filtering for peptides confidence levels (O-glycopeptides only): ", paste(confidenceLevels, collapse= ";")))
      }
      if(!identical(deltaModCutoff,FALSE)){
        fmessage(paste0("Filtering for Delta Mod >= ", deltaModCutoff))
      }
    if("PSMScore" %in% existingColsPSMTable){
      df <- df %>%
        dplyr::filter(.data$PSMScore >= peptideScoreCutoff | is.na(.data$PSMScore))
    }
    if("GlycanQValue" %in% existingColsPSMTable){
      df <- df %>%
        dplyr::filter(.data$GlycanQValue <= glycanScoreCutoff | is.na(.data$GlycanQValue))
    }
    if("HasNSequon" %in% existingColsPSMTable && filterForNoNSequon){
      df <- df %>%
        dplyr::filter(!.data$HasNSequon | is.na(.data$HasNSequon))
    }

    if("ConfidenceLevel" %in% existingColsPSMTable && !identical(confidenceLevels,FALSE)){
      df <- df %>%
        dplyr::filter(.data$ConfidenceLevel %in% confidenceLevels | is.na(.data$HasNSequon))
    }

    if("DeltaMod" %in% existingColsPSMTable && !identical(deltaModCutoff,FALSE)){
      df <- df %>%
        dplyr::filter(.data$DeltaMod >= deltaModCutoff | is.na(.data$DeltaMod))
    }
  }else if(searchEngine %in% c("pGlyco")){

    fmessage(paste0("Filtering for PSMScore <= ", peptideScoreCutoff, " and glycan score <= ", glycanScoreCutoff))

    df <- df %>%
      dplyr::filter(.data$PSMScore <= peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue <= glycanScoreCutoff | is.na(.data$GlycanQValue))
  }else if(searchEngine %in% c("pGlyco")){

    fmessage(paste0("Filtering for PSMScore >= ", peptideScoreCutoff, " and glycan score >= ", glycanScoreCutoff))

    df <- df %>%
      dplyr::filter(.data$PSMScore >= peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue >= glycanScoreCutoff | is.na(.data$GlycanQValue))
  }else if(searchEngine %in% c("GlycanFinder")){
    fmessage(paste0("Filtering for PSMScore >= ", peptideScoreCutoff, " and glycan score >= ", glycanScoreCutoff))

    df <- df %>%
      dplyr::filter(.data$PSMScore >= peptideScoreCutoff | is.na(.data$PSMScore)) %>%
      dplyr::filter(.data$GlycanQValue >= glycanScoreCutoff | is.na(.data$GlycanQValue))
  }

  fmessage(paste0("Filtered PSM table: ", preFilterPSMRows, " to ", nrow(df), " rows."))
  if("SScore" %in% existingColsPSMTable & !identical(SScoreCutoff,FALSE)){
    fmessage(paste0("Filtering for S Score >= ", SScoreCutoff))

    df <- df %>%
      dplyr::filter(.data$SScore >= SScoreCutoff | is.na(.data$SScore))
  }
  if("StructureConfidence" %in% existingColsPSMTable & !identical(structureConfidenceCutoff,FALSE)){
    fmessage(paste0("Filtering for structure Confidence == ", structureConfidenceCutoff))

    df <- df %>%
      dplyr::filter(.data$StructureConfidence %in% structureConfidenceCutoff | is.na(.data$StructureConfidence))
  }
  if("AScore" %in% existingColsPSMTable & !identical(AScoreCutoff,FALSE)){
    fmessage(paste0("Filtering for AScore >= ", AScoreCutoff))

    df <- df %>%
      dplyr::filter(.data$AScore >= AScoreCutoff | is.na(.data$AScore))
  }

  return(df)
}

GetGlycoSitesPerProtein <- function(IDVec, fastaFile){
  ids <- sub(",.*", "", IDVec)

  seqs <- fastaFile[ids]

  count_NX <- stringr::str_count(seqs, "N[^P][ST]")
  count_ST <- stringr::str_count(seqs, "[ST]")

  rslt <- paste0(count_NX, ";", count_ST)
  return(rslt)
}

fmessage <- function(m){
  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: ", m, "\033[0m")
}

GetRawUniprotInfo <- function(accVec, size = 100, silent = FALSE){
  #https://www.uniprot.org/help/return_fields
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  fieldsUrl <- "&format=tsv&fields=accession,cc_subcellular_location,ft_intramem,ft_topo_dom,ft_transmem"
  sizeUrl <- paste0("&size=", size)
  errorIDs <- c()
  scrape_df <- data.frame()

  tempdf <- data.frame(UniprotIDs = accVec) %>%
    dplyr::distinct(.data$UniprotIDs) %>%
    dplyr::mutate(Entry = gsub(",.*", "", .data$UniprotIDs))

  if(!silent){
    fmessage(paste0("Now connecting to Uniprot for ", nrow(tempdf), " proteins...
                    Use scrape = FALSE using the importer function to skip this step."))
  }

  totalNum <- nrow(tempdf)
  for (i in seq(1, nrow(tempdf), by = size)) {
    cat("\rGetting information from Uniprot. Now at protein ", i, "of ", totalNum, "...")
    rows <- tempdf[i:min(i + size - 1, nrow(tempdf)), "Entry"]

    fullUrl <- paste0(baseUrl, paste(rows, collapse = "%20OR%20accession:"), fieldsUrl, sizeUrl)

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
        if(length(rows) == 0) {return(data.frame())}

        fullUrl <- paste0(baseUrl, paste(rows, collapse = "%20OR%20accession:"), fieldsUrl, sizeUrl)

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
    if(!grepl("..", rawStringVec[i], fixed=TRUE)){
      next
    }
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
  output <- intensityVec
  intensityVec_log2 <- log2(intensityVec[nonzero])
  localMedian <- stats::median(intensityVec_log2, na.rm = TRUE)
  deltaMedian <- log2(globalMedian) - localMedian
  output[nonzero] <- 2^(intensityVec_log2 + deltaMedian)
  return(output)
}

FPModCodeToModMass <- function(modifiedPep, assignedMods){
  #1. Generate a dataframe with the modification codes and modification masses
  #2. Look over to rows to replace the codes with the masses
  # FPModCodeToModMass(modifiedPep = mydata$PSMTable$ModifiedPeptide, assignedMods = mydata$PSMTable$AssignedModifications)
  modLookupTable <- data.frame(modCode = as.character(),
                               modMass = as.character())

  tempdf <- data.frame(modifiedPep = modifiedPep, assignedMods = assignedMods)

  modTable <- tempdf %>%
    dplyr::filter(!is.na(.data$assignedMods) & .data$assignedMods != "") %>%
    tidyr::separate_longer_delim(cols = assignedMods, delim = ",") %>%
    dplyr::filter(!grepl("C\\(57\\.02", .data$assignedMods)) %>%
    dplyr::distinct()

  modTable$Mass <- modMass <- stringr::str_extract(modTable$assignedMods, "(?<=\\()[^)]*(?=\\))")

  is_nterm <- stringr::str_detect(modTable$assignedMods, "N-term")

  modTable$splitIdx <- ifelse(
    is_nterm,
    3,
    as.numeric(stringr::str_extract(modTable$assignedMods, "^[0-9]+"))
  )

  letterPos <- gregexpr("[A-Z]", modTable$modifiedPep)

  modTable$splitAA <- mapply(function(posVec, idx) {
    if (is.na(idx) || idx < 1 || idx > length(posVec)) return(NA_integer_)
    posVec[idx] + 2
  }, letterPos, modTable$splitIdx)
  modTable <- modTable %>%
    dplyr::mutate(
      splitAA = dplyr::case_when(
        stringr::str_detect(.data$assignedMods, "N-term") ~ 3L,
        TRUE ~ as.integer(.data$splitAA)))

  modTable$getString <- substr(modTable$modifiedPep, modTable$splitAA, nchar(modTable$modifiedPep))
  modTable$modCode <- stringr::str_extract(modTable$getString, "^[0-9]+")

  modTable <- modTable %>%
    dplyr::mutate(modCode = dplyr::if_else(
      stringr::str_detect(.data$assignedMods, "C-term"),
      stringr::str_extract(.data$modifiedPep, "(?<=\\[)[0-9]+(?=\\][^\\[]*$)"),
      .data$modCode)) %>%
    dplyr::select("modMass" = "Mass", "modCode") %>%
    dplyr::distinct(.data$modMass, .data$modCode, .keep_all = TRUE)

  modTable$modCode <- paste("[", modTable$modCode, "]", sep = "")
  modTable$modMass <- paste("[", modTable$modMass, "]", sep = "")

  tempdf$correctedPep <- tempdf$modifiedPep
  modTable$modMassClean <- gsub("\\[|\\]", "", modTable$modMass)

  for (k in seq_len(nrow(modTable))) {

    hit <- grepl(modTable$modMassClean[k], tempdf$assignedMods)

    if (any(hit)) {
      tempdf$correctedPep[hit] <- gsub(
        modTable$modCode[k],
        modTable$modMass[k],
        tempdf$correctedPep[hit],
        fixed = TRUE
      )
    }
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
  if(identical(whichPeptides, NULL)){
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
    return(returnVec)
  }else{
    return(returnVec)
  }
}

FilterForProteins <- function(rawdf, whichProtein, exactProteinMatch = TRUE){
  if(identical(whichProtein, NULL)){
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
    return(TRUE)
  }

  return(FALSE)
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

ComputeGlycanMass <- function(glycanComposition) {
  if(is.na(glycanComposition) || glycanComposition == "") return(NA)

  # 1. Extract Name(Count) pairs
  matches <- stringr::str_match_all(glycanComposition, "([A-Za-z0-9]+)\\(([0-9]+)\\)")[[1]]

  if (length(matches) == 0) return(0)

  # 2. Create data frame of found glycans
  found_counts <- data.frame(
    ShortName = matches[,2],
    count = as.numeric(matches[,3]),
    stringsAsFactors = FALSE
  )

  # 3. Join with your database
  final_df <- found_counts %>%
    dplyr::left_join(.modEnv$GlycanDatabase, by = "ShortName")

  # 4. Check for missing monosaccharides
  missing_names <- final_df$ShortName[is.na(final_df$GMass)]

  if (length(missing_names) > 0) {
    warning(paste("The following monosaccharides were not found in GlycanDatabase:",
                  paste(unique(missing_names), collapse = ", ")))
  }

  # 5. Calculate total mass (using na.rm = TRUE to ignore missing database entries)
  totalMass <- sum(final_df$count * final_df$GMass, na.rm = TRUE)

  return(totalMass)
}

ComputeGlycanMass_legacy <- function(glycanComposition){
  if(is.na(glycanComposition) || glycanComposition == ""){
    return(NA)
  }
  glycanMass_df <- .modEnv$GlycanDatabase %>%
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
    uniqueRundf$colName <- paste(uniqueRundf$Run, ".Intensity", sep = "")
  }else if(normalization == "FP_MaxLFQ"){
    uniqueRundf$colName <- paste(uniqueRundf$Run, ".MaxLFQ.Intensity", sep = "")
  }else{
    stop("The normalization was not recognized: ", normalization)
  }

  expected <- uniqueRundf$colName
  have_exact  <- expected %in% names(quantdata)
  have_noX    <- expected %in% sub("^X", "", names(quantdata))

  if (all(have_exact)) {

  } else if (all(have_noX)) {
    names(quantdata) <- sub("^X", "", names(quantdata))
  } else {
    notfound <- uniqueRundf$Run[!have_exact]

    warning(
      "Did not find the following runs in combined_peptide.tsv:\n",
      paste(" -", notfound, collapse = "\n")
    )
  }

  quantdata <- quantdata %>%
    dplyr::select("ModifiedPeptide" = "Modified.Sequence",
                  dplyr::any_of(uniqueRundf$colName)) %>%
    dplyr::mutate(ModifiedPeptide = gsub("\\[57\\.0214\\]|\\[57\\.0215\\]", "", .data$ModifiedPeptide),
                  dplyr::across(tidyselect::all_of(uniqueRundf$colName), as.numeric)) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(uniqueRundf$colName), names_to = "colName", values_to = "Intensity") %>%
    dplyr::filter(!is.na(.data$Intensity)) %>%
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

UpdateFPIntensitiesTMT <- function(clean_df, raw_df, quant_data){
  #Get a clean PTM table for matching
  clean_df_PTM <- PSMToPTMTable(clean_df, silent = TRUE)

  clean_df_PTM <- clean_df_PTM %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::distinct(.data$ModifiedPeptide, .data$ModificationID, .data$TotalGlycanComposition, .keep_all = TRUE) %>%
    dplyr::select(c("ModifiedPeptide", "ModificationID", "TotalGlycanComposition", "ID")) %>%
    dplyr::mutate(glycan_mass = purrr::map_dbl(.data$TotalGlycanComposition, ~ ComputeGlycanMass(.x)),
                  modification_site = as.numeric(gsub("\\D", "", .data$ModificationID))) %>%
    dplyr::mutate(glycan_mass = format(.data$glycan_mass, nsmall = 4)) %>%
    dplyr::left_join(raw_df[c("ID", "MSFragger.Localization", "Protein.ID")], by = "ID") %>%
    dplyr::arrange(.data$modification_site) %>%
    dplyr::summarise(.by = "ModifiedPeptide",
                     match_column = paste(unique(.data$Protein.ID),
                                          paste0(.data$ModificationID),
                                          paste0(.data$glycan_mass),
                                          sep = "-")) %>%
    dplyr::mutate(match_column = gsub(" ", "", .data$match_column))

  #First Clean the index column
  quant_data <- quant_data %>%
    dplyr::mutate(
      match_column = sapply(
        strsplit(.data$Index, "_"),
        function(x) paste(x[c(1,6,8)], collapse = "-")))

  #This shennanigans is to correct for potential subtle mass differences
  #missing_keys <- tibble(match_column = setdiff(quant_data$match_column, clean_df_PTM$match_column))
  missing_keys <- data.frame(match_column = setdiff(quant_data$match_column,
                                                    clean_df_PTM$match_column))

  if (nrow(missing_keys) > 0) {

    # Parse quant missing keys into clean structure
    quant_parsed <- missing_keys %>%
      tidyr::separate_wider_delim(
        cols = "match_column", delim = "-", names = c("Protein", "Site", "Mass"), cols_remove = FALSE
      ) %>%
      dplyr::mutate(Mass_num = as.numeric(.data$Mass))

    # Parse ALL unique clean_df_PTM keys into a reference structure
    ptm_parsed <- data.frame(match_column = unique(clean_df_PTM$match_column)) %>%
      tidyr::separate_wider_delim(
        cols = "match_column", delim = "-", names = c("Protein", "Site", "Mass"), cols_remove = FALSE
      ) %>%
      dplyr::mutate(Mass_num = as.numeric(.data$Mass))

    # Map them using a numeric proximity loop
    repair_lookup <- quant_parsed %>%
      dplyr::inner_join(ptm_parsed, by = c("Protein", "Site"), suffix = c("_old", "_new")) %>%
      dplyr::mutate(mass_dev = abs(.data$Mass_num_old - .data$Mass_num_new)) %>%
      dplyr::filter(.data$mass_dev <= 0.01) %>%
      dplyr::group_by(.data$match_column_old) %>%
      dplyr::slice_min(.data$mass_dev, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select("match_column_old", "match_column_new")

    quant_data <- quant_data %>%
      dplyr::left_join(repair_lookup, by = c("match_column" = "match_column_old")) %>%
      dplyr::mutate(
        match_column = dplyr::coalesce(.data$match_column_new, .data$match_column)
      ) %>%
      dplyr::select(-c("match_column_new"))
  }

  start_data <- which(names(quant_data) == "ReferenceIntensity") + 1
  quant_data <- quant_data[,start_data:ncol(quant_data)]

  # Start with the longest file name and match as greedy as possible
  # Then go with the next. This is to ensure the script works for names
  # containing -
  re_escape <- function(string) {
    gsub("([.\\|+*?\\[^\\]$(){}=!<>|:-])", "\\\\\\1", string)
  }

  filenames <- unique(clean_df$Run)
  search_ids <- names(quant_data)
  search_ids <- search_ids[search_ids != "match_column"]

  search_ids <- search_ids[order(nchar(search_ids), decreasing = TRUE)]

  mapping_results <- data.frame(
    Filename = character(),
    Matched_ID = character(),
    stringsAsFactors = FALSE
  )

  available_filenames <- filenames

  for (id in search_ids) {
    if (length(available_filenames) == 0) break

    matches <- stringr::str_ends(available_filenames, stringr::fixed(id))

    if (any(matches)) {
      matched_files <- available_filenames[matches]

      new_rows <- data.frame(Filename = matched_files, Matched_ID = id)
      mapping_results <- rbind(mapping_results, new_rows)

      available_filenames <- available_filenames[!matches]
    }
  }

  quant_runs <- length(search_ids)
  mapping_runs <- nrow(mapping_results)
  quant_ids <- length(unique(quant_data$match_column))
  psm_ids <- length(unique(clean_df_PTM$match_column))

  if(mapping_runs != quant_runs) {
    warning("Different number of psm.tsv runs and abundance_multi-mass_MD runs found.")
  }

  fmessage(paste0("Intensity values found for ", quant_ids, " of the ",
                  psm_ids, " peptides."))

  #Now replace the column names and do the actual mapping
  rename_vector <- stats::setNames(mapping_results$Filename, mapping_results$Matched_ID)

  cols_to_rename <- colnames(quant_data) %in% names(rename_vector)

  colnames(quant_data)[cols_to_rename] <- rename_vector[colnames(quant_data)[cols_to_rename]]

  quant_data <- quant_data %>%
    tidyr::pivot_longer(cols = colnames(quant_data)[colnames(quant_data) != "match_column"],
                        names_to = "Run", values_to = "Intensity") %>%
    dplyr::mutate(Intensity = ifelse(.data$Intensity > 0 & is.finite(.data$Intensity),
                                     2^.data$Intensity, .data$Intensity))

  number_valid_quant <- sum(is.finite(quant_data$Intensity))
  nrow_clean_df <- nrow(clean_df)

  clean_df <- clean_df %>%
    dplyr::left_join(clean_df_PTM[c("ModifiedPeptide", "match_column")] %>% dplyr::distinct(),
                     by = "ModifiedPeptide") %>%
    dplyr::select(-c("Intensity")) %>%
    dplyr::left_join(quant_data, by = c("Run", "match_column")) %>%
    dplyr::select(-c("match_column"))

  number_valid_quant_post <- sum(is.finite(quant_data$Intensity))
  nrow_clean_df_post <- nrow(clean_df)

  if(number_valid_quant != number_valid_quant_post) {
    warning(paste0("Number of valid values before and after import: ", number_valid_quant, " - ", number_valid_quant_post))
  }
  if(nrow_clean_df != nrow_clean_df_post) {
    warning(paste0("Number of rows before and after import: ", nrow_clean_df, " - ", nrow_clean_df_post))
  }

  return(clean_df)
}

GetPSMGlycanCategory <- function(GType, color_df = .modEnv$GlycanColors) {
  categories <- color_df$GlycanType[color_df$GlycanType != "Multi"]
  regex_pattern <- paste(categories, collapse = "|")

  found_matches <- stringr::str_extract_all(GType, regex_pattern)
  match_counts <- sapply(found_matches, function(x) length(x))

  as.vector(dplyr::case_when(
    match_counts > 1 ~ "Multi",
    match_counts == 0 ~ "nonGlyco",
    TRUE ~ sapply(found_matches, function(x) x[1])
  ))
}

FilterForMinPeptides <- function(df, minPeptideCoverage, thresholdMode = c("group", "total"), silent = FALSE){
  uniquePeptidesStart <- length(unique(df$ModifiedPeptide))
  thresholdMode <- match.arg(thresholdMode)
  df_work <- GetMeanTechReps(df)

  if(thresholdMode == "total"){
    if(length(minPeptideCoverage) == 1 && is.finite(minPeptideCoverage) && minPeptideCoverage == as.numeric(minPeptideCoverage)){

      maxValuesInGroup <- length(unique(df_work$Run))

      if(minPeptideCoverage < 1){
        minPeptideCoverage <- minPeptideCoverage * maxValuesInGroup}

      if(minPeptideCoverage > maxValuesInGroup){
        minPeptideCoverage <- maxValuesInGroup
        fmessage(paste0("The minimum value for the minPeptideCoverage is larger than the largest group. New value: ", maxValuesInGroup))
      }

      tokeep <- df_work %>%
        dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all=TRUE) %>%
        dplyr::summarise(.by = c("ModifiedPeptide"), count = dplyr::n()) %>%
        dplyr::filter(.data$count >= minPeptideCoverage) %>%
        dplyr::pull(.data$ModifiedPeptide)

      df <- df %>%
        dplyr::filter(.data$ModifiedPeptide %in% tokeep)

      if(!identical(silent, TRUE)){
        fmessage(paste0("Minimum found values filter from ", uniquePeptidesStart, " to ", length(unique(df$ModifiedPeptide)), " unique peptides."))
      }
      return(df)
    }else{stop("Please provide a single numeric value for minPeptideCoverage (absolute or fraction).")}
  }else{
    if(length(minPeptideCoverage) == 2 &&
       is.numeric(minPeptideCoverage) &&
       all(is.finite(minPeptideCoverage))){

      numberOfGroups <- length(unique(df_work$Condition))
      maxValuesInGroup <- df_work %>%
        dplyr::distinct(.data$Condition, .data$BioReplicate) %>%
        dplyr::summarise(.by = "Condition", total = dplyr::n()) %>%
        dplyr::pull(.data$total) %>%
        max()

      if(minPeptideCoverage[1] < 1){
        minPeptideCoverage[1] <- minPeptideCoverage[1] * maxValuesInGroup}

      if(minPeptideCoverage[1] > maxValuesInGroup){
        minPeptideCoverage[1] <- maxValuesInGroup
        fmessage(paste0("The minimum value for the minPeptideCoverage is larger than the largest group. New value: ", maxValuesInGroup))
      }

      if(minPeptideCoverage[2] > numberOfGroups){
        minPeptideCoverage[2] <- numberOfGroups
        fmessage(paste0("The minimum number of groups for the minPeptideCoverage is more than the nr of groups. New value: ", numberOfGroups))
      }

      tokeep <- df_work %>%
        dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all=TRUE) %>%
        dplyr::summarise(.by = c("Condition", "ModifiedPeptide"), count = dplyr::n()) %>%
        dplyr::filter(.data$count >= minPeptideCoverage[1]) %>%
        dplyr::summarise(.by = c("ModifiedPeptide"),
                         newCount = dplyr::n()) %>%
        dplyr::filter(.data$newCount >= minPeptideCoverage[2]) %>%
        dplyr::pull(.data$ModifiedPeptide)

      df <- df %>%
        dplyr::filter(.data$ModifiedPeptide %in% tokeep)

      if(!identical(silent, TRUE)){
        fmessage(paste0("Minimum found values filter from ", uniquePeptidesStart, " to ", length(unique(df$ModifiedPeptide)), " peptides."))
      }
      return(df)
    }else{
      stop("Please provide a vector containing two numeric values for filtering")
    }
  }
}

GetIntersections <- function(inputList, nintersects){
  sets <- names(inputList)
  allCombs <- list()

  for(k in 1:length(sets)){
    comb_k <- utils::combn(sets, k, simplify = FALSE)
    allCombs <- c(allCombs, comb_k)
  }

  result <- list()

  for(combo in allCombs){
    in_all <- Reduce(intersect, inputList[combo])

    outside_sets <- setdiff(sets, combo)
    if(length(outside_sets) > 0){
      not_in_others <- in_all[!in_all %in% unlist(inputList[outside_sets])]
    } else {
      not_in_others <- in_all
    }

    result[[paste(combo, collapse = "&")]] <- not_in_others
  }

  returndf <- utils::stack(result) %>%
    dplyr::rename("ModifiedPeptide" = "values",
                  "Intersect" = "ind")

  tokeep <- returndf %>%
    dplyr::summarise(.by = "Intersect",
                     count = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(.data$count)) %>%
    dplyr::slice(1:nintersects) %>%
    dplyr::pull(.data$Intersect) %>%
    as.character()

  returndf <- returndf %>%
    dplyr::mutate(Intersect = as.character(.data$Intersect)) %>%
    dplyr::filter(.data$Intersect %in% tokeep)

  return(returndf)
}

CompToGlyToucan <- function(compositionVec, returnType = "GlyToucan") {
  compositionDf <- data.frame(TotalGlycanComposition = compositionVec) %>%
    dplyr::distinct(.data$TotalGlycanComposition) %>%
    tidyr::separate_longer_delim(cols = "TotalGlycanComposition", delim = ",") %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$TotalGlycanComposition != "" & !is.na(.data$TotalGlycanComposition)) %>%
    dplyr::mutate(WURCS = NA, GlyToucan = NA)

  compositionDf$WURCS <- sapply(compositionDf$TotalGlycanComposition, function(x) getWURCS(x))
  compositionDf$GlyToucan <- sapply(compositionDf$WURCS, function(x) getGlyToucan(x))
  compositionDf$GlyToucan <- stringr::str_extract(compositionDf$GlyToucan, "(?<=\"id\":\\s\")[A-Z0-9]+")

  if(returnType == "GlyToucan") {
    return_df <- data.frame(TotalGlycanComposition = compositionVec) %>%
      dplyr::distinct(.data$TotalGlycanComposition) %>%
      dplyr::filter(.data$TotalGlycanComposition != "" & !is.na(.data$TotalGlycanComposition))

    lookup_vec <- stats::setNames(compositionDf$GlyToucan, compositionDf$TotalGlycanComposition)
    testcase <- strsplit(return_df$TotalGlycanComposition, ",")
    replaced_list <- lapply(testcase, function(x) {
      matches <- lookup_vec[x]
      ifelse(is.na(matches), x, matches)
    })

    return_df$GlyToucan <- sapply(replaced_list, paste, collapse = ",")

    return(return_df)
  }
}

getWURCS <- function(comp) {
  baseUrl <- "https://api.glycosmos.org/glycancompositionconverter/1.0.0/composition2wurcs/"

  match_df <- data.frame(name = c("hex", "hexnac", "dhex", "neu5ac", "neu5gc", "P",
                                  "S", "Ac"),
                         dbName =  c("H", "N", "F", "A", "G", "Phospho",
                                     "Sulfo", "Acetyl"))

  matches <- stringr::str_match_all(comp, "([A-Za-z0-9]+)\\(([0-9]+)\\)")[[1]]
  if (length(matches) == 0) return(NA)

  # 2. Create data frame of found glycans
  found_counts <- data.frame(
    dbName = matches[,2],
    count = as.numeric(matches[,3]),
    stringsAsFactors = FALSE
  )

  # Using left_join with the mapping table you provided
  found_counts <- found_counts %>%
    dplyr::left_join(match_df, by = "dbName") %>%
    dplyr::mutate(ShortName = ifelse(is.na(.data$dbName), .data$name, .data$dbName))

  # 5. Check for missing database entries
  missing_names <- found_counts$dbName[is.na(found_counts$name)]
  if (length(missing_names) > 0) {
    return(NA)
  }

  #Now construct the url
  endUrl <- found_counts %>%
    dplyr::mutate(se = paste(.data$name, .data$count, sep = ":")) %>%
    dplyr::summarise(endUrl = paste(.data$se, collapse = "%7C")) %>%
    dplyr::pull("endUrl")

  finalUrl <- paste0(baseUrl, endUrl)

  response <- tryCatch({
    httr::GET(finalUrl, httr::timeout(10)) # Added a 10-second timeout
  }, error = function(e) {
    message("Connection error: ", e$message)
    return(NULL)
  })

  if (is.null(response)) return(NA)

  if (httr::status_code(response) != 200) {
    warning(paste("API Error:", httr::status_code(response), "at URL:", finalUrl))
    return(NA)
  }

  wurcs_result <- httr::content(response, as = "text", encoding = "UTF-8")

  wurcs_result
}

getGlyToucan <- function(wurcs) {
  baseUrl <- "https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs="
  fullUrl <- paste0(baseUrl, wurcs)

  response <- tryCatch({
    httr::GET(fullUrl, httr::timeout(10)) # Added a 10-second timeout
  }, error = function(e) {
    message("Connection error: ", e$message)
    return(NULL)
  })

  if (is.null(response)) return(NA)

  if (httr::status_code(response) != 200) {
    warning(paste("API Error:", httr::status_code(response), "at URL:", fullUrl))
    return(NA)
  }

  GlyToucan_result <- httr::content(response, as = "text", encoding = "UTF-8")

  GlyToucan_result
}

getModifiedPeptideForGlycanFinder <- function(rawModPep){
  formatted_column <-length(rawModPep)

  for(i in 1:length(rawModPep)){
    mod_in_string <- rawModPep[i]
    clean_mod <- c()
    while(grepl("\\[", mod_in_string)) {
      this_mod <- stringr::str_extract(mod_in_string, "\\[.*?\\]")
      site_number <- as.integer(regexpr(this_mod, mod_in_string, perl = TRUE)) - 2
      aa <- substring(mod_in_string, site_number, site_number)
      mass <- gsub("\\[|\\]", "", this_mod)
      mod_in_string <- sub(this_mod, "", mod_in_string, fixed = TRUE)

      if(length(clean_mod) == 0) {
        clean_mod <- paste0(site_number, aa, "(", mass, ")")
      } else {
        clean_mod <- paste0(clean_mod, ",", paste0(site_number, aa, "(", mass, ")"))
      }
    }
    formatted_column[i] <- clean_mod
  }
  formatted_column
}

GetProteinAndGeneGlycanFinder <- function(raw_col){
  raw_col <- strsplit(raw_col, ";")
  raw_col <- lapply(raw_col, function(x) sub("_.*", "", x))
  genes <- unlist(lapply(raw_col, function(x) paste(sub(".*\\|", "", x), collapse = ",")))
  prots <- unlist(lapply(raw_col, function(x) paste(sub("\\|.*", "", x), collapse = ",")))

  combined <- paste(prots, genes, sep = ";")
  combined
}

GetProteinStartGlycanFinder <- function(glyco_in_peptide, glyco_in_protein) {
  glyco_in_protein <- sub(";.*", "",  glyco_in_protein)
  glyco_in_protein <- sub("_.*", "",  glyco_in_protein)
  glyco_in_protein <- as.integer(gsub("\\D", "", glyco_in_protein))
  glyco_in_peptide <- as.integer(sub("_.*", "",  glyco_in_peptide))

  glyco_in_protein - glyco_in_peptide + 1
}

GetAScore <- function(Ascore) {
  Ascore <- strsplit(Ascore, ";")
  Ascore <- lapply(Ascore, function(x) as.double(sub(".*:", "", x)))
  highest <- sapply(Ascore, max)
  highest
}

Databases <- function(){
  GlycanDatabase <- data.frame(
    FullName = c("HexNAc", "Hex", "NeuAc", "Fuc", "NeuGc", "Pent",
                 "KDN", "HexA", "pseudaminic",
                 "Phospho", "Sulfo", "Acetyl", "NH4", "Na", "Fe", "Ca",
                 "Al", "K", "Formyl", "Succinyl"),
    ShortName = c("N", "H", "A", "F", "G", "Pen", "Kdn", "HexA", "p",
                  "P", "S", "C", "NH4", "Na", "Fe", "Ca",
                  "Al", "K", "M", "U"),
    GMass = c(203.07937, 162.05282, 291.09542, 146.05791, 307.09033, 132.0423,
              250.06889, 176.03209, 232.10592,
              79.96633, 79.95682, 42.01056, 17.02655, 21.98194, 52.91146,
              37.94694, 23.95806, 38.96371, 27.99491, 100.01608),
    stringsAsFactors = FALSE)

  ModificationDatabase <- data.frame(
    FullName = c("Oxidation", "Oxidation2", "CCarbamidomethylation1", "CCarbamidomethylation2",
                 "CCarbamidomethylation3",
                 "NAcetylation", "Carbamidomethyl", "TMT0",
                 "TMT2", "TMT6", "TMT10", "TMT11", "TMT16", "TMT18", "TMT35"),
    ModificationMass = c("15.9949", "15.99", "57.0214", "57.0215", "57.02", "42.0106", "57.0215", "295.18959",
                         "225.1558", "229.1629", "229.1629", "229.1629", "304.2071", "304.2071", "304.2071")
  )

  GlycanColors = data.frame(GlycanType = c("Complex/Hybrid", "Sialofucosylated", "Sialylated",
                                           "Fucosylated", "Oligomannose", "Truncated",
                                           "Paucimannose", "Phosphomannose", "OGlycan",
                                           "NonCanonicalGlyco", "Multi"),
                            color = c("#A1CAE8", "#FFE8A2", "#CC82C3",
                                      "#FFA1A1", "#A0E2BE", "#686963",
                                      "#0072BC", "#EE8866", "#BF5A6B",
                                      "#FABC3C", "#664C43"))

#   original <- c(
#     "#baa5cc", "#9adcee", "#44aa99", "#ddcc77",
#     "#ffaabb", "#cc6677", "#882255", "#32006e"
#   )
#    n_new <- 100
#    interp_fun <- colorRampPalette(original)
#    extra <- interp_fun(length(original) + n_new)[-(1:length(original))]
#    set.seed(3)
#    #scales::show_col(sample(extra))
#    extra <- sample(extra)
#    colorScheme <- c(original, extra)
#
#   usethis::use_data(GlycanDatabase, colorScheme, ModificationDatabase, GlycanColors, internal = TRUE, overwrite = TRUE)
}
