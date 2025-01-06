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
      dplyr::filter(!is.na(AssignedModifications) & AssignedModifications != "") %>%
      tidyr::separate_rows(AssignedModifications, sep = ",")

  tempdf$AssignedModifications <- gsub("N-term", "1", tempdf$AssignedModifications)

  tempdf <- tempdf %>%
      dplyr::rowwise() %>%
      dplyr::mutate(PeptidePTMLocalization = as.numeric(regmatches(AssignedModifications, gregexpr("^[0-9]+", AssignedModifications))[[1]]),
                    ProteinPTMLocalization = PeptidePTMLocalization + ProteinStart,
                    ModificationSite = sub(".*([A-Za-z])\\(.*", "\\1", AssignedModifications),
                    ModificationID = paste0(ModificationSite,ProteinPTMLocalization))

  tempdf$GlycanType <- apply(tempdf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))

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
        hexNAc_count <- suppressWarnings(as.numeric(sub(".*HexNAc\\(([0-9]+)\\).*", "\\1", glycanComp)))
        hex_count <- suppressWarnings(as.numeric(sub(".*Hex\\(([0-9]+)\\).*", "\\1", glycanComp)))

        glycanCat <- dplyr::case_when(
          grepl("A", glycanComp) & grepl("F", glycanComp) ~ "Sialyl+Fucose",
          grepl("A", glycanComp) ~ "Sialyl",
          grepl("F", glycanComp) ~ "Fucose",
          !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 3 & hex_count > 4 ~ "High Mannose",
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
    dplyr::arrange(desc(Intensity)) %>%
    dplyr::distinct(Run, AssignedModifications,ModifiedPeptide, Condition, BioReplicate, TechReplicate, .keep_all = TRUE) %>%
    ungroup()

  #Take median of technical reps
  df <- df %>%
    dplyr::group_by(ModifiedPeptide, AssignedModifications, Condition, BioReplicate) %>%
    dplyr::mutate(Intensity = stats::median(Intensity, na.rm = TRUE)) %>%
    distinct(ModifiedPeptide, Condition, BioReplicate, .keep_all = TRUE)

  return(df)
}

CleanGlycanNames <- function(glycan){
  glycan <- gsub("HexNAc", "N", glycan)
  glycan <- gsub("Hex", "H", glycan)
  glycan <- gsub("NeuAc", "A", glycan)
  glycan <- gsub("Fuc", "F", glycan)
  glycan <- gsub("NeuGc", "G", glycan)
  return(glycan)
}

FilterForCutoffs <- function(input){
  if(input$searchEngine %in% c("MSFragger")){
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Filtering for PSMScore >= ", test$peptideScoreCutoff, " and glycan score <= ", test$glycanScoreCutoff, ".\033[0m")
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Prefilter number of rows PSM table: ", nrow(input$PSMTable), ". Prefilter number of rows PTM table: ", nrow(input$PTMTable), ".\033[0m")
    input$PSMTable <- subset(input$PSMTable, (PSMScore >= test$peptideScoreCutoff & GlycanQValue <= test$glycanScoreCutoff) |
                               (PSMScore >= test$peptideScoreCutoff & is.na(GlycanQValue)))
    input$PTMTable <- subset(input$PTMTable, (PSMScore >= test$peptideScoreCutoff & GlycanQValue <= test$glycanScoreCutoff) |
                               (PSMScore >= test$peptideScoreCutoff & is.na(GlycanQValue)))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Postfilter number of rows PSM table: ", nrow(input$PSMTable), ". Postfilter number of rows PTM table: ", nrow(input$PTMTable), ".\033[0m")
    return(input)
  }else{
    warning("No search engine recognized, returning unfiltered dataframe.")
    return(input)
  }
}
