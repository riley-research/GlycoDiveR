CombineTwoCols <- function(row){
  if (is.na(row[2]) | row[2] == ""){return(as.character(row[1]))}
  else{return(as.character(paste(row[1], row[2], sep = ",")))}}

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

  tempdf$NGlycanType <- apply(tempdf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))

  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Generated PTM table.\033[0m")

  return(tempdf)
}

GlycanComptToGlycanType <- function(mod, glycanComp){
  if(is.na(mod) | mod == ""){return("NotGlyco")}
  else{
    #glycanComp = "HexNAc(4)Hex(5)NeuAc(1) % 1913.6771"
    glycanMass = strsplit(glycanComp, "%")[[1]][2]
    glycanMass = substring(glycanMass, 2, nchar(glycanMass) - 1 )

    if(TRUE %in% grepl(glycanMass, mod)){
      hexNAc_count <- as.numeric(sub(".*HexNAc\\(([0-9]+)\\).*", "\\1", glycanComp))
      hex_count <- as.numeric(sub(".*Hex\\(([0-9]+)\\).*", "\\1", glycanComp))

      glycanCat <- case_when(
        grepl("NeuAc", glycanComp) & grepl("Fuc", glycanComp) ~ "Sialyl+Fucose",
        grepl("NeuAc", glycanComp) ~ "Sialyl",
        grepl("Fuc", glycanComp) ~ "Fucose",
        !is.na(hexNAc_count) & !is.na(hex_count) & hexNAc_count < 3 & hex_count > 4 ~ "High Mannose",
        TRUE ~ "Complex/Hybrid"
      )
      return(glycanCat)
    }else{
        return("NotGlyco")}
  }
}
