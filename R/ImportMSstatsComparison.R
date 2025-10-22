ImportMSstatsComparison <- function(path, cleanCCarbamidomethylation = TRUE){
  comparison_df <- data.frame(UniprotIDs = rep(NA_character_, nrow(MSstats_raw)))

  MSstats_raw <- read.csv(path)
  existingCols <- names(MSstats_raw)

  #UniprotIDs####
  if("Protein" %in% existingCols){
    comparison_df$UniprotIDs <- MSstats_raw %>%
      dplyr::rowwise() %>%
      dplyr::mutate(returnVal = stringr::str_extract(.data$Protein,  "(?<=\\|)[^|]+(?=\\|)")) %>%
      dplyr::pull(.data$returnVal)
    fmessage("Successfully imported Assigned Modifications column.")
  }else{stop("The column Protein was not found in the input dataframe.")}

  #Proteins####
  if("proteinName" %in% existingCols){
    comparison_df$Proteins <- as.character(MSstats_raw$proteinName)
    fmessage("Successfully imported Protein column.")
  }else{warning("The column Protein was not found in the input dataframe.")}

  #ModificationID####
  if("Protein" %in% existingCols){
    comparison_df$ModificationID <- MSstats_raw %>%
      dplyr::mutate(returnVec = sub(".*_(.*)$", "\\1", .data$Protein)) %>%
      dplyr::pull(.data$returnVec)
    fmessage("Successfully imported ModificationID column.")
  }else{warning("The column Protein was not found in the input dataframe.")}

  if("PeptideSequence" %in% existingCols){
    comparison_df$ModifiedPeptide <- as.character(MSstats_raw$PeptideSequence)
    fmessage("Successfully imported ModifiedPeptide column.")
  }else{stop("The column PeptideSequence was not found in the input dataframe.")}

  if(cleanCCarbamidomethylation){
    comparison_df$ModifiedPeptide <- gsub("[57.0215]", "", comparison_df$ModifiedPeptide,
                                          fixed = TRUE)
    fmessage("Successfully removed the string [57.0215] from the modified peptide sequence.")
  }

  if("Label" %in% existingCols){
    comparison_df$Label <- as.character(MSstats_raw$Label)
    fmessage("Successfully imported Label column.")
  }else{warning("The column Label was not found in the input dataframe.")}

  if("log2FC" %in% existingCols){
    comparison_df$log2FC <- as.double(MSstats_raw$log2FC)
    fmessage("Successfully imported log2FC column.")
  }else{warning("The column log2FC was not found in the input dataframe.")}

  if("pvalue" %in% existingCols){
    comparison_df$pvalue <- as.double(MSstats_raw$pvalue)
    fmessage("Successfully imported pvalue column.")
  }else{warning("The column pvalue was not found in the input dataframe.")}

  if("adj.pvalue" %in% existingCols){
    comparison_df$adjpvalue <- as.double(MSstats_raw$adj.pvalue)
    fmessage("Successfully imported adjpvalue column.")
  }else{warning("The column adj.pvalue was not found in the input dataframe.")}

  return(comparison_df)
}
