#' ImportPerseusAbundances
#'
#' @param input
#' @param PerseusPath
#' @param cleanCCarbamidomethylation
#'
#' @returns
#' @export
#'
#' @examples
ImportPerseusAbundances <- function(input, PerseusPath, cleanCCarbamidomethylation = TRUE){
  columnNames <- names(utils::read.delim(PerseusPath))

  # Read the data using read.delim from text
  PerseusData <- utils::read.delim(PerseusPath, skip = 2)
  names(PerseusData) <- columnNames

  AliasInData <- unique(as.character(input$PSMTable$Alias))

  if(all(AliasInData %in% names(PerseusData))){
    fmessage("Found all Alias values in Perseus data.")
  }else{
    missing_aliases <- AliasInData[!AliasInData %in% names(PerseusData)]
    warning("The following Alias values were not found in Perseus data: ",
            paste(missing_aliases, collapse = ", "))
  }

  if(cleanCCarbamidomethylation){
    PerseusData$Modified.Sequence <- gsub("\\[57\\.0215\\]", "", PerseusData$Modified.Sequence)
  }

  if(!all(unique(PerseusData$Modified.Sequence) %in% unique(input$PSMTable$ModifiedPeptide))){
    missing_peptides <- unique(PerseusData$Modified.Sequence)[
      !unique(PerseusData$Modified.Sequence) %in% unique(input$PSMTable$ModifiedPeptide)]
    warning("The following Alias values were not found in Perseus data: ",
            paste(missing_peptides, collapse = ", "))
  }else{
    fmessage("All Modified Sequences in Perseus data are found in the GlycoDiveR dataset.")
  }

  PerseusData <- PerseusData %>%
    dplyr::select("ModifiedPeptide" = "Modified.Sequence",
                  dplyr::all_of(AliasInData)) %>%
    tidyr::pivot_longer(cols = AliasInData, names_to = "Alias", values_to = "Intensity") %>%
    dplyr::mutate(Intensity = 2^.data$Intensity,
                  Alias = factor(.data$Alias, levels = levels(input$PSMTable$Alias)))

  input$PSMTable <- input$PSMTable %>%
    dplyr::select(-"Intensity")

  input$PSMTable <- input$PSMTable %>%
    dplyr::full_join(PerseusData, by = c("Alias", "ModifiedPeptide"))

  #Fill the missing data with the rest of the data
  input$PSMTable <- input$PSMTable %>%
    dplyr::group_by(.data$ModifiedPeptide) %>%
    tidyr::fill(c("ModifiedPeptide", "AssignedModifications", "TotalGlycanComposition",
                  "IsUnique", "UniprotIDs", "Genes", "ProteinLength", "NumberOfNSites",
                  "NumberOfOSites", "ProteinStart", "GlycanType", "SubcellularLocalization",
                  "Domains", "RetentionTime", "ID"),
                .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Alias) %>%
    tidyr::fill(c("Condition", "BioReplicate", "TechReplicate", "Run")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(GlycanQValue = 0,
                  PSMScore = 1000)

  nrowPSMTable <- input$PSMTable %>%
    dplyr::filter(!is.na(.data$Intensity)) %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
    nrow()

  nrowPerseusData <- PerseusData %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
    nrow()

  if(nrowPSMTable != nrowPerseusData){
    warning("Some peptides were not detected after combining.")
  }

  PTMdf <- PSMToPTMTable(input$PSMTable)
  PTMdf$TotalGlycanComposition <- sapply(PTMdf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])

  input$PTMTable <- PTMdf

  return(input)
}
