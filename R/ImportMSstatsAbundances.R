#' ImportMSstatsAbundances
#'
#' This script uses MSstats generated quantitative data and uses this data to fill
#' in the existing quantitative values. It will only retain quantitative values
#' present in the MSstats data. Requirements are 1) matching Alias column with the
#' same values as in your GlycoDiveR data, 2) MSstats data has a LogIntensities
#' column, 3) MSstats data has a ModifiedPeptide column with values correspoding
#' to the ModifiedPeptide column in your data. This script is for MSstats data,
#' but should work for other inputs too.
#'
#' @param input This is your formatted GlycoDiveR data
#' @param MSstatsPath The filepath to the MSstats generated csv.
#' @param cleanCCarbamidomethylation Removes \\[57.0215\\] from the MSstats
#' ModifiedPeptide column (default = TRUE)
#'
#' @returns your input data with updated quantitative values.
#' @export
#'
#' @examples \dontrun{ImportMSstatsComparison(mydata, "C:/User/Pathtocsv.csv")}
ImportMSstatsAbundances <- function(input, MSstatsPath, cleanCCarbamidomethylation = TRUE){
  #1. Import MSstats files
    #Requirements are 1) matching Alias column, 2) LogIntensities column, 3) ModifiedPeptide column
  #2. Check validity with input
  #3. Remove existing intensity values in PSMTable
  #4. Use the MSStats intensity values
  #5. Recalculate PTMTable
  MSstats_raw <- utils::read.csv(MSstatsPath) %>%
    dplyr::select("Alias" = "originalRUN", "ModifiedPeptide", "LogIntensities") %>%
    dplyr::mutate(Intensity = 2 ^ .data$LogIntensities,
                  Alias = factor(Alias, levels = levels(input$PSMTable$Alias))) %>%
    dplyr::select(-"LogIntensities")

  if(cleanCCarbamidomethylation){
    MSstats_raw$ModifiedPeptide <- gsub("[57.0215]", "", MSstats_raw$ModifiedPeptide,
                                          fixed = TRUE)
    fmessage("Successfully removed the string [57.0215] from the modified peptide sequence.")
  }

  if(!all(MSstats_raw$ModifiedPeptide %in% input$PSMTable$ModifiedPeptide)){
    warning("Not all ModifiedPeptide IDs in the MSstats file are found in your GlycoDiveR data.")
  }

  input$PSMTable <- input$PSMTable %>%
    dplyr::select(-"Intensity")

  input$PSMTable <- input$PSMTable %>%
    dplyr::full_join(MSstats_raw, by = c("Alias", "ModifiedPeptide"))

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
    dplyr::mutate(GlycanQValue = 0)

  nrowPSMTable <- input$PSMTable %>%
    dplyr::filter(!is.na(.data$Intensity)) %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
    nrow()

  nrowMSstats <- MSstats_raw %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
    nrow()

  if(nrowPSMTable != nrowMSstats){
    warning("Some peptides were not detected after combining.")
  }

  PTMdf <- PSMToPTMTable(input$PSMTable)
  PTMdf$TotalGlycanComposition <- sapply(PTMdf$TotalGlycanComposition, function(x) strsplit(x, " % ")[[1]][1])

  input$PTMTable <- PTMdf

  return(input)
}
