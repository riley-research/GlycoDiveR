#' ComputeMissingness
#'
#' This computes the percentage of missing values in your GlycoDiveR uploaded data.
#' The computation is done on a (glyco)peptide level. Peptides with a quantitative
#' value of 0 are considered missing.
#'
#' @param input Your formatted GlycoDiveR data.
#' @param type Choose between "glyco" or "all"
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param silent silence printed information (default = FALSE)
#'
#' @returns The percentage of missing data.
#' @export
#'
#' @examples \dontrun{ComputeMissingness(mydata, type = "glyco")}
ComputeMissingness <- function(input, type = "glyco", whichPeptide = NA,
                                    whichAlias = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  df <- FilterForPeptides(input$PSMTable, whichPeptide)
  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(type == "glyco"){
    glycoPSMTypes <- c("Sialyl", "Complex/Hybrid", "Sialyl+Fucose",
                       "Fucose", "Truncated", "High Mannose", "Paucimannose",
                       "OGlycan")
    df <- input$PSMTable %>%
      dplyr::mutate(PSMType = dplyr::case_when(stringr::str_count(.data$GlycanType, paste(glycoPSMTypes, collapse = "|")) > 1 &
                                                 stringr::str_count(.data$GlycanType, "Sialyl") == 2 &
                                                 grepl("Sialyl+Fucose", .data$GlycanType) == 1 ~ "Multi",
                                               stringr::str_count(.data$GlycanType, paste(glycoPSMTypes, collapse = "|")) == 0 ~ "nonGlyco",
                                               grepl("Sialyl", .data$GlycanType) ~ "Sialyl",
                                               grepl("Complex/Hybrid", .data$GlycanType) ~ "Complex/Hybrid",
                                               grepl("Sialyl+Fucose", .data$GlycanType) ~ "Sialyl+Fucose",
                                               grepl("Fucose", .data$GlycanType) ~ "Fucose",
                                               grepl("Truncated", .data$GlycanType) ~ "Truncated",
                                               grepl("High Mannose", .data$GlycanType) ~ "High Mannose",
                                               grepl("Paucimannose", .data$GlycanType) ~ "Paucimannose",
                                               grepl("OGlycan", .data$GlycanType) ~ "OGlycan",
                                               TRUE ~ "ERROR")) %>%
      dplyr::filter(.data$PSMType != "nonGlyco")
  }

  df <- df %>%
    dplyr::filter(!is.na(.data$Intensity)) %>%
    dplyr::distinct(.data$ModifiedPeptide, .data$Alias, .keep_all=TRUE) %>%
    dplyr::select("ModifiedPeptide", "Alias", "Intensity") %>%
    tidyr::complete(.data$Alias, .data$ModifiedPeptide)

  missingValues <- sum(is.na(df$Intensity))
  foundValues <- sum(!is.na(df$Intensity))

  return(paste0("Percentage of missing values: ", missingValues/(foundValues + missingValues) *100,"%"))
}
