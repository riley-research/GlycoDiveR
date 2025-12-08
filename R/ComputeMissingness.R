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
#' @param removeAllNA Remove peptides with only NA values. This can, for example
#' happen if you upload abundance values using Perseus or MSstats.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param silent silence printed information (default = FALSE)
#'
#' @returns The percentage of missing data.
#' @export
#'
#' @examples \dontrun{ComputeMissingness(mydata, type = "glyco")}
ComputeMissingness <- function(input, type = "glyco", removeAllNA = TRUE,
                               whichPeptide = NULL, whichAlias = NULL,
                               whichProtein = NULL, exactProteinMatch = TRUE,
                               silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  df <- FilterForPeptides(input$PSMTable, whichPeptide)
  df <- FilterForProteins(df, whichProtein, exactProteinMatch)

  if(removeAllNA){
    df <- df %>%
      dplyr::filter(!is.na(.data$Intensity))
  }
  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(type == "glyco"){
    glycoPSMTypes <- c("Sialyl", "Complex/Hybrid", "Sialyl+Fucose",
                       "Fucose", "Truncated", "Oligomannose", "Paucimannose",
                       "OGlycan", "Phosphomannose", "NonCanonicalGlyco")
    df <- df %>%
      dplyr::mutate(PSMType = GetPSMGlycanCategory(.data$GlycanType)) %>%
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
