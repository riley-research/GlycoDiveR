#' ComputeNumberOfGlycoforms
#'
#' @param input Formatted data
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
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
#' @returns The number of potential glycoforms.
#' @export
#'
#' @examples \dontrun{
#' ComputeNumberOfGlycoforms(mydata, whichProtein = "Q9NZQ7")}
ComputeNumberOfGlycoforms <- function(input, whichProtein = NULL,
                                      exactProteinMatch = TRUE, whichAlias = NULL,
                                      whichPeptide = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

  df <- input$PTMTable %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df <- df %>%
    dplyr::filter(.data$UniprotIDs == whichProtein) %>%
    dplyr::distinct(.data$ModificationID, .data$TotalGlycanComposition) %>%
    dplyr::summarise(.by = "ModificationID", count = dplyr::n())

  if(nrow(df) == 0 | sum(df$count, na.rm = TRUE) == 0){
    return("No glycoforms after filtering")
  }else{
    df$count <- df$count + 1
    rslt = prod(df$count, na.rm = TRUE)
    return(rslt)
  }
}
