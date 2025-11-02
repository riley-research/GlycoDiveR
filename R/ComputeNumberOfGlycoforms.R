#' ComputeNumberOfGlycoforms
#'
#' @param input Formatted data
#' @param UniprotID The UniprotID as in the UniprotIDs column
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#'
#' @returns The number of potential proteoforms
#' @export
#'
#' @examples \dontrun{ComputeNumberOfGlycoforms(mydata, "Q9NZQ7")}
ComputeNumberOfGlycoforms <- function(input, UniprotID, whichAlias = NULL){
  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df <- df %>%
    dplyr::filter(.data$UniprotIDs == UniprotID) %>%
    dplyr::distinct(.data$ModificationID, .data$TotalGlycanComposition) %>%
    dplyr::summarise(.by = .data$ModificationID, count = dplyr::n())

  if(nrow(df) == 0 | sum(df$count, na.rm = TRUE) == 0){
    return("None found")
  }else{
    df$count <- df$count + 1
    rslt = prod(df$count, na.rm = TRUE)
    return(rslt)
  }
}
