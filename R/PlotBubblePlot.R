#' PlotBubblePlot
#'
#' This is an explanation.
#'
#' @param input input
#' @param grouping  the goruping
#' @param whichAlias  the alias
#' @param whichPeptide  the peptide
#' @param whichProtein  the protein
#' @param exactProteinMatch  my match
#' @param silent is silent
#'
#' @returns a graph
#' @export
#'
#' @examples \dontrun{PlotBubblePlot(mydata)}
PlotBubblePlot <- function(input, grouping = "condition", whichAlias = NULL,
                              whichPeptide = NULL, whichProtein = NULL,
                              exactProteinMatch = TRUE, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  input$PSMTable <- input$PSMTable %>% dplyr::filter(!is.na(.data$Intensity))

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(input$PSMTable) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  my_Variable = 2

  View(input$PSMTable)
  View(input$PTMTable)

  p <- input$PTMTable |>
    dplyr::distinct(.data$Genes, .keep_all = TRUE) |>
    dplyr::mutate(PSMScore = .data$PSMScore * 100) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$PSMScore)) +
    ggplot2::geom_histogram() +
    ggplot2::coord_flip()

  return(p)


}
