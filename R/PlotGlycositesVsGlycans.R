#' PlotGlycositesVsGlycans
#'
#' Plots the proteins by the number of glycosites on the proteins to the
#' number of identified glycans identified on that protein.
#'
#' @param input The formatted GlycoDiveR data.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param labelGlycositeCutoff The minimal number of glycosites found to label
#' a protein (default = 0)
#' @param labelGlycanCutoff The minimal number of glycans found to label a protein
#' (default = 0)
#' @param maxOverlaps The maximum number allowed label overlaps (default = 10)
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#'
#' @returns A scatter plot comparing the number of glycosites to the number of glycans
#' for each protein.
#' @export
#'
#' @examples \dontrun{PlotGlycositesVsGlycans(mydata)}
PlotGlycositesVsGlycans <- function(input, whichAlias = NULL,
                                    labelGlycositeCutoff = 0, labelGlycanCutoff = 0,
                                    maxOverlaps = 10, whichPeptide = NA){
  input <- FilterForCutoffs(input)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

  if(!is.null(whichAlias)){
    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  #Prepare data
  df <- input$PTMTable %>%
    dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "" &
                    .data$GlycanType != "NonGlyco") %>%
    dplyr::mutate(.by = .data$UniprotIDs,
                  count = dplyr::n_distinct(.data$TotalGlycanComposition)) %>%
    dplyr::summarise(.by = c(.data$UniprotIDs, .data$count, .data$Genes),
                             numberOfGlycosites = dplyr::n_distinct(.data$ModificationID))

  #Get labels
  df <- df %>%
    dplyr::mutate(label = dplyr::case_when(.data$numberOfGlycosites > labelGlycositeCutoff &
                                      .data$count > labelGlycanCutoff ~ .data$Genes,
                                    TRUE ~ NA))
  # Plot
  df <- df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$numberOfGlycosites, y =.data$count)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(color = "#7b5799", alpha = 0.5) +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data$label), max.overlaps = 10,
                              fill = NA, label.size = NA) +
    ggplot2::annotate("text", x= 0.95 * max(df$numberOfGlycosites), y= 1.1 * max(df$numberOfGlycosites),
                      label= "y=x", color = "grey30") +
    ggplot2::labs(x = "Number of glycosites", y = "Number of glycans")

  return(df)
}
