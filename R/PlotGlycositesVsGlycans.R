#' PlotGlycositesVsGlycans
#'
#' Plots the proteins by the number of glycosites on the proteins to the
#' number of identified glycans identified on that protein.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param pointSize The size of the points.
#' @param labelGlycositeCutoff The minimal number of glycosites found to label
#' a protein (default = 0)
#' @param labelGlycanCutoff The minimal number of glycans found to label a protein
#' (default = 0)
#' @param maxOverlaps The maximum number allowed label overlaps (default = 10)
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
#' @param plotColor The color of plotted points.
#' @param alpha The alpha of the plotted points.
#' @param silent silence printed information (default = TRUE)
#'
#' @returns A scatter plot comparing the number of glycosites to the number of glycans
#' for each protein.
#' @export
#'
#' @examples \dontrun{
#' PlotGlycositesVsGlycans(mydata)
#'
#' PlotGlycositesVsGlycans(mydata, plotColor = "black", alpha = 0.1, maxOverlaps = Inf)}
PlotGlycositesVsGlycans <- function(input, whichAlias = NULL, pointSize = 3,
                                    labelGlycositeCutoff = 0, labelGlycanCutoff = 0,
                                    maxOverlaps = 10, whichPeptide = NULL,
                                    whichProtein = NULL, exactProteinMatch = TRUE,
                                    plotColor = c("#BAA5CC", "#32006e"), alpha = 0.5,
                                    silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

  if(!is.null(whichAlias)){
    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  #Prepare data
  df <- input$PTMTable %>%
    dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "" &
                    .data$GlycanType != "NonGlyco") %>%
    dplyr::mutate(.by = "UniprotIDs",
                  count = dplyr::n_distinct(.data$TotalGlycanComposition)) %>%
    dplyr::summarise(.by = c("UniprotIDs", "count", "Genes"),
                             numberOfGlycosites = dplyr::n_distinct(.data$ModificationID))

  if(nrow(df) == 0){
    return(fmessage("No data is left after filtering."))
  }

  #Get labels
  df <- df %>%
    dplyr::mutate(label = dplyr::case_when(.data$numberOfGlycosites > labelGlycositeCutoff &
                                      .data$count > labelGlycanCutoff ~ .data$Genes,
                                    TRUE ~ NA))
  # Plot
  df <- df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$numberOfGlycosites, y =.data$count)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(fill = plotColor[1], color = plotColor[2], shape = 21, alpha = alpha,
                        size = pointSize) +
    suppressWarnings(ggrepel::geom_label_repel(ggplot2::aes(label = .data$label), max.overlaps = maxOverlaps,
                              fill = NA, label.size = NA)) +
    ggplot2::annotate("text", x= 0.95 * max(df$numberOfGlycosites), y= 1.1 * max(df$numberOfGlycosites),
                      label= "y=x", color = "grey30") +
    ggplot2::labs(x = "Number of glycosites", y = "Number of glycans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

  return(df)
}
