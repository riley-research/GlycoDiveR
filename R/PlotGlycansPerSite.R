#' PlotGlycanPerSite
#'
#' @param input Formatted data
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#'
#' @returns The number of glycans per site
#' @export
#'
#' @examples \dontrun{PlotGlycanPerSite(myData)}
PlotGlycansPerSite <- function(input, whichAlias = NULL, whichPeptide = NA){
  input <- FilterForCutoffs(input)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

  df <- GetMeanTechReps(input$PTMTable)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df <- df %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::summarise(.by = c(.data$UniprotIDs, .data$ModificationID),
                     GlycansPerSite = dplyr::n_distinct(.data$TotalGlycanComposition))

  df$GlycansPerSite <- sapply(df$GlycansPerSite, function(x) ifelse(x > 9, ">10", toString(x)))
  df$GlycansPerSite <- factor(df$GlycansPerSite, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", ">10")))

  df <- df %>%
    dplyr::summarise(.by = .data$GlycansPerSite, Count = dplyr::n())

  df <- df %>%
    dplyr::arrange(dplyr::desc(.data$GlycansPerSite)) %>%
    dplyr::mutate(prop = .data$Count / sum(.data$Count) * 100) %>%
    dplyr::mutate(ypos = cumsum(.data$prop) - 0.5 * .data$prop)

  p <- ggplot2::ggplot(df, ggplot2::aes(x="", y=.data$Count, fill=.data$GlycansPerSite)) +
    ggplot2::geom_bar(stat="identity", width=1, color= "black") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Count),
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::labs(fill = "Glycans\nper site")

  print(p + ggplot2::theme_void())
}
