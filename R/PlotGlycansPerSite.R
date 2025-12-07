#' PlotGlycanPerSite
#'
#' This plot shows how many glycans were detected at each site—
#' that is, how many sites have 1 glycan, 2 glycans, and so on.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
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
#' @returns The number of glycans per site
#' @export
#'
#' @examples \dontrun{
#' PlotGlycanPerSite(myData)}
PlotGlycansPerSite <- function(input, whichAlias = NULL, whichPeptide = NULL,
                               whichProtein = NULL, exactProteinMatch = TRUE,
                               silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

  df <- GetMeanTechReps(input$PTMTable)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  df <- df %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::summarise(.by = c("UniprotIDs", "ModificationID"),
                     GlycansPerSite = dplyr::n_distinct(.data$TotalGlycanComposition))

  df$GlycansPerSite <- sapply(df$GlycansPerSite, function(x) ifelse(x > 9, ">10", toString(x)))
  df$GlycansPerSite <- factor(df$GlycansPerSite, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", ">10")))

  df <- df %>%
    dplyr::summarise(.by = "GlycansPerSite", Count = dplyr::n())

  df <- df %>%
    dplyr::arrange(dplyr::desc(.data$GlycansPerSite)) %>%
    dplyr::mutate(prop = .data$Count / sum(.data$Count) * 100) %>%
    dplyr::mutate(ypos = cumsum(.data$prop) - 0.5 * .data$prop)

  nColNeeded <- length(unique(as.character(df$GlycansPerSite)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x="", y=.data$Count, fill=.data$GlycansPerSite)) +
    ggplot2::geom_bar(stat="identity", linewidth=1.25, color= "white") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::scale_fill_manual(values = .modEnv$colorScheme[nColNeeded:1], guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Count),
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::labs(fill = "Glycans\nper site") +
    ggplot2::theme_void()

  return(p)
}
