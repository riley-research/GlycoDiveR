#' PlotGlycanSubcellularLocation
#'
#' Visualize the count or intensity of glycopeptides per subcellular localization.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param summaryFunction Either "count" to visualize the number of unique
#' glycopeptides, or "intensity" to summarise by the log2 intensity of the unique
#' glycopeptides.
#' @param zscoreCutoff An absolute z-score cutoff.
#' @param type Either "all" to plot for all glycopeptides, or specify a vector of
#' glycotypes as in the GlycanType column (e.g., c("Sialylated", "Fucosylated"))
#' @param pointSize What is the minimum and maximum point size.
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
#' @param silent silence printed information
#'
#' @returns a point plot.
#' @export
#'
#' @examples \dontrun{
#' PlotGlycanSubcellularLocation(mydata, summaryFunction = "intensity", zscoreCutoff = 2,
#' type = "Fucosylated", pointSize = c(1,10))
#' }
PlotGlycanSubcellularLocation <- function(input, summaryFunction = "count", zscoreCutoff = 1.5,
                                          type = "all", pointSize = c(3, 10), whichAlias = NULL,
                                          whichPeptide = NULL, whichProtein = NULL,
                                          exactProteinMatch = TRUE, silent = FALSE){
  summary_fun <- function(x = NULL) {
    if (summaryFunction == "count") {
      dplyr::n()
    } else {
      log(sum(x, na.rm = TRUE), 2)
    }
  }

  col_fun <- circlize::colorRamp2(
    breaks = c(-3, -2, -0.25, 0, 0.25, 2, 3),
    colors = c("#2166AC", "#88CCEE", "white", "white", "white", "#8877A1", "#882255")
  )

  legendTitle <- ifelse(summaryFunction == "count", "Glycopeptide count", "Glycopeptide intensity\n(log2)")

  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein,exactProteinMatch)
  input$PTMTable <- input$PTMTable %>% dplyr::filter(!is.na(.data$Intensity))

  df <- GetMeanTechReps(input$PTMTable) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if (length(type) == 1 && !identical(type, "all")) {
    if(!silent){
      fmessage(paste("Filtering for: ", type, collapse = ", "))
    }
    df <- df %>%
      dplyr::filter(.data$GlycanType %in% type)
  }

  if(nrow(df) == 0){return("Nothing left after filtering.")}

  df <- df %>%
    dplyr::filter(!is.na(.data$SubcellularLocalization) & .data$SubcellularLocalization != "") %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .keep_all = TRUE) %>%
    tidyr::separate_longer_delim(cols = "SubcellularLocalization", delim = ";")

  p <- df %>%
    dplyr::summarise(.by = c("Alias", "Condition", "SubcellularLocalization"),
                     summarized = summary_fun(.data$Intensity)) %>%
    dplyr::summarise(.by = c("Condition", "SubcellularLocalization"),
                     summarized = stats::median(.data$summarized)) %>%
    dplyr::mutate(.by = "SubcellularLocalization",
      z_score = (.data$summarized - mean(.data$summarized, na.rm = TRUE)) /
        stats::sd(.data$summarized, na.rm = TRUE),
      z_score = dplyr::if_else(is.na(.data$z_score), 0, .data$z_score),
      fill_color = col_fun(.data$z_score)) %>%
    dplyr::mutate(.by = "SubcellularLocalization",
                  toKeep = dplyr::if_else(max(abs(.data$z_score)) >= zscoreCutoff, "keep", "remove")) %>%
    dplyr::filter(.data$toKeep == "keep") %>%
    dplyr::arrange(.data$summarized) %>%
    dplyr::mutate(SubcellularLocalization = factor(.data$SubcellularLocalization, levels = unique(.data$SubcellularLocalization))) %>%
    ggplot2::ggplot(ggplot2::aes(x=.data$Condition,y= .data$SubcellularLocalization,
                                 size = .data$summarized, fill = .data$z_score)) +
    ggplot2::geom_point(stat="identity", shape = 21, color = "black") +
    ggplot2::labs(x=NULL, y=NULL, size = legendTitle) +
    ggplot2::scale_size(range = pointSize) +
    ggplot2::scale_fill_gradientn(
      colors = col_fun(seq(-3, 3, length.out = 100)),
      limits = c(-3, 3),
      oob = scales::squish,
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "+1", "+2", "+3")
    )

  return(p)
}
