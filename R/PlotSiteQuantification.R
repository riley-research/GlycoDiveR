#' Barplot quantification per site
#'
#' A barplot showing site glycan quantification.
#'
#' @param input The input data as imported through one of the GlycoDiveR importers.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param site The site as in the ModificationID column.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param cutoff An Intensity column, either in percentage or as an absolute number.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param silent TRUE if you want info to be printed, FALSE if not
#'
#' @returns A graph.
#' @export
#'
#' @details The mean of technical replicates is calculated and the input data filtered for the user specific glycan and peptide cutoffs.
#'
#' @examples
#' \dontrun{
#' PlotSiteQuantification(mydata, whichProtein = "P04004", site = "N243", cutoff = 1e9)
#'
#' PlotSiteQuantification(mydata, whichProtein = "P04004", site = "N243", cutoff = "10%")}
PlotSiteQuantification <- function(input, whichProtein, site, whichPeptide = NULL,
                                   exactProteinMatch = TRUE, cutoff = NA,
                                   whichAlias = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

  df <- GetMeanTechReps(input$PTMTable)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df <- df %>%
    dplyr::filter(.data$ModificationID == site) %>%
    dplyr::summarise(.by = c("UniprotIDs", "ModificationID", "Condition",
                             "BioReplicate", "TechReplicate", "TotalGlycanComposition"),
                     Intensity = sum(.data$Intensity, na.rm = TRUE))

  if(CheckForQuantitativeValues(df$Intensity)){
    if(!silent){
      return(fmessage("No quantitative data found."))
    }else{
      return()
    }
  }

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if(!is.na(cutoff)){
    if(substr(cutoff, nchar(cutoff), nchar(cutoff)) == "%"){
      ctfp <- as.numeric(substr(cutoff, 1, nchar(cutoff) -1))
      dftemp <- df %>%
        dplyr::filter(.data$Intensity > (ctfp/100) * max(.data$Intensity, na.rm = TRUE))
      df <- df %>%
        dplyr::filter(.data$TotalGlycanComposition %in% dftemp$TotalGlycanComposition)
    }else{
      dftemp <- df %>%
        dplyr::filter(.data$Intensity > cutoff)
      df <- df %>%
        dplyr::filter(.data$TotalGlycanComposition %in% dftemp$TotalGlycanComposition)
    }
  }

  dfsum <- df %>%
    dplyr::summarise(.by = c("Condition", "TotalGlycanComposition"),
                     mean = mean(.data$Intensity, na.rm = TRUE),
                     sd = stats::sd(.data$Intensity, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(mean))

  dfsum$TotalGlycanComposition <- factor(dfsum$TotalGlycanComposition, levels = unique(dfsum$TotalGlycanComposition))

  p <- ggplot2::ggplot(data = dfsum) +
    ggplot2::geom_bar(data = dfsum, ggplot2::aes(x = .data$TotalGlycanComposition,
                                                 y = .data$mean, fill = .data$Condition),
                      stat = "identity", color = "black") +
    ggplot2::geom_errorbar(data = dfsum, ggplot2::aes(x = .data$TotalGlycanComposition,
                                                      ymin = .data$mean, ymax = .data$mean+.data$sd)) +
    ggplot2::geom_point(data = df, ggplot2::aes(x = .data$TotalGlycanComposition, y = .data$Intensity)) +
    ggplot2::labs(x = "", y = "Intensity (a.u.)") +
    ggplot2::guides(fill="none") +
    ggplot2::facet_wrap(~.data$Condition) +
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, max(df$Intensity)*1.05)) +
    ggplot2::scale_fill_manual(values = .modEnv$colorScheme)

  return(p)
}
