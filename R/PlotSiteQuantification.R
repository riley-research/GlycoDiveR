#' Barplot quantification per site
#'
#' @param input The input data as imported through one of the GlycoDiveR importers.
#' @param protein The protein defined in the UniprotIDs column.
#' @param site The site as in the ModificationID column.
#' @param cutoff An Intensity column, either in percentage or as an absolute number.
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
#' \dontrun{PlotSiteQuantification(inputdata, protein = "P04004", site = "N243", cutoff = 1e9)}
#'
#' \dontrun{PlotSiteQuantification(inputdata, protein = "P04004", site = "N243", cutoff = "10%")}
PlotSiteQuantification <- function(input, protein, site, cutoff = NA,
                                   whichAlias = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)

  df <- GetMeanTechReps(input$PTMTable)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df <- df %>%
    dplyr::filter(.data$UniprotIDs == protein & .data$ModificationID == site) %>%
    dplyr::summarise(.by = c("UniprotIDs", "ModificationID", "Condition",
                             "BioReplicate", "TechReplicate", "TotalGlycanComposition"),
                     Intensity = sum(.data$Intensity, na.rm = TRUE))

  if(nrow(df) == 0){
    if(!silent){
      fmessage("Nothing left after filtering")
    }
    return(NULL)
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
    ggplot2::scale_fill_manual(values = colorScheme)

  return(p)
}
