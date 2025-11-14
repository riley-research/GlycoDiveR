#' PlotQuantificationQC
#'
#' @param input Formatted data
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param whichQuantification choose "both", "corrected", or "raw". Determines
#' what quantification values to plot.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param silent silence printed information (default = TRUE)
#'
#' @returns Log2 intensity boxplots
#' @export
#'
#' @examples \dontrun{PlotQuantificationQC(mydata, whichQuantification = "both")}
PlotQuantificationQC <- function(input, whichAlias = NULL, whichQuantification = "both",
                                 whichPeptide = NA, silent = FALSE){
  df <- FilterForPeptides(input$PSMTable, whichPeptide)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(whichQuantification == "corrected"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$Alias,y=log(.data$Intensity,2), fill = .data$Condition))+
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = colorScheme)
  }else if(whichQuantification == "raw"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$Alias,y=log(.data$RawIntensity,2), fill = .data$Condition))+
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = colorScheme)
  }else if(whichQuantification == "both"){
    p <- df %>%
      tidyr::pivot_longer(cols = c("Intensity", "RawIntensity"),
                          names_to = "QuantificationType", values_to = "Quantification") %>%
      ggplot2::ggplot(df, mapping = ggplot2::aes(x=.data$Alias,y=log(.data$Quantification,2), fill = .data$Condition))+
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = colorScheme) +
      ggplot2::facet_wrap(~.data$QuantificationType)
  }else{
    fmessage("Did not recognize whichQuantification parameter. Defaulting to 'both'")
    p <- df %>%
      tidyr::pivot_longer(cols = c("Intensity", "RawIntensity"),
                          names_to = "QuantificationType", values_to = "Quantification") %>%
      ggplot2::ggplot(df, mapping = ggplot2::aes(x=.data$Alias,y=log(.data$Quantification,2), fill = .data$Condition))+
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = colorScheme) +
      ggplot2::facet_wrap(~.data$QuantificationType)
  }

  return(p)
}
