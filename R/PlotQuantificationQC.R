#' PlotQuantificationQC
#'
#' A boxplot showing Log2 intensity values either before or after normalization, or both.
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
#' @returns Log2 intensity boxplots
#' @export
#'
#' @examples \dontrun{
#' PlotQuantificationQC(mydata, whichQuantification = "both")}
PlotQuantificationQC <- function(input, whichAlias = NULL, whichQuantification = "both",
                                 whichPeptide = NULL, whichProtein = NULL,
                                 exactProteinMatch = TRUE, silent = FALSE){
  df <- FilterForPeptides(input$PSMTable, whichPeptide)
  df <- FilterForProteins(df, whichProtein, exactProteinMatch)
  df <- df %>% dplyr::filter(!is.na(.data$Intensity))

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

  if(whichQuantification == "corrected"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$Alias,y=log(.data$Intensity,2), fill = .data$Condition))+
      ggplot2::geom_boxplot(outlier.fill = NULL, outlier.color = "grey25", outlier.shape = 21) +
      ggplot2::labs(x = NULL, y = "Intensity (log2)", title = "Normalized intensity") +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold", hjust=0.5))
  }else if(whichQuantification == "raw"){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$Alias,y=log(.data$RawIntensity,2), fill = .data$Condition))+
      ggplot2::geom_boxplot(outlier.fill = NULL, outlier.color = "grey25", outlier.shape = 21) +
      ggplot2::labs(x = NULL, y = "Intensity (log2)", title = "Raw intensity") +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold", hjust=0.5))
  }else if(whichQuantification == "both"){
    p <- df %>%
      tidyr::pivot_longer(cols = c("Intensity", "RawIntensity"),
                          names_to = "QuantificationType", values_to = "Quantification") %>%
      dplyr::mutate(QuantificationType = ifelse(.data$QuantificationType == "Intensity", "Normalized intensity", "Raw intensity")) %>%
      ggplot2::ggplot(df, mapping = ggplot2::aes(x=.data$Alias,y=log(.data$Quantification,2), fill = .data$Condition))+
      ggplot2::geom_boxplot(outlier.fill = NULL, outlier.color = "grey25", outlier.shape = 21) +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::facet_wrap(~.data$QuantificationType) +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face="bold", size = 12))
  }else{
    fmessage("Did not recognize whichQuantification parameter. Defaulting to 'both'")
    p <- df %>%
      tidyr::pivot_longer(cols = c("Intensity", "RawIntensity"),
                          names_to = "QuantificationType", values_to = "Quantification") %>%
      dplyr::mutate(QuantificationType = ifelse(.data$QuantificationType == "Intensity", "Normalized intensity", "Raw intensity")) %>%
      ggplot2::ggplot(df, mapping = ggplot2::aes(x=.data$Alias,y=log(.data$Quantification,2), fill = .data$Condition))+
      ggplot2::geom_boxplot(outlier.fill = NULL, outlier.color = "grey25", outlier.shape = 21) +
      ggplot2::labs(x = NULL, y = "Intensity (log2)") +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::facet_wrap(~.data$QuantificationType) +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face="bold", size = 12))
  }

  return(p)
}
