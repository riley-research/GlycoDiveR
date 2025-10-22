#' PlotComparisonVolcano
#'
#' @param input formatted comparison dataframe (e.g. using GlycoDiveR's
#' ImportMSstatsComparison function)
#' @param whichComparison which comparison as in the Label column
#' @param statisticalCutoff what statistical cutoff (default = 0.05)
#' @param log2FCCutoff what log2 fold-change cutoff (default = 1)
#' @param statistic choose "adjpvalue" or "pvalue"
#' @param whichLabel what to label in the plot.
#' choose "significant", "none", or supply a vector with the points to label.
#' For example c("HPT-149", "SPA3K-N271"). These are internally generated using
#' the comparison dataframe: paste(Proteins, ModificationID, sep = "-")
#' @param maxOverlaps the maximum number of overlapping labels in the volcanoplot
#'
#' @returns A volcanoplot
#' @export
#'
#' @examples \dontrun{PlotComparisonVolcano(comparison, whichComparison = "Sal-PBS",
#' statistic = "pvalue", statisticalCutoff = 0.01)}
PlotComparisonVolcano <- function(input, whichComparison,
                                  statisticalCutoff = 0.05, log2FCCutoff = 1,
                                  statistic = "adjpvalue", whichLabel = "significant",
                                  maxOverlaps = 10){

  if(statistic == "pvalue"){
    df <- input %>%
      dplyr::mutate(statistic = .data$pvalue)
    ylabel = "-Log10 (p-value)"
  }else if(statistic == "adjpvalue"){
    df <- input %>%
      dplyr::mutate(statistic = .data$adjpvalue)
    ylabel = "-Log10 (FDR)"
  }else{
    warning("Did not recognize statistic parameter. Defaulting to 'adjpvalue'")
    df <- input %>%
      dplyr::mutate(statistic = .data$adjpvalue)
    ylabel = "-Log10 (FDR)"
  }

  df <- df %>%
    dplyr::filter(.data$Label == whichComparison & !is.na(.data$statistic)) %>%
    dplyr::mutate(col = dplyr::case_when(.data$statistic < statisticalCutoff & .data$log2FC < -log2FCCutoff ~ colorScheme[1],
                           .data$statistic < statisticalCutoff & .data$log2FC > log2FCCutoff ~ colorScheme[2],
                           TRUE ~ "black"))

  if(identical(whichLabel, "significant")){
    df <- df %>%
      dplyr::mutate(plotLabel = ifelse(.data$col != "black",
                                   paste(.data$Proteins, .data$ModificationID, sep = "-"),
                                   NA_character_))
  }else if(identical(whichLabel, "none")){
    df <- df %>%
      dplyr::mutate(plotLabel = NA_character_)
  }else{
    df <- df %>%
      dplyr::mutate(plotLabel = ifelse(paste(.data$Proteins, .data$ModificationID, sep = "-") %in% whichLabel,
                                       paste(.data$Proteins, .data$ModificationID, sep = "-"),
                                       NA_character_))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$log2FC, y = -log(.data$statistic, 10))) +
    ggplot2::geom_point(color = df$col, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log(statisticalCutoff, 10), linetype="dashed", color="grey60") +
    ggplot2::geom_vline(xintercept = c(-log2FCCutoff, log2FCCutoff), linetype="dashed", color="grey60") +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data$plotLabel), fill = NA,
                              label.size = NA, max.overlaps = maxOverlaps) +
    ggplot2::labs(y = ylabel, x = "Log2 fold change", title = whichComparison) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(p)

}
