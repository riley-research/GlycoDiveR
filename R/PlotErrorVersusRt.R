#' PlotErrorVersusRt
#'
#' A scatter plot showing the retention time versus the m/z error. The three
#' dotted lines represent the median, and the median plus or minus one standard
#' deviation. The median and standard deviation calculations will use only points
#' within the specified upper and lower limits.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param type Choose "both", "glyco", or "nonGlyco".
#' @param upperLimit The upper m/z error limit.
#' @param lowerLimit The lower m/z error limit.
#' @param alpha The transparency of the dots. Set between 0 and 1.
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
#' @param silent silence printed information.
#'
#' @returns A ggplot graph.
#' @export
#'
#' @examples \dontrun{
#' PlotErrorVersusRt(mydata, type = "both")
#'
#' PlotErrorVersusRt(mydata, type = "glyco", lowerLimit = -20, upperLimit = 20)
#' }
PlotErrorVersusRt <- function(input, type = c("both", "glyco", "nonGlyco"),
                              upperLimit = NA, lowerLimit = NA, alpha = 0.5,
                              whichAlias = NULL, whichPeptide = NULL,
                              whichProtein = NULL, exactProteinMatch = TRUE,
                              silent = FALSE) {
  type <- match.arg(type)
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  input$PSMTable <- input$PSMTable %>% dplyr::filter(!is.na(.data$Intensity))

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(input$PSMTable) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if(type == "both"){
    colH <- c("Glycosylated" = .modEnv$colorScheme[1],
              "nonGlycosylated" = .modEnv$colorScheme[2])

    sdValues <- input$PSMTable$ppmError[!is.na(input$PSMTable$ppmError)]
    if (!identical(lowerLimit, NA)) {sdValues <- sdValues[sdValues >= lowerLimit]}
    if (!identical(upperLimit, NA)) {sdValues <- sdValues[sdValues <= upperLimit]}
    sd <- stats::sd(sdValues, na.rm = TRUE)
    medianVal <- stats::median(sdValues, na.rm = TRUE)
    sd <- c(medianVal - sd, medianVal + sd)
    sd <- sd[sd >= min(sdValues, na.rm = TRUE) & sd <= max(sdValues, na.rm = TRUE)]


    input$PSMTable %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$RetentionTime,
                                   y = .data$ppmError,
                                   fill = .data$Glycan)) +
      ggplot2::geom_hline(yintercept = c(medianVal, sd), color = "grey50", linetype = "dashed") +
      ggplot2::geom_point(color = scales::alpha("black", alpha), shape = 21, size = 3) +
      ggplot2::labs(x = "Retention time (min)", y = "m/z error (ppm)", fill = NULL) +
      ggplot2::scale_fill_manual(values = scales::alpha(colH, alpha)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05), limits = c(lowerLimit, upperLimit))
  }else if(type == "glyco") {
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated")

    sdValues <- input$PSMTable$ppmError[!is.na(input$PSMTable$ppmError)]
    if (!identical(lowerLimit, NA)) {sdValues <- sdValues[sdValues >= lowerLimit]}
    if (!identical(upperLimit, NA)) {sdValues <- sdValues[sdValues <= upperLimit]}
    sd <- stats::sd(sdValues, na.rm = TRUE)
    medianVal <- stats::median(sdValues, na.rm = TRUE)
    sd <- c(medianVal - sd, medianVal + sd)
    sd <- sd[sd >= min(sdValues, na.rm = TRUE) & sd <= max(sdValues, na.rm = TRUE)]

    input$PSMTable %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$RetentionTime,
                                   y = .data$ppmError)) +
      ggplot2::geom_hline(yintercept = c(medianVal, sd), color = "grey50", linetype = "dashed") +
      ggplot2::geom_point(color = scales::alpha("black", alpha),
                          shape = 21, size = 3,
                          fill = scales::alpha(.modEnv$colorScheme[1], alpha)) +
      ggplot2::labs(x = "Retention time (min)", y = "m/z error (ppm)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05), limits = c(lowerLimit, upperLimit))
  }else if (type == "nonGlyco") {
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "nonGlycosylated")

    sdValues <- input$PSMTable$ppmError[!is.na(input$PSMTable$ppmError)]
    if (!identical(lowerLimit, NA)) {sdValues <- sdValues[sdValues >= lowerLimit]}
    if (!identical(upperLimit, NA)) {sdValues <- sdValues[sdValues <= upperLimit]}
    sd <- stats::sd(sdValues, na.rm = TRUE)
    medianVal <- stats::median(sdValues, na.rm = TRUE)
    sd <- c(medianVal - sd, medianVal + sd)
    sd <- sd[sd >= min(sdValues, na.rm = TRUE) & sd <= max(sdValues, na.rm = TRUE)]

    input$PSMTable %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$RetentionTime,
                                   y = .data$ppmError)) +
      ggplot2::geom_hline(yintercept = c(medianVal, sd), color = "grey50", linetype = "dashed") +
      ggplot2::geom_point(color = scales::alpha("black", alpha),
                          shape = 21, size = 3,
                          fill = scales::alpha(.modEnv$colorScheme[2], alpha)) +
      ggplot2::labs(x = "Retention time (min)", y = "m/z error (ppm)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05), limits = c(lowerLimit, upperLimit))
  }

}
