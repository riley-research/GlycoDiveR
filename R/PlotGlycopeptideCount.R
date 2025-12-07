#' PlotGlycoPeptideCount
#'
#' Visualize the number of unique glycopeptides identified using a barplot.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition".
#' @param whichAlias provide a vector of Aliases to only select these aliases
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
#' @returns A graph showing the number of unique glycopeptides
#' @export
#'
#' @examples \dontrun{
#' PlotGlycoPSMCount(mydata)
#'
#' PlotGlycoPSMCount(mydata, grouping = "technicalReps",
#'                   whichAlias = c("Alias1", "Alias2"))}
PlotGlycopeptideCount <- function(input, grouping = "condition", whichAlias = NULL,
                                  whichPeptide = NULL, whichProtein = NULL,
                                  exactProteinMatch = TRUE, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)

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

  if(grouping == "technicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select("Run", "Alias", "Condition", "Glycan", "Genes") %>%
      dplyr::summarise(.by = c("Run", "Alias", "Glycan", "Condition"),
                       PSMCount = dplyr::n())

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    p <- ggplot2::ggplot(tempdf, ggplot2::aes(x=.data$Alias, y = .data$PSMCount, fill = .data$Condition)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "", y = "Unique glycopeptides (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.05)) +
      ggplot2::scale_fill_manual(values = c(.modEnv$colorScheme)) +
      ggplot2::theme(axis.ticks.x= ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(hjust = 2))

    return(p)
  }else if(grouping == "biologicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select("Run", "Alias", "Condition", "BioReplicate",
                    "TechReplicate", "Glycan", "Genes") %>%
      dplyr::summarise(.by = c("Alias", "Glycan", "Condition",
                               "BioReplicate", "TechReplicate"),
                       PSMCount = dplyr::n()) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate))

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::summarise(.by = c("Condition", "BioReplicate", "x"),
                       mean = mean(.data$PSMCount, na.rm = TRUE),
                       sd = stats::sd(.data$PSMCount, na.rm = TRUE))

    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$x, y = .data$mean, fill = .data$Condition),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$x, ymin = .data$mean-.data$sd,
                                                            ymax = .data$mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "Unique glycopeptides (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$x, y = .data$PSMCount),
                          size = 2) +
      ggplot2::scale_fill_manual(values = c(.modEnv$colorScheme)) +
      ggplot2::theme(axis.ticks.x= ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(hjust = 2))

    return(p)
  }else if(grouping == "condition"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select("Run", "Alias", "Condition", "BioReplicate",
                    "Glycan", "Genes") %>%
      dplyr::summarise(.by = c("Alias", "Glycan", "Condition", "BioReplicate"),
                       PSMCount = dplyr::n())

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::summarise(.by = "Condition",
                       mean = mean(.data$PSMCount, na.rm = TRUE),
                       sd = stats::sd(.data$PSMCount, na.rm = TRUE))

    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$Condition, y = .data$mean, fill = .data$Condition),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$Condition,
                                                            ymin = .data$mean-.data$sd,
                                                            ymax = .data$mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "Unique glycopeptides (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$Condition, y = .data$PSMCount),
                          size = 2) +
      ggplot2::scale_fill_manual(values = c(.modEnv$colorScheme)) +
      ggplot2::theme(axis.ticks.x= ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(hjust = 2))

    return(p)
  }else{
    warning("Unidentified grouping: ", grouping)
  }

}
