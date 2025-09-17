#' PlotGlycoPeptideCount
#'
#' @param input Formatted data
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition"
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#'
#' @returns A graph showing the number of unique glycopeptides
#' @export
#'
#' @examples \dontrun{PlotGlycoPSMCount(mydata, grouping = "condition")}
PlotGlycopeptideCount <- function(input, grouping, whichAlias = NULL){
  input <- FilterForCutoffs(input)

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(grouping == "technicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select(.data$Run, .data$Alias, .data$Condition, .data$Glycan, .data$Genes) %>%
      dplyr::summarise(.by = c(.data$Run, .data$Alias, .data$Glycan, .data$Condition),
                       PSMCount = dplyr::n())

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    p <- ggplot2::ggplot(tempdf, ggplot2::aes(x=.data$Alias, y = .data$PSMCount, fill = .data$Condition)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.05)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "biologicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select(.data$Run, .data$Alias, .data$Condition, .data$BioReplicate,
                    .data$TechReplicate, .data$Glycan, .data$Genes) %>%
      dplyr::summarise(.by = c(.data$Alias, .data$Glycan, .data$Condition,
                               .data$BioReplicate, .data$TechReplicate),
                       PSMCount = dplyr::n()) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate))

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::summarise(.by = c(.data$Condition, .data$BioReplicate, .data$x),
                       mean = mean(.data$PSMCount, na.rm = TRUE),
                       sd = stats::sd(.data$PSMCount, na.rm = TRUE))

    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$x, y = .data$mean, fill = .data$Condition),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$x, ymin = .data$mean-.data$sd,
                                                            ymax = .data$mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$x, y = .data$PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "condition"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Run, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::select(.data$Run, .data$Alias, .data$Condition, .data$BioReplicate,
                    .data$Glycan, .data$Genes) %>%
      dplyr::summarise(.by = c(.data$Alias, .data$Glycan, .data$Condition, .data$BioReplicate),
                       PSMCount = dplyr::n())

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::summarise(.by = .data$Condition,
                       mean = mean(.data$PSMCount, na.rm = TRUE),
                       sd = stats::sd(.data$PSMCount, na.rm = TRUE))

    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$Condition, y = .data$mean, fill = .data$Condition),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$Condition,
                                                            ymin = .data$mean-.data$sd,
                                                            ymax = .data$mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$Condition, y = .data$PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else{
    warning("Unidentified grouping: ", grouping)
  }

}
