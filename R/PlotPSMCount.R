#' PlotPSMCount
#'
#' @param input Formatted data
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition"
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param silent silence printed information (default = TRUE)
#'
#' @returns the PSM count
#' @export
#'
#' @examples \dontrun{PlotPSMCount(myData)}
PlotPSMCount <- function(input, grouping = "technicalReps", whichAlias = NULL,
                         whichPeptide = NA, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(grouping == "technicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c("Run", "Alias", "Glycan")) %>%
      dplyr::mutate(Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")),
                    Alias = factor(.data$Alias, levels = levels(input$PSMTable$Alias)))

    p <- ggplot2::ggplot(tempdf, ggplot2::aes(x=.data$Alias, y = .data$PSMCount, fill = .data$Glycan)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, NA)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    return(p)
  }

  if(grouping == "biologicalReps"){
    #Get separate points
    tempdf <- input$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c("Condition", "BioReplicate",
                               "TechReplicate", "Glycan")) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate),
                    Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(.data$PSMCount, na.rm = TRUE),
        sd   = stats::sd(.data$PSMCount, na.rm = TRUE),
        .by = c("x", "Glycan")
      ) %>%
      dplyr::mutate(.by = "x",
                    sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_mean, .data$mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")],
                       by = c("Glycan", "x")) %>%
      dplyr::mutate(.by = "x", mean = mean[which(.data$Glycan == "Glycosylated")[1]]) %>%
      dplyr::mutate(sum_count = .data$mean + .data$PSMCount) %>%
      dplyr::mutate(corrected_count = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_count, .data$PSMCount))

    #Get highest and lowest values to calculate y-axis limits
    minVal <- min(c(tempdfsum$corrected_mean - tempdfsum$sd, tempdf$corrected_count), na.rm = TRUE)
    maxVal <- max(c(tempdfsum$corrected_mean + tempdfsum$sd, tempdf$corrected_count), na.rm = TRUE)

    #Plot graph
    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$x, y = .data$mean, fill = .data$Glycan),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$x, ymin = .data$corrected_mean-.data$sd,
                                                            ymax = .data$corrected_mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) ggplot2::expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$x, y = .data$corrected_count, shape = .data$Glycan)) +
      ggplot2::scale_fill_manual(values = c(colorScheme)) +
      ggplot2::scale_shape_manual(values = c(15, 17))

    return(p)
  }
  if(grouping == "condition"){
    #Get separate points
    tempdf <- input$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c("Condition", "BioReplicate", "Glycan")) %>%
      dplyr::mutate(x = paste0(.data$Condition),
                    Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(.data$PSMCount, na.rm = TRUE),
        sd = stats::sd(.data$PSMCount, na.rm = TRUE),
        .by = c("x", "Glycan")
      ) %>%
      dplyr::mutate(.by = "x",
                    sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_mean, .data$mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")],
                       by = c("Glycan", "x")) %>%
      dplyr::mutate(.by = "x", mean = mean[which(.data$Glycan == "Glycosylated")[1]]) %>%
      dplyr::mutate(sum_count = .data$mean + .data$PSMCount) %>%
      dplyr::mutate(corrected_count = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_count, .data$PSMCount))

    #Get highest and lowest values to calculate y-axis limits
    minVal <- min(c(tempdfsum$corrected_mean - tempdfsum$sd, tempdf$corrected_count), na.rm = TRUE)
    maxVal <- max(c(tempdfsum$corrected_mean + tempdfsum$sd, tempdf$corrected_count), na.rm = TRUE)

    #Plot graph
    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, ggplot2::aes(x=.data$x, y = .data$mean, fill = .data$Glycan),
                        stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, ggplot2::aes(x = .data$x, ymin = .data$corrected_mean-.data$sd,
                                                            ymax = .data$corrected_mean+.data$sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) ggplot2::expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, ggplot2::aes(x=.data$x, y = .data$corrected_count, shape = .data$Glycan)) +
      ggplot2::scale_fill_manual(values = c(colorScheme)) +
      ggplot2::scale_shape_manual(values = c(15, 17))

    return(p)
  }

}
