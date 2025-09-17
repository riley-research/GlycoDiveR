#' PlotPSMCount
#'
#' @param inputData Formatted data
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition"
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#'
#' @returns the PSM count
#' @export
#'
#' @examples \dontrun{PlotPSMCount(myData)}
PlotPSMCount <- function(inputData, grouping = "technicalReps", whichAlias = NULL){
  inputData <- FilterForCutoffs(inputData)

  inputData$PSMTable$Glycan <- sapply(inputData$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(!is.null(whichAlias)){
    inputData$PSMTable <- inputData$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(grouping == "technicalReps"){
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c(.data$Run, .data$Alias, .data$Glycan)) %>%
      dplyr::mutate(Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")),
                    Alias = factor(.data$Alias, levels = levels(inputData$PSMTable$Alias)))

    p <- ggplot2::ggplot(tempdf, ggplot2::aes(x=.data$Alias, y = .data$PSMCount, fill = .data$Glycan)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, NA)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }

  if(grouping == "biologicalReps"){
    #Get separate points
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c(.data$Condition, .data$BioReplicate,
                               .data$TechReplicate, .data$Glycan)) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate),
                    Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(.data$PSMCount, na.rm = TRUE),
        sd   = stats::sd(.data$PSMCount, na.rm = TRUE),
        .by = c(.data$x, .data$Glycan)
      ) %>%
      dplyr::mutate(.by = .data$x,
                    sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_mean, .data$mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")],
                       by = c("Glycan", "x")) %>%
      dplyr::mutate(.by = .data$x, mean = mean[which(.data$Glycan == "Glycosylated")[1]]) %>%
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

    print(p)
  }
  if(grouping == "condition"){
    #Get separate points
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = dplyr::n(),
                       .by = c(.data$Condition, .data$BioReplicate, .data$Glycan)) %>%
      dplyr::mutate(x = paste0(.data$Condition),
                    Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(.data$PSMCount, na.rm = TRUE),
        sd = stats::sd(.data$PSMCount, na.rm = TRUE),
        .by = c(.data$x, .data$Glycan)
      ) %>%
      dplyr::mutate(.by = .data$x,
                    sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = dplyr::if_else(.data$Glycan == "nonGlycosylated", .data$sum_mean, .data$mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")],
                       by = c("Glycan", "x")) %>%
      dplyr::mutate(.by = .data$x, mean = mean[which(.data$Glycan == "Glycosylated")[1]]) %>%
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

    print(p)
  }

}
