#' PlotPSMCount
#'
#' @param inputData Formatted data
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition"
#'
#' @returns the PSM count
#' @export
#'
#' @examples PlotPSMCount(myData)
PlotPSMCount <- function(inputData, grouping = "technicalReps"){
  inputData <- FilterForCutoffs(inputData)

  inputData$PSMTable$Glycan <- sapply(inputData$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(grouping == "technicalReps"){
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = n(), .by = c(Run, Alias, Glycan)) %>%
      dplyr::mutate(Glycan = factor(Glycan, levels = c("nonGlycosylated", "Glycosylated")),
                    Alias = factor(Alias, levels = levels(inputData$PSMTable$Alias)))

    p <- ggplot(tempdf, aes(x=Alias, y = PSMCount, fill = Glycan)) +
      geom_bar(stat = "identity", position = "stack", color = "black") +
      labs(x = "", y = "PSM (count)") +
      scale_y_continuous(expand=c(0,0), limits = c(0, NA)) +
      scale_fill_manual(values = c(colorScheme))

    print(p)
  }

  if(grouping == "biologicalReps"){
    #Get separate points
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = n(), .by = c(Condition, BioReplicate, TechReplicate, Glycan)) %>%
      dplyr::mutate(x = paste0(Condition, BioReplicate),
                    Glycan = factor(Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(PSMCount, na.rm = TRUE),
        sd   = sd(PSMCount, na.rm = TRUE),
        .by = c(x, Glycan)
      ) %>%
      dplyr::mutate(.by = x, sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = if_else(Glycan == "nonGlycosylated", sum_mean, mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")], by = join_by(Glycan, x)) %>%
      dplyr::mutate(.by = x, mean = mean[which(Glycan == "Glycosylated")[1]]) %>%
      dplyr::mutate(sum_count = mean + PSMCount) %>%
      dplyr::mutate(corrected_count = if_else(Glycan == "nonGlycosylated", sum_count, PSMCount))

    #Get highest and lowest values to calculate y-axis limits
    minVal <- min(c(tempdfsum$corrected_mean - tempdfsum$sd, tempdf$corrected_count))
    maxVal <- max(c(tempdfsum$corrected_mean + tempdfsum$sd, tempdf$corrected_count))

    #Plot graph
    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=x, y = mean, fill = Glycan), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = x, ymin = corrected_mean-sd, ymax = corrected_mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, aes(x=x, y = corrected_count, shape = Glycan)) +
      ggplot2::scale_fill_manual(values = c(colorScheme)) +
      ggplot2::scale_shape_manual(values = c(15, 17))

    print(p)
  }
  if(grouping == "condition"){
    #Get separate points
    tempdf <- inputData$PSMTable %>%
      dplyr::summarise(PSMCount = n(), .by = c(Condition, BioReplicate, Glycan)) %>%
      dplyr::mutate(x = paste0(Condition),
                    Glycan = factor(Glycan, levels = c("nonGlycosylated", "Glycosylated")))

    #Calculate mean and CV
    tempdfsum <- tempdf %>%
      dplyr::summarise(
        mean = mean(PSMCount, na.rm = TRUE),
        sd   = sd(PSMCount, na.rm = TRUE),
        .by = c(x, Glycan)
      ) %>%
      dplyr::mutate(.by = x, sum_mean = sum(mean, na.rm = TRUE)) %>%
      dplyr::mutate(corrected_mean = if_else(Glycan == "nonGlycosylated", sum_mean, mean))

    #Add mean of glycosylated peaks to nonGlycosylated PSMcount values
    tempdf <- tempdf %>%
      dplyr::left_join(tempdfsum[c("x", "Glycan", "mean")], by = join_by(Glycan, x)) %>%
      dplyr::mutate(.by = x, mean = mean[which(Glycan == "Glycosylated")[1]]) %>%
      dplyr::mutate(sum_count = mean + PSMCount) %>%
      dplyr::mutate(corrected_count = if_else(Glycan == "nonGlycosylated", sum_count, PSMCount))

    #Get highest and lowest values to calculate y-axis limits
    minVal <- min(c(tempdfsum$corrected_mean - tempdfsum$sd, tempdf$corrected_count))
    maxVal <- max(c(tempdfsum$corrected_mean + tempdfsum$sd, tempdf$corrected_count))

    #Plot graph
    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=x, y = mean, fill = Glycan), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = x, ymin = corrected_mean-sd, ymax = corrected_mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, aes(x=x, y = corrected_count, shape = Glycan)) +
      ggplot2::scale_fill_manual(values = c(colorScheme)) +
      ggplot2::scale_shape_manual(values = c(15, 17))

    print(p)
  }

}
