#' PlotGlycoPSMCount
#'
#' @param input Formatted data
#' @param grouping grouping is "technicalReps", "biologicalReps", or "condition"
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#'
#' @returns A GlycoPSM graph
#' @export
#'
#' @examples \dontrun{PlotGlycoPSMCount(mydata, grouping = "condition")}
PlotGlycoPSMCount <- function(input, grouping, whichAlias = NULL){
  input <- FilterForCutoffs(input)

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(Alias %in% whichAlias)
  }

  if(grouping == "technicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::summarise(PSMCount = n(), .by = c(Run, Alias, Condition, Glycan)) %>%
      dplyr::mutate(Alias = factor(Alias, levels = levels(input$PSMTable$Alias)))

    p <- ggplot2::ggplot(tempdf, aes(x=Alias, y = PSMCount, fill = Condition)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.05)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "biologicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::summarise(PSMCount = n(), .by = c(Alias, Glycan, Condition, BioReplicate, TechReplicate)) %>%
      dplyr::mutate(x = paste0(Condition, BioReplicate),
                    Alias = factor(Alias, levels = levels(input$PSMTable$Alias)))

    tempdfsum <- tempdf %>%
      dplyr::group_by(Condition, BioReplicate) %>%
      dplyr::reframe(x = x, mean = mean(PSMCount, na.rm = TRUE),
                     sd = sd(PSMCount, na.rm = TRUE)) %>%
      dplyr::distinct(Condition, BioReplicate, .keep_all = TRUE)

    minVal <- min(c(tempdfsum$mean - tempdfsum$sd, tempdf$PSMCount), na.rm = TRUE)
    maxVal <- max(c(tempdfsum$mean + tempdfsum$sd, tempdf$PSMCount), na.rm = TRUE)

    p <- ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=x, y = mean, fill = Condition), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = x, ymin = mean-sd, ymax = mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, aes(x=x, y = PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "condition"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::summarise(PSMCount = n(), .by=c(Alias, Glycan, Condition, BioReplicate)) %>%
      dplyr::mutate(Alias = factor(Alias, levels = levels(input$PSMTable$Alias)))

    tempdfsum <- tempdf %>%
      dplyr::group_by(Condition) %>%
      dplyr::reframe(Condition = Condition, mean = mean(PSMCount, na.rm = TRUE),
                     sd = sd(PSMCount, na.rm = TRUE)) %>%
      dplyr::distinct(Condition, .keep_all = TRUE)

    minVal <- min(c(tempdfsum$mean - tempdfsum$sd, tempdf$PSMCount), na.rm = TRUE)
    maxVal <- max(c(tempdfsum$mean + tempdfsum$sd, tempdf$PSMCount), na.rm = TRUE)

    p <- ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=Condition, y = mean, fill = Condition), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = Condition, ymin = mean-sd, ymax = mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand = if (minVal < 0) expansion(0.01, 0) else c(0, 0),
                                  limits = if (minVal < 0) c(NA, maxVal * 1.05) else c(0, maxVal * 1.05)) +
      ggplot2::geom_point(data = tempdf, aes(x=Condition, y = PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else{
    warning("Unidentified grouping: ", grouping)
  }

}
