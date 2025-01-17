PlotGlycoPSMCount <- function(input, grouping){
  input <- FilterForCutoffs(input)

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  if(grouping == "technicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::select(Run, Alias, Condition, Glycan, Genes) %>%
      group_by(Run, Alias, Glycan) %>%
      reframe(Run = Run, Alias = Alias, Condition = Condition, PSMCount = n()) %>%
      distinct()

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    p <- ggplot(tempdf, aes(x=Alias, y = PSMCount, fill = Condition)) +
      geom_bar(stat = "identity", position = "stack", color = "black") +
      labs(x = "", y = "PSM (count)") +
      scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.05)) +
      scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "biologicalReps"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::select(Run, Alias, Condition, BioReplicate, TechReplicate, Glycan, Genes) %>%
      dplyr::group_by(Alias, Glycan, Condition, BioReplicate, TechReplicate) %>%
      dplyr::reframe(Alias = Alias, Condition = Condition, BioReplicate = BioReplicate, TechReplicate = TechReplicate, PSMCount = n()) %>%
      dplyr::distinct() %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(x = paste0(Condition, BioReplicate))

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::group_by(Condition, BioReplicate) %>%
      dplyr::reframe(x = x, mean = mean(PSMCount, na.rm = TRUE),
                     sd = sd(PSMCount, na.rm = TRUE)) %>%
      dplyr::distinct(Condition, BioReplicate, .keep_all = TRUE)

    p <- ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=x, y = mean, fill = Condition), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = x, ymin = mean-sd, ymax = mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, aes(x=x, y = PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else if(grouping == "condition"){
    tempdf <- input$PSMTable %>%
      dplyr::filter(Glycan == "Glycosylated") %>%
      dplyr::select(Run, Alias, Condition, BioReplicate, Glycan, Genes) %>%
      dplyr::group_by(Alias, Glycan, Condition, BioReplicate) %>%
      dplyr::reframe(Alias = Alias, Condition = Condition, BioReplicate = BioReplicate, PSMCount = n()) %>%
      dplyr::distinct() %>%
      dplyr::ungroup() %>%
      dplyr::rowwise()

    tempdf$Alias <- factor(tempdf$Alias, levels = levels(input$PSMTable$Alias))

    tempdfsum <- tempdf %>%
      dplyr::group_by(Condition) %>%
      dplyr::reframe(Condition = Condition, mean = mean(PSMCount, na.rm = TRUE),
                     sd = sd(PSMCount, na.rm = TRUE)) %>%
      dplyr::distinct(Condition, .keep_all = TRUE)

    p <- ggplot() +
      ggplot2::geom_bar(data = tempdfsum, aes(x=Condition, y = mean, fill = Condition), stat = "identity", position = "stack", color = "black") +
      ggplot2::geom_errorbar(data = tempdfsum, aes(x = Condition, ymin = mean-sd, ymax = mean+sd), width = 0.2) +
      ggplot2::labs(x = "", y = "PSM (count)") +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = c(0, max(tempdf$PSMCount) * 1.10)) +
      ggplot2::geom_point(data = tempdf, aes(x=Condition, y = PSMCount)) +
      ggplot2::scale_fill_manual(values = c(colorScheme))

    print(p)
  }else{
    warning("Unidentified grouping: ", grouping)
  }

}
