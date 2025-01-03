PlotPSMCount <- function(inputData){
  inputData$PSMTable$Glycan <- sapply(inputData$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  tempdf <- inputData$PSMTable %>%
    dplyr::select(Run, Alias, Glycan, Genes) %>%
    group_by(Run, Alias, Glycan) %>%
    reframe(Run = Run, Alias = Alias, PSMCount = n()) %>%
    distinct() %>%
    mutate(Glycan = factor(Glycan, levels = c("nonGlycosylated", "Glycosylated")))

  tempdf$Alias <- factor(tempdf$Alias, levels = levels(inputData$PSMTable$Alias))

  p <- ggplot(tempdf, aes(x=Alias, y = PSMCount, fill = Glycan)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    labs(x = "", y = "PSM (count)") +
    scale_y_continuous(expand=c(0,0), limits = c(0, NA)) +
    scale_fill_manual(values = c(colorScheme[1:2]))

  p <- p + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid = element_blank())

  print(p)
}
