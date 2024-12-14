PlotPSMCount <- function(inputData){
  inputData$PSMTable$Glycan <- sapply(inputData$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

  tempdf <- inputData$PSMTable %>%
    dplyr::select(Run, Alias, Glycan, Genes) %>%
    group_by(Run, Alias, Glycan) %>%
    reframe(Run = Run, Alias = Alias, PSMCount = n()) %>%
    distinct() %>%
    mutate(Glycan = factor(Glycan, levels = c("nonGlycosylated", "Glycosylated")))

  p <- ggplot(tempdf, aes(x=Alias, y = PSMCount, fill = Glycan)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "", y = "PSMs (count)")

  p
}
