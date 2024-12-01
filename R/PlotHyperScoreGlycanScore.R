PlotHyperScoreGlycanScore <- function(filteredData){
  if(filteredData$searchEngine %in% c("MSFragger")){
    filteredData$PSMTable$col <- ifelse(filteredData$PSMTable$PeptideScore > filteredData$peptideScoreCutoff &
                                          filteredData$PSMTable$GlycanScore > filteredData$glycanScoreCutoff , "green", "black")
  }

  p <- ggplot2::ggplot(filteredData$PSMTable, aes(x = GlycanScore, y = PeptideScore)) +
    geom_point(color = filteredData$PSMTable$col) +
    facet_wrap(~Alias)

  print(p)
}
