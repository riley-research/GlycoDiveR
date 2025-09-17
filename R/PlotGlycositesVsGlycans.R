PlotGlycositesVsGlycans <- function(input, whichAlias = NULL){
  input <- FilterForCutoffs(input)

  if(!is.null(whichAlias)){
    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  input$PTMTable %>%
    dplyr::filter(!is.na(TotalGlycanComposition) & TotalGlycanComposition != "") %>%
    dplyr::distinct(UniprotIDs, TotalGlycanComposition, ModificationID, .keep_all = TRUE) %>%
    dplyr::summarise(.by = c(UniprotIDs, Genes, NumberOfNSites), count = dplyr::n()) %>%
    ggplot(mapping = ggplot2::aes(x=NumberOfNSites, y =count)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::geom_point(fill = "red", color = "blue", alpha = 0.8) +
    ggplot2::labs(x = "Number of glycosites", y = "Number of glycans")

}
