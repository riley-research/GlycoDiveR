PlotGlycositesVsGlycans <- function(input, whichAlias = NULL,
                                    labelGlycositeCutoff = 0, labelGlycanCutoff = 0,
                                    maxOverlaps = 10){
  input <- FilterForCutoffs(input)

  if(!is.null(whichAlias)){
    input$PTMTable <- input$PTMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  #Prepare data
  df <- input$PTMTable %>%
    dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "") %>%
    dplyr::distinct(.data$UniprotIDs, .data$TotalGlycanComposition, .data$ModificationID, .keep_all = TRUE) %>%
    dplyr::summarise(.by = c(.data$UniprotIDs, .data$Genes, .data$NumberOfNSites),
                     count = dplyr::n())

  #Get labels
  df <- df %>%
    dplyr::mutate(label = dplyr::case_when(.data$NumberOfNSites > labelGlycositeCutoff &
                                      .data$count > labelGlycanCutoff ~ .data$Genes,
                                    TRUE ~ NA))
  # Plot
  df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$NumberOfNSites, y =.data$count)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(color = "#7b5799", alpha = 0.5) +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data$label), max.overlaps = 10,
                              fill = NA, label.size = NA) +
    ggplot2::annotate("text", x= 0.95 * max(df$NumberOfNSites), y= 1.1 * max(df$NumberOfNSites),
                      label= "y=x", color = "grey30") +
    ggplot2::labs(x = "Number of glycosites", y = "Number of glycans")

}
