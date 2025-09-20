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
    dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "" &
                    GlycanType != "NonGlyco") %>%
    dplyr::mutate(.by = .data$UniprotIDs,
                  count = dplyr::n_distinct(.data$TotalGlycanComposition)) %>%
    dplyr::summarise(.by = c(.data$UniprotIDs, .data$count, .data$Genes),
                             numberOfGlycosites = dplyr::n_distinct(.data$ModificationID))

  #Get labels
  df <- df %>%
    dplyr::mutate(label = dplyr::case_when(.data$numberOfGlycosites > labelGlycositeCutoff &
                                      .data$count > labelGlycanCutoff ~ .data$Genes,
                                    TRUE ~ NA))
  # Plot
  df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$numberOfGlycosites, y =.data$count)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(color = "#7b5799", alpha = 0.5) +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data$label), max.overlaps = 10,
                              fill = NA, label.size = NA) +
    ggplot2::annotate("text", x= 0.95 * max(df$numberOfGlycosites), y= 1.1 * max(df$numberOfGlycosites),
                      label= "y=x", color = "grey30") +
    ggplot2::labs(x = "Number of glycosites", y = "Number of glycans")

}
