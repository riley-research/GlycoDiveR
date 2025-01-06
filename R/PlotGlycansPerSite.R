PlotGlycansPerSite <- function(input){
  input <- FilterForCutoffs(input)

  df <- GetMeanTechReps(input$PTMTable)

  df <- subset(df, GlycanType != "NonGlyco")

  df <- df %>%
    dplyr::group_by(UniprotIDs, ModificationID) %>%
    dplyr::reframe(UniprotIDs = UniprotIDs, ModificationID = ModificationID,
            GlycansPerSite = dplyr::n_distinct(TotalGlycanComposition)) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  df$GlycansPerSite <- sapply(df$GlycansPerSite, function(x) ifelse(x > 9, ">10", toString(x)))
  df$GlycansPerSite <- factor(df$GlycansPerSite, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", ">10"))

  df <- df %>%
    dplyr::group_by(GlycansPerSite) %>%
    dplyr::reframe(GlycansPerSite = GlycansPerSite, Count = n()) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  df <- df %>%
    dplyr::arrange(desc(GlycansPerSite)) %>%
    dplyr::mutate(prop = Count / sum(Count) * 100) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop)

  p <- ggplot2::ggplot(df, aes(x="", y=Count, fill=GlycansPerSite)) +
    ggplot2::geom_bar(stat="identity", width=1, color= "black") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::scale_fill_manual(values = colorScheme) +
    ggplot2::geom_text(aes(label = Count),
                       position = position_stack(vjust = 0.5)) +
    ggplot2::labs(fill = "Glycans\nper site")

  print(p + ggplot2::theme_void())
}
