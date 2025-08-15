#' PlotGlycanPerSite
#'
#' @param input Formatted data
#'
#' @returns The number of glycans per site
#' @export
#'
#' @examples \dontrun{PlotGlycanPerSite(myData)}
PlotGlycansPerSite <- function(input){
  input <- FilterForCutoffs(input)

  df <- GetMeanTechReps(input$PTMTable)

  df <- df %>%
    dplyr::filter(GlycanType != "NonGlyco") %>%
    dplyr::summarise(.by = c(UniprotIDs, ModificationID), GlycansPerSite = dplyr::n_distinct(TotalGlycanComposition))

  df$GlycansPerSite <- sapply(df$GlycansPerSite, function(x) ifelse(x > 9, ">10", toString(x)))
  df$GlycansPerSite <- factor(df$GlycansPerSite, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", ">10")))

  df <- df %>%
    dplyr::group_by() %>%
    dplyr::summarise(.by = GlycansPerSite, Count = n())

  df <- df %>%
    dplyr::arrange(desc(GlycansPerSite)) %>%
    dplyr::mutate(prop = Count / sum(Count) * 100) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop)

  p <- ggplot2::ggplot(df, aes(x="", y=Count, fill=GlycansPerSite)) +
    ggplot2::geom_bar(stat="identity", width=1, color= "black") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::scale_fill_manual(values = colorScheme, guide = guide_legend(reverse = TRUE)) +
    ggplot2::geom_text(aes(label = Count),
                       position = position_stack(vjust = 0.5)) +
    ggplot2::labs(fill = "Glycans\nper site")

  print(p + ggplot2::theme_void())
}
