#' PlotGlycositeIdsPerRun
#'
#' @param input Formatted data
#' @param protein Protein listed as in the UniprotIDs column
#'
#' @returns A plot displaying the glycosites identified across different runs
#' @export
#'
#' @examples \dontrun{PlotGlycositeIdsPerRun(mydata, protein = "P17047")}
PlotGlycositeIdsPerRun <- function(input, protein){
  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", AssignedModifications)) %>%
    dplyr::filter(GlycanType != "NonGlyco") %>%
    dplyr::filter(UniprotIDs == protein)

  if(nrow(df) == 0){
    return(fmessage("No glycosylation found on the protein"))
  }

  df$GlycanIdentifier <- apply(df[,c("ModificationID", "TotalGlycanComposition")], 1, function(x) paste(x[1], x[2], sep = "-"))

  #Generate lineplot with PTM annotation
  labeldf <- distinct(df[c("ProteinPTMLocalization", "ModificationID")])
  labeldf2 <- data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
                         "ModificationID" = c(1, df$ProteinLength[1]))

  foundMods <- distinct(df[c("ProteinPTMLocalization", "Alias")]) %>%
    tidyr::complete(Alias = factor(c(levels(df$Alias)), levels = c(levels(df$Alias)))) %>%
    dplyr::mutate(y = (length(levels(Alias)) + 1 - as.integer(Alias)) / 2)

  yVal <- max(foundMods$y)+0.5

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = yVal),
                       aes(x = x, y = y), linewidth = 6, color = "#6dc381") +
    ggplot2::geom_point(data = df, aes(x= ProteinPTMLocalization, y = yVal),
                        fill = "pink", color = "#6761A8", size = 8, shape = 21) +
    ggplot2::geom_point(data = foundMods, mapping = aes(x=ProteinPTMLocalization, y = y),
                        shape = 22, size = 5, na.rm = TRUE, fill = "grey50", color = "black") +
    ggplot2::geom_label(data = labeldf2, aes(x =ProteinPTMLocalization, y = yVal, label = ModificationID), label.size = NA) +
    ggrepel::geom_label_repel(data = labeldf, aes(x =ProteinPTMLocalization, y = yVal, label = ModificationID),
                              max.overlaps = Inf, nudge_y = 0.4, label.size = NA, fill = NA) +
    ggplot2::theme_void() +
    ggplot2::scale_y_continuous(breaks = unique(foundMods$y), labels = unique(foundMods$Alias)) +
    ggplot2::theme(axis.text.y = element_text())

  print(p)
  return(p)
}
