PlotGlycansPerProtein <- function(input, protein){
  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", AssignedModifications)) %>%
    dplyr::filter(GlycanType != "NonGlyco") %>%
    dplyr::filter(UniprotIDs == protein)

  df$GlycanIdentifier <- apply(df[,c("ModificationID", "TotalGlycanComposition")], 1, function(x) paste(x[1], x[2], sep = "-"))

  if(nrow(df) > 0){
    #Generate lineplot with PTM annotation
    labeldf <- distinct(df[c("ProteinPTMLocalization", "ModificationID")])
    #labeldf <- rbind(labeldf, data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
    #                                     "ModificationID" = c(1, df$ProteinLength[1])))
    labeldf2 <- data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
                           "ModificationID" = c(1, df$ProteinLength[1]))

    numberOfGroups <- levels(input$PTMTable$Alias)

    df

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = 1), aes(x = x, y = y), size = 6, color = "#009DDC") +
      ggplot2::geom_point(data = df, aes(x= ProteinPTMLocalization, y = 1),
                          fill = "pink", color = "#6761A8", size = 8, shape = 21) +
      ggplot2::geom_label(data = labeldf2, aes(x =ProteinPTMLocalization, y = 1, label = ModificationID), label.size = NA) +
      ggrepel::geom_label_repel(data = labeldf, aes(x =ProteinPTMLocalization, y = 1, label = ModificationID), max.overlaps = Inf) +
      theme_void()

    }
}

#PlotGlycansPerProtein(mydata, "P16144")
