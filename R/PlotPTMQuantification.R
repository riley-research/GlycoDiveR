PlotPTMQuantification <- function(input, gene){
  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", AssignedModifications)) %>%
    dplyr::filter(Genes == gene)

  if(nrow(df) > 0){
    #Generate lineplot with PTM annotation
    labeldf <- distinct(df[c("ProteinPTMLocalization", "ModificationID")])
    labeldf <- rbind(labeldf, data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
                                         "ModificationID" = c(1, df$ProteinLength[1])))

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = 1), aes(x = x, y = y), size = 2) +
      ggplot2::geom_point(data = df, aes(x= ProteinPTMLocalization, y = 1), color = "white") +
      ggrepel::geom_label_repel(data = labeldf, aes(x =ProteinPTMLocalization, y = 1, label = ModificationID), max.overlaps = Inf) +
      theme_void()

    if(input$quantAvailable){
      dfQuant <- GetMeanTechReps(df)

      hmdf <- dfQuant[c("Alias", "ModificationID", "Intensity", "ProteinPTMLocalization", "")]

      hmdf <- hmdf %>%
        dplyr::arrange(ProteinPTMLocalization) %>%
        dplyr::select(!ProteinPTMLocalization) %>%
        dplyr::group_by(Alias, ModificationID) %>%
        dplyr::mutate(Intensity = log(sum(Intensity, na.rm = TRUE)), 2) %>%
        dplyr::distinct(Alias, ModificationID, Intensity) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = Alias, values_from = Intensity)

      hmdf[is.na(hmdf)] <- 0

      mtrx <- data.matrix(hmdf[,2:ncol(hmdf)])
      rownames(mtrx) <- hmdf$ModificationID

      p2 <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mtrx, cluster_columns = FALSE, cluster_rows = FALSE)))


      `%/%` <- patchwork:::`/.ggplot`

      print(p1 %/% p2)

    }


  }else{message("No glyco on this protein: ", gene)}

}
