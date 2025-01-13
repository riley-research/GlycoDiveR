PlotPTMQuantification <- function(input, protein){
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

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = 1), aes(x = x, y = y), size = 6, color = "#009DDC") +
      ggplot2::geom_point(data = df, aes(x= ProteinPTMLocalization, y = 1),
                          fill = "pink", color = "#6761A8", size = 8, shape = 21) +
      ggplot2::geom_label(data = labeldf2, aes(x =ProteinPTMLocalization, y = 1, label = ModificationID), label.size = NA) +
      ggrepel::geom_label_repel(data = labeldf, aes(x =ProteinPTMLocalization, y = 1, label = ModificationID), max.overlaps = Inf) +
      theme_void()

    if(input$quantAvailable){
      dfQuant <- GetMeanTechReps(df)

      hmdf <- dfQuant[c("Alias", "GlycanIdentifier", "Intensity", "ProteinPTMLocalization")]
      hmdf <- rbind(data.frame(Alias = levels(input$PTMTable$Alias),
                               GlycanIdentifier = "filler",
                               Intensity = NA,
                               ProteinPTMLocalization = NA),
                    hmdf)

      hmdf <- hmdf %>%
        dplyr::arrange(ProteinPTMLocalization) %>%
        dplyr::select(!ProteinPTMLocalization) %>%
        dplyr::group_by(Alias, GlycanIdentifier) %>%
        dplyr::mutate(Intensity = log(sum(Intensity, na.rm = TRUE)), 2) %>%
        dplyr::distinct(Alias, GlycanIdentifier, Intensity) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = Alias, values_from = Intensity)

      row_to_remove <- which(apply(hmdf, 1, function(row) "filler" %in% row))
      hmdf <- hmdf[-row_to_remove,]

      mtrx <- data.matrix(hmdf[,2:ncol(hmdf)])
      rownames(mtrx) <- hmdf$GlycanIdentifier
      mtrx <- mtrx[,levels(input$PTMTable$Alias), drop = FALSE]

      mtrx[is.infinite(mtrx)] <- 0
      mtrx[is.na(mtrx)] <- 0

      row_indices <- rowSums(mtrx != 0) > 0
      mtrx <- mtrx[row_indices, , drop = FALSE]

      #Get the color scheme
      valuesInMtrx <- sort(unique(as.vector(mtrx)))

      if(sum(valuesInMtrx != 0) > 1){
        lowestVal <- valuesInMtrx[1]
        secondLowestVal <- valuesInMtrx[2]
        highestVal <- valuesInMtrx[length(valuesInMtrx)]

        col_fun = circlize::colorRamp2(c(secondLowestVal-0.0001, secondLowestVal - 0.001 ,secondLowestVal, highestVal), c("darkgrey","darkgrey", "blue", "red"))
        col_fun(seq(-3, 3))

        mtrx[is.na(mtrx)] <- secondLowestVal - 0.001
        mtrx[is.infinite(mtrx)] <- secondLowestVal - 0.001
        mtrx[mtrx == 0] <- secondLowestVal - 0.001
      }else if(sum(valuesInMtrx == 0) == length(valuesInMtrx)){
        col_fun = circlize::colorRamp2(c(-1, 1), c("darkgrey", "darkgrey"))
      }else{
        lowestVal <- valuesInMtrx[1]
        secondLowestVal <- valuesInMtrx[2]

        col_fun = circlize::colorRamp2(c(lowestVal, secondLowestVal), c("darkgrey", "red"))
        col_fun(seq(-3, 3))
      }

      #Get the right split
      splitVec <- sapply(rownames(mtrx), function(x) strsplit(x, "-")[[1]][1])

      #Get the legend right
      lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "log2 Intensity", direction = "horizontal")

      p2 <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mtrx, cluster_columns = FALSE, cluster_rows = TRUE,
                                                                             rect_gp = grid::gpar(col = "black", lwd = 2),
                                                                             column_title = protein,
                                                                             col = col_fun,
                                                                             show_heatmap_legend = FALSE,
                                                                             row_split = splitVec,
                                                                             cluster_row_slices = TRUE),
                                                     heatmap_legend_list = lgd, heatmap_legend_side = "top"))

      `%/%` <- patchwork:::`/.ggplot`

      print(p1 %/% p2 + patchwork::plot_layout(heights = c(1,8)))
    }

  }else{message("No glyco on this protein: ", protein)}

}
