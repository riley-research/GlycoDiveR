#' PlotPTMQuantification
#'
#' @param input Formatted data
#' @param protein What protein to quantify the PTM for
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param lineWidth defines the black line around the heatmap cells. Use NA
#' to remove the line
#' @param rowFontSize font size for the row labels
#' @param showRowNames set to TRUE or FALSE
#'
#' @returns The PTM quantification for a specific protein
#' @export
#'
#' @examples \dontrun{PlotPTMQuantification(mydata, "P07361")}
PlotPTMQuantification <- function(input, protein, whichAlias = NULL, lineWidth = 2,
                                  rowFontSize = 12, showRowNames = TRUE){
  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", .data$AssignedModifications)) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::filter(.data$UniprotIDs == protein)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  df$GlycanIdentifier <- apply(df[,c("ModificationID", "TotalGlycanComposition")], 1, function(x) paste(x[1], x[2], sep = "-"))

  if(nrow(df) > 0){
    #Generate lineplot with PTM annotation
    labeldf <- dplyr::distinct(df[c("ProteinPTMLocalization", "ModificationID")])
    #labeldf <- rbind(labeldf, data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
    #                                     "ModificationID" = c(1, df$ProteinLength[1])))
    labeldf2 <- data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
                                                                "ModificationID" = c(1, df$ProteinLength[1]))

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = 1),
                         ggplot2::aes(x = .data$x, y = .data$y), size = 6, color = "#009DDC") +
      ggplot2::geom_point(data = df, ggplot2::aes(x= .data$ProteinPTMLocalization, y = 1),
                          fill = "pink", color = "#6761A8", size = 8, shape = 21) +
      ggplot2::geom_label(data = labeldf2, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                        y = 1, label = .data$ModificationID),
                          label.size = NA) +
      ggrepel::geom_label_repel(data = labeldf, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                             y = 1, label = .data$ModificationID),
                                max.overlaps = Inf) +
      ggplot2::theme_void()

    if(input$quantAvailable){
      dfQuant <- GetMeanTechReps(df)

      #Get levels right
      if(!is.null(whichAlias)){
        levels_mtrx <- levels(input$PTMTable$Alias)
        levels_mtrx <- levels_mtrx[levels_mtrx %in% df$Alias]
      }else{
        levels_mtrx <- levels(input$PTMTable$Alias)
      }

      hmdf <- dfQuant[c("Alias", "GlycanIdentifier", "Intensity", "ProteinPTMLocalization")]

      hmdf <- rbind(data.frame(Alias = levels_mtrx,
                               GlycanIdentifier = "filler",
                               Intensity = NA,
                               ProteinPTMLocalization = NA),
                    hmdf)

      hmdf <- hmdf %>%
        dplyr::arrange(.data$ProteinPTMLocalization) %>%
        dplyr::select(!.data$ProteinPTMLocalization) %>%
        dplyr::summarise(.by = c(.data$Alias, .data$GlycanIdentifier),
                         Intensity = log(sum(.data$Intensity, na.rm = TRUE), base = 2)) %>%
        tidyr::pivot_wider(names_from = .data$Alias, values_from = .data$Intensity)

      row_to_remove <- which(apply(hmdf, 1, function(row) "filler" %in% row))
      hmdf <- hmdf[-row_to_remove,]

      mtrx <- data.matrix(hmdf[,2:ncol(hmdf)])
      rownames(mtrx) <- hmdf$GlycanIdentifier
      mtrx <- mtrx[,levels_mtrx, drop = FALSE]

      mtrx[is.infinite(mtrx)] <- 0
      mtrx[is.na(mtrx)] <- 0

      row_indices <- rowSums(mtrx != 0) > 0
      mtrx <- mtrx[row_indices, , drop = FALSE]

      #Get the color scheme
      valuesInMtrx <- sort(unique(as.vector(mtrx)))
      valuesInMtrx <- valuesInMtrx[valuesInMtrx != 0]

      if(sum(valuesInMtrx != 0) > 1){
        lowestVal <- floor(valuesInMtrx[1])
        secondLowestVal <- valuesInMtrx[2]
        highestVal <- ceiling(valuesInMtrx[length(valuesInMtrx)])

        col_fun = circlize::colorRamp2(c(lowestVal-1, lowestVal - 0.0001 ,lowestVal, mean(c(lowestVal, highestVal), na.rm =TRUE), highestVal),
                                       c("grey99","grey99", "#4575B4", "#FEE08B", "#D73027"))
        #c("#D3D3D3","#D3D3D3", "#4575B4", "#74C476", "#D73027"))


      }else if(length(valuesInMtrx) == 0){
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
      lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "log2 Intensity", direction = "horizontal",
                                   at = seq(lowestVal, highestVal, length.out= 4))

      p2 <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mtrx, cluster_columns = FALSE, cluster_rows = TRUE,
                                                                             rect_gp = grid::gpar(col = "black", lwd = lineWidth),
                                                                             column_title = protein,
                                                                             col = col_fun,
                                                                             show_heatmap_legend = FALSE,
                                                                             row_split = splitVec,
                                                                             cluster_row_slices = TRUE,
                                                                             row_names_gp = grid::gpar(fontsize = rowFontSize),
                                                                             show_row_names = showRowNames),
                                                     heatmap_legend_list = lgd, heatmap_legend_side = "top"))

      print(patchwork::wrap_plots(p1, p2) + patchwork::plot_layout(heights = c(1, 8)))
    }

  }else{message("No glyco on this protein: ", protein)}

}
