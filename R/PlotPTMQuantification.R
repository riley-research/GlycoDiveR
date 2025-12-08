#' PlotPTMQuantification
#'
#' Visualize the glycan quantification of a single protein using a heatmap.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param whichProtein Filter what protein to plot (make sure to only provide one).
#' These should be the IDs as presented in the UniprotIDs column in your GlycoDiveR data.
#' This can either be a dataframe with a UniprotIDs column, or a vector with the
#' UniprotID you want to keep.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param normalization "none".
#' @param plotColors Defines the colors of the barplot.
#' Default: plotColors = c("#00394c", "#27b56e", "white").
#' @param heatmapColors Defines the colors of the heatmap.
#' Default: heatmapColors = c("white", "#27b56e", "#fdd835", "#d84315").
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param lineWidth defines the black line around the heatmap cells. Use NA
#' to remove the line
#' @param rowFontSize font size for the row labels
#' @param showRowNames set to TRUE or FALSE
#' @param silent silence printed information
#'
#' @returns The PTM quantification for a specific protein
#' @export
#'
#' @examples \dontrun{
#' PlotPTMQuantification(mydata, whichProtein = "P07361")
#'
#' PlotPTMQuantification(mydata, whichProtein = "P07361", rowFontSize = 10,
#'                       showRowNames = FALSE, lineWidth = 0)}
PlotPTMQuantification <- function(input, whichProtein = NULL, whichPeptide = NULL, normalization = "none",
                                  collapseTechReps = FALSE, plotColors = c("#BAA5CC", "#32006e", "white"),
                                  heatmapColors = c("white", "#88CCEE", "#8877A1", "#882255"),
                                  exactProteinMatch = TRUE, whichAlias = NULL, lineWidth = 2,
                                  rowFontSize = 12, showRowNames = TRUE, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", .data$AssignedModifications)) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if(CheckForQuantitativeValues(df$Intensity)){
    if(!silent){
      return(fmessage("No quantitative data found."))
    }else{
      return()
    }
  }

  df$GlycanIdentifier <- apply(df[,c("ModificationID", "TotalGlycanComposition")], 1, function(x) paste(x[1], x[2], sep = "-"))

  plotTitle = paste0(df$UniprotIDs[1], ";", df$Genes[1])

  #Generate lineplot with PTM annotation
  labeldf <- dplyr::distinct(df[c("ProteinPTMLocalization", "ModificationID")])

  labeldf2 <- data.frame(ProteinPTMLocalization = c(1 - (df$ProteinLength[1] * 0.03),
                                                    df$ProteinLength[1] + (df$ProteinLength[1] * 0.05)),
                         "ModificationID" = c(1, df$ProteinLength[1]))

  p1 <- ggplot2::ggplot(df) +
    ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = 1),
                       ggplot2::aes(x = .data$x, y = .data$y), linewidth = 6, color = plotColors[1]) +
    ggplot2::geom_point(data = df, ggplot2::aes(x= .data$ProteinPTMLocalization, y = 1),
                        fill = plotColors[2], color = plotColors[3], size = 8, shape = 21) +
    ggplot2::geom_label(data = labeldf2, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                      y = 1, label = .data$ModificationID),
                        linewidth = NA) +
    ggrepel::geom_label_repel(data = labeldf, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                           y = 1, label = .data$ModificationID),
                              max.overlaps = Inf, nudge_y = 0.01, label.size = NA, fill = NA,
                              color = plotColors[2], show.legend = FALSE, verbose = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(b = 10)) +
    ggplot2::coord_cartesian(clip = "off")

    if(collapseTechReps){
      dfQuant <- GetMeanTechReps(df)
      dfQuant$Alias <- droplevels(dfQuant$Alias)

      if(!is.null(whichAlias)){
        levels_mtrx <- levels(input$PTMTable$Alias)
        levels_mtrx <- levels_mtrx[levels_mtrx %in% df$Alias]
      }else{
        levels_mtrx <- levels(dfQuant$Alias)
      }
    }else{
      dfQuant <- df

      if(!is.null(whichAlias)){
        levels_mtrx <- levels(input$PTMTable$Alias)
        levels_mtrx <- levels_mtrx[levels_mtrx %in% df$Alias]
      }else{
        levels_mtrx <- levels(input$PTMTable$Alias)
      }
    }

    hmdf <- dfQuant[c("Alias", "GlycanIdentifier", "Intensity", "ProteinPTMLocalization")]

    hmdf <- rbind(data.frame(Alias = levels_mtrx,
                             GlycanIdentifier = "filler",
                             Intensity = NA,
                             ProteinPTMLocalization = NA),
                  hmdf)

    hmdf <- hmdf %>%
      dplyr::arrange("ProteinPTMLocalization") %>%
      dplyr::select(!"ProteinPTMLocalization") %>%
      dplyr::summarise(.by = c("Alias", "GlycanIdentifier"),
                       Intensity = log(sum(.data$Intensity, na.rm = TRUE), base = 2)) %>%
      tidyr::pivot_wider(names_from = "Alias", values_from = "Intensity")

    row_to_remove <- which(apply(hmdf, 1, function(row) "filler" %in% row))
    hmdf <- hmdf[-row_to_remove,]

    mtrx <- data.matrix(hmdf[,2:ncol(hmdf)])
    rownames(mtrx) <- hmdf$GlycanIdentifier
    mtrx <- mtrx[,levels_mtrx, drop = FALSE]

    mtrx[is.infinite(mtrx)] <- 0
    mtrx[is.na(mtrx)] <- 0

    row_indices <- rowSums(mtrx != 0) > 0
    mtrx <- mtrx[row_indices, , drop = FALSE]

    if(normalization == "ZScore"){
      mtrx <- t(apply(mtrx, 1, function(row) {
        nz <- row != 0
        m  <- mean(row[nz])
        s  <- stats::sd(row[nz])
        row[nz] <- (row[nz] - m) / s
        row
      }))
      legendTitle <- "Z-score normalized\nintensity"
      notFoundVal <- min(mtrx, na.rm = TRUE) - 100
      mtrx[is.na(mtrx)] <- notFoundVal
      mtrx[mtrx == 0] <- notFoundVal
    }else{
      legendTitle <- "log2 Intensity"
      notFoundVal <- 0
    }

    #Get the color scheme
    valuesInMtrx <- sort(unique(as.vector(mtrx)))
    valuesInMtrx <- valuesInMtrx[valuesInMtrx != notFoundVal]

    if(sum(valuesInMtrx != notFoundVal) > 1){
      lowestVal <- floor(valuesInMtrx[1])
      secondLowestVal <- valuesInMtrx[2]
      highestVal <- ceiling(valuesInMtrx[length(valuesInMtrx)])

      col_fun = circlize::colorRamp2(c(lowestVal-1, lowestVal - 0.0001 ,lowestVal, mean(c(lowestVal, highestVal), na.rm =TRUE), highestVal),
                                     c(heatmapColors[1], heatmapColors))

      #Get the legend right
      lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = legendTitle, direction = "horizontal",
                                   at = round(seq(lowestVal, highestVal, length.out= 4), 2))
    }else if(length(valuesInMtrx) == notFoundVal){
      col_fun = circlize::colorRamp2(c(-1, 1), c(heatmapColors[1], heatmapColors[1]))
    }else{
      col_fun <- circlize::colorRamp2(c(valuesInMtrx[1]-0.01, valuesInMtrx[1]), c(heatmapColors[1], heatmapColors[4]))
      #col_fun(seq(low, high, length.out = 5))

      lgd <- ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = "log2 Intensity",
        direction = "horizontal",
        at = round(seq(valuesInMtrx[1]-0.01, valuesInMtrx, length.out = 2), 2))
    }

    #Get the right split
    splitVec <- sapply(rownames(mtrx), function(x) strsplit(x, "-")[[1]][1])

    p2 <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mtrx, cluster_columns = FALSE, cluster_rows = TRUE,
                                                                           rect_gp = grid::gpar(col = "black", lwd = lineWidth),
                                                                           column_title = plotTitle,
                                                                           col = col_fun,
                                                                           show_heatmap_legend = FALSE,
                                                                           row_split = splitVec,
                                                                           cluster_row_slices = TRUE,
                                                                           row_names_gp = grid::gpar(fontsize = rowFontSize),
                                                                           show_row_names = showRowNames, row_title_rot = 0),
                                                   heatmap_legend_list = lgd, heatmap_legend_side = "top"))

    return(patchwork::wrap_plots(p1, p2) + patchwork::plot_layout(heights = c(1, 8)))

}
