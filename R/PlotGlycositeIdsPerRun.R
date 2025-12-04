#' PlotGlycositeIdsPerRun
#'
#' Visualize site coverage per protein for different runs.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param whichProtein Filter what protein to plot (make sure to only provide one).
#' These should be the IDs as presented in the UniprotIDs column in your GlycoDiveR data.
#' This can either be a dataframe with a UniprotIDs column, or a vector with the
#' UniprotID you want to keep.
#' @param exactProteinMatch When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' @param colorVec The color scheme used in the plot. Default = c("#6dc381", "pink", "#6761A8")
#' @param silent silence printed information (default = FALSE).
#'
#' @returns A plot displaying the glycosites identified across different runs
#' @export
#'
#' @examples \dontrun{
#' PlotGlycositeIdsPerRun(mydata, whichProtein = "P17047")}
PlotGlycositeIdsPerRun <- function(input, whichProtein = NULL, exactProteinMatch = TRUE,
                                   whichAlias = NULL, colorVec = c("#00394c", "#27b56e", "white"),
                                   silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

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

  df$GlycanIdentifier <- apply(df[,c("ModificationID", "TotalGlycanComposition")], 1, function(x) paste(x[1], x[2], sep = "-"))

  #Generate lineplot with PTM annotation
  labeldf <- dplyr::distinct(df[c("ProteinPTMLocalization", "ModificationID")])
  labeldf2 <- data.frame(ProteinPTMLocalization = c(1, df$ProteinLength[1]),
                         "ModificationID" = c(1, df$ProteinLength[1]))

  #If filter, use only the filtered levels, otherwise use all runs
  if(!is.null(whichAlias)){
    factors_to_use <- levels(input$PTMTable$Alias)
    factors_to_use <- factors_to_use[factors_to_use %in% df$Alias]

    foundMods <- dplyr::distinct(df[c("ProteinPTMLocalization", "Alias")]) %>%
      tidyr::complete(Alias = factor(factors_to_use, levels = factors_to_use)) %>%
      dplyr::mutate(y = (length(levels(.data$Alias)) + 1 - as.integer(.data$Alias)) / 2)
  }else{
    foundMods <- dplyr::distinct(df[c("ProteinPTMLocalization", "Alias")]) %>%
      tidyr::complete(Alias = factor(c(levels(df$Alias)), levels = c(levels(df$Alias)))) %>%
      dplyr::mutate(y = (length(levels(.data$Alias)) + 1 - as.integer(.data$Alias)) / 2)
  }

  yVal <- max(foundMods$y)+0.5

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_line(data = data.frame(x = seq(1,df$ProteinLength[1]), y = yVal),
                       ggplot2::aes(x = .data$x, y = .data$y), linewidth = 6, color = colorVec[1]) +
    ggplot2::geom_point(data = df, ggplot2::aes(x= .data$ProteinPTMLocalization, y = yVal),
                        fill = colorVec[2], color = colorVec[3], size = 8, shape = 21) +
    ggplot2::geom_point(data = foundMods, mapping = ggplot2::aes(x=.data$ProteinPTMLocalization, y = .data$y),
                        shape = 22, size = 5, na.rm = TRUE, fill = "grey50", color = "black") +
    ggplot2::geom_label(data = labeldf2, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                      y = yVal, label = .data$ModificationID), label.size = NA) +
    ggrepel::geom_label_repel(data = labeldf, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                           y = yVal, label = .data$ModificationID),
                              max.overlaps = Inf, nudge_y = 0.4, label.size = NA, fill = NA,
                              color=colorVec[2]) +
    ggplot2::theme_void() +
    ggplot2::scale_y_continuous(breaks = unique(foundMods$y), labels = unique(foundMods$Alias)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text())

  return(p)
}
