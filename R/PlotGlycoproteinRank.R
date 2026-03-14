#' PlotGlycoproteinRank
#'
#' Plot the glycoprotein ranks using a scatter plot. The protein intensity is
#' the summed glycopeptide intensities of that proteins and then log10 transformed.
#' The histogram shows the total intensity of all glycoproteins per bin.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping Grouping is "all" "technicalReps", "biologicalReps",
#' or "condition".
#' @param bins The number of histogram bars.
#' @param pointSize The point size in the scatter plot.
#' @param geneLabel Specifu what genes to label. The values should correspond
#' to the values in your imported data's Genes column. Supply like this:
#' c("ORM1", "IGHM").
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param histogramColor The color of the histogram bars.
#' @param silent silence printed information
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples \dontrun{
#' PlotGlycoproteinRank(mydata, grouping = "condition")
#'
#' PlotGlycoproteinRank(mydata, pointSize = 6,
#' bins = 20, histogramColor = "#FFFFFF", grouping = "condition")
#' }
PlotGlycoproteinRank <- function(input, grouping = c("all", "technicalReps",
                                                     "biologicalReps", "condition"),
                                 bins = 10, pointSize = 4, geneLabel = NULL,
                                 whichProtein = NULL, whichPeptide = NULL,
                                 exactProteinMatch = TRUE, whichAlias = NULL,
                                 histogramColor = "grey70", silent = FALSE){
  grouping <- match.arg(grouping)
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  df <- FilterForPeptides(input$PSMTable, whichPeptide)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  #keep only glyco and remove 0 values, and keep one peptide per Alias
  df <- df %>%
    dplyr::mutate(PSMType = GetPSMGlycanCategory(.data$GlycanType)) %>%
    dplyr::filter(.data$PSMType != "nonGlyco" &
                    .data$Intensity != 0 &
                    !is.na(.data$Intensity)) %>%
    dplyr::arrange(dplyr::desc(.data$Intensity)) %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .keep_all = TRUE)


  if(grouping != "technicalReps"){
    df <- GetMeanTechReps(df)
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

  if (grouping == "all"){
    #Get the scatter plot
    df_plot <- df %>%
      dplyr::summarise(.by = c("UniprotIDs", "Genes"),
                       ProteinIntensity = log10(sum(.data$Intensity, na.rm = TRUE))) %>%
      dplyr::arrange(dplyr::desc(.data$ProteinIntensity)) %>%
      dplyr::mutate(Rank = dplyr::row_number(),
                    label = ifelse(.data$Genes%in% geneLabel, .data$Genes, NA))

    #Get the histogram
    df_hist <- df_plot %>%
      dplyr::mutate(Bin = cut(.data$ProteinIntensity, breaks = bins),
                    Bin_Mean = purrr::map_dbl(.data$Bin, function(x) {
                      nums <-
                        as.numeric(unlist(strsplit(gsub("\\(|\\]", "", as.character(x)), ",")))
                      mean(nums)})) %>%
      dplyr::summarise(.by = "Bin_Mean",
                       TotalProteinIntensity = sum(.data$ProteinIntensity))

    #Get the axis the same
    minY <- min(df_plot$ProteinIntensity, na.rm = TRUE)
    maxY <- max(df_plot$ProteinIntensity, na.rm = TRUE)
    bin_width_val <- (maxY - minY) / bins
    minY <- min(df_hist$Bin_Mean, na.rm = TRUE) - 0.55 * bin_width_val
    maxY <- max(df_hist$Bin_Mean, na.rm = TRUE) + 0.55 * bin_width_val

    #Now the plotting
    scatter_p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$Rank, y = .data$ProteinIntensity)) +
      ggplot2::geom_point(color = .modEnv$colorScheme[1], size = pointSize) +
      ggrepel::geom_text_repel(data = df_plot %>% dplyr::filter(!is.na(.data$label)),
                               ggplot2::aes(label = .data$label),
                               min.segment.length = 0, max.overlaps = Inf) +
      ggplot2::labs(x = "Protein rank", y = "Intensity (log10)") +
      ggplot2::scale_y_continuous(limits = c(minY, maxY)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5))

    hist_p <- ggplot2::ggplot(df_hist, ggplot2::aes(x = .data$Bin_Mean, y = .data$TotalProteinIntensity)) +
      ggplot2::geom_bar(stat = "identity", color = "black",
                        fill = histogramColor, width = bin_width_val) +
      ggplot2::labs(x = NULL, y = "Summed intensity") +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_x_continuous(limits = c(minY, maxY)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5),
                     plot.margin = ggplot2::margin(t = 5, r = 20, b = 5, l = 0))

    return(patchwork::wrap_plots(scatter_p, hist_p, nrow = 1, widths = c(6, 1)))
  }else if (grouping == "technicalReps") {
    df_col <- df %>%
      dplyr::arrange(.data$Alias) %>%
      dplyr::summarise(.by = "Condition", Alias = list(.data$Alias)) %>%
      dplyr::mutate(col = .modEnv$colorScheme[1:dplyr::n()]) %>%
      tidyr::unnest(cols = "Alias") %>%
      dplyr::distinct(.data$Alias, .data$col) %>%
      dplyr::rename("group" = "Alias")

    df_plot <- df %>%
      dplyr::mutate(group = .data$Alias) %>%
      dplyr::summarise(.by = c("UniprotIDs", "Genes", "group"),
                       ProteinIntensity = log10(sum(.data$Intensity, na.rm = TRUE))) %>%
      dplyr::arrange(.by = "group",
                     dplyr::desc(.data$ProteinIntensity)) %>%
      dplyr::mutate(.by = "group",
                    Rank = dplyr::row_number(),
                    label = ifelse(.data$Genes%in% geneLabel, .data$Genes, NA)) %>%
      dplyr::left_join(df_col, dplyr::join_by("group"))

  }else if (grouping == "biologicalReps") {
    df_col <- df %>%
      dplyr::mutate(group = paste0(.data$Condition, "-", .data$BioReplicate)) %>%
      dplyr::arrange(.data$Alias) %>%
      dplyr::mutate(group = factor(.data$group, levels = unique(.data$group))) %>%
      dplyr::summarise(.by = "Condition", group = list(.data$group)) %>%
      dplyr::mutate(col = .modEnv$colorScheme[1:dplyr::n()]) %>%
      tidyr::unnest(cols = "group") %>%
      dplyr::distinct(.data$group, .data$col)

    df_plot <- df %>%
      dplyr::mutate(group = paste0(.data$Condition, "-", .data$BioReplicate)) %>%
      dplyr::arrange(.data$Alias) %>%
      dplyr::mutate(group = factor(.data$group, levels = unique(.data$group))) %>%
      dplyr::summarise(.by = c("UniprotIDs", "Genes","group"),
                       ProteinIntensity = log10(sum(.data$Intensity, na.rm = TRUE))) %>%
      dplyr::arrange(.by = "group",
                     dplyr::desc(.data$ProteinIntensity)) %>%
      dplyr::mutate(.by = "group",
                    Rank = dplyr::row_number(),
                    label = ifelse(.data$Genes%in% geneLabel, .data$Genes, NA)) %>%
      dplyr::left_join(df_col, dplyr::join_by("group"))

  }else if (grouping == "condition") {
    df_col <- df %>%
      dplyr::mutate(group = .data$Condition) %>%
      dplyr::arrange(.data$Alias) %>%
      dplyr::mutate(group = factor(.data$group, levels = unique(.data$group))) %>%
      dplyr::summarise(.by = "Condition", group = list(.data$group)) %>%
      dplyr::mutate(col = .modEnv$colorScheme[1:dplyr::n()]) %>%
      tidyr::unnest(cols = "group") %>%
      dplyr::distinct(.data$group, .data$col)

    df_plot <- df %>%
      dplyr::mutate(group = .data$Condition) %>%
      dplyr::arrange(.data$Alias) %>%
      dplyr::mutate(group = factor(.data$group, levels = unique(.data$group))) %>%
      dplyr::summarise(.by = c("UniprotIDs", "Genes", "group"),
                       ProteinIntensity = log10(sum(.data$Intensity, na.rm = TRUE))) %>%
      dplyr::arrange(.by = "group",
                     dplyr::desc(.data$ProteinIntensity)) %>%
      dplyr::mutate(.by = "group",
                    Rank = dplyr::row_number(),
                    label = ifelse(.data$Genes%in% geneLabel, .data$Genes, NA)) %>%
      dplyr::left_join(df_col, dplyr::join_by("group"))
  }

  #Get the axis and breaks the same
  minY <- min(df_plot$ProteinIntensity, na.rm = TRUE)
  maxY <- max(df_plot$ProteinIntensity, na.rm = TRUE)
  predefined_breaks <- seq(minY, maxY, length.out = bins + 1)
  bin_width_val <- (maxY - minY) / bins
  minY <- mean(predefined_breaks[1:2]) - 0.55 * bin_width_val
  maxY <- mean(predefined_breaks[bins:bins+1]) + 0.55 * bin_width_val

  max_hist <- df_plot %>%
    dplyr::mutate(.by = "group",
                  Bin = cut(.data$ProteinIntensity, breaks = predefined_breaks, include.lowest = TRUE),
                  Bin_Mean = purrr::map_dbl(.data$Bin, function(x) {
                    nums <-
                      as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", as.character(x)), ",")))
                    mean(nums)})) %>%
    dplyr::summarise(.by = c("group", "Bin_Mean"),
                     TotalProteinIntensity = sum(.data$ProteinIntensity)) %>%
    dplyr::slice_max(.data$TotalProteinIntensity) %>%
    dplyr::pull(.data$TotalProteinIntensity)
  max_hist <- 10^floor(log10(max_hist)) * ceiling(max_hist / 10^floor(log10(max_hist)))

  allGroups <- unique(levels(droplevels(df_plot$group)))
  numOfGroups <- length(allGroups)
  plotList <- vector("list", numOfGroups * 2)

  for (i in seq_along(1:numOfGroups)) {
    df_this_plot <- df_plot %>%
      dplyr::filter(.data$group == allGroups[i])

    #Get the histogram
    df_hist <- df_this_plot %>%
      dplyr::mutate(Bin = cut(.data$ProteinIntensity, breaks = predefined_breaks, include.lowest = TRUE),
                    Bin_Mean = purrr::map_dbl(.data$Bin, function(x) {
                      nums <-
                        as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", as.character(x)), ",")))
                      mean(nums)})) %>%
      dplyr::summarise(.by = "Bin_Mean",
                       TotalProteinIntensity = sum(.data$ProteinIntensity))

    #Now the plotting
    scatter_p <- ggplot2::ggplot(df_this_plot, ggplot2::aes(x = .data$Rank, y = .data$ProteinIntensity)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$col),
                         size = pointSize) +
      ggrepel::geom_text_repel(data = df_this_plot %>% dplyr::filter(!is.na(.data$label)),
                               ggplot2::aes(label = .data$label),
                               min.segment.length = 0, max.overlaps = Inf) +
      ggplot2::scale_color_identity() +
      ggplot2::labs(x = "Protein rank", y = "Intensity (log10)",
                    title = allGroups[i]) +
      ggplot2::scale_y_continuous(limits = c(minY, maxY)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5),
                     plot.title = ggplot2::element_text(hjust = 0.5))

    hist_p <- ggplot2::ggplot(df_hist, ggplot2::aes(x = .data$Bin_Mean, y = .data$TotalProteinIntensity)) +
      ggplot2::geom_bar(stat = "identity", color = "black",
                        fill = histogramColor, width = bin_width_val) +
      ggplot2::labs(x = NULL, y = "Summed intensity") +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_x_continuous(limits = c(minY, maxY)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max_hist),
                                  breaks = c(0, max_hist)) +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5,
                                                         vjust = 0.5),
                     plot.margin = ggplot2::margin(t = 5, r = 20, b = 5, l = 0))

    plotList[[i*2-1]] <- scatter_p
    plotList[[i*2]] <- hist_p
  }

  #Get the dimensions, then wrap
  col <- floor(sqrt(numOfGroups)) * 2

  patchwork::wrap_plots(plotList, ncol = col, widths = rep(c(6, 1), col / 2))
}
