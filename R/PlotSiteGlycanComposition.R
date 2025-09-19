#' PlotSiteGlycanComposition
#'
#' @param input your formatted data
#' @param protein The protein you wish to plot, corresponding to the UniprotIDs
#' column of your dataframe
#' @param whichAlias a vector of aliases you want to include, corresponding
#' to the aliases in your annotation input
#' @param horizontalPoints specify the point size. Larger number is smaller points.
#' @param yCorrection move the points closer or further away. smaller values is
#' less spacing between the points
#' @param yNudge nudge the points up
#' @param boxSpacing defines the spacing of the box around the points. Larger numbers
#' mean more spacing
#' @param silent TRUE if you want info to be printed, FALSE if not
#'
#' @returns A plot showing glyco site
#' @export
#'
#' @examples \dontrun{PlotSiteGlycanComposition(input = mydata, protein = "P01042",
#' whichAlias = c("firstSample", "secondSample"))
#' }
PlotSiteGlycanComposition <- function(input, protein, whichAlias = NULL,
                                      horizontalPoints = 50, yCorrection = 0.25,
                                      yNudge = 1, boxSpacing = 0.05, silent = FALSE){

  set.seed(1)

  input <- FilterForCutoffs(input, silent)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", .data$AssignedModifications)) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::filter(.data$UniprotIDs == protein)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No glycosylation found on the protein"))
    }else{
      return()
    }

  }

  plotTitle <- paste0(protein, " - ", df[df$UniprotIDs == protein, "Genes"][[1]][1])

  df_sum <- df %>%
    dplyr::distinct(.data$ModificationID, .data$TotalGlycanComposition, .keep_all = TRUE) %>%
    dplyr::summarise(.by = c(.data$ModificationID, .data$GlycanType, .data$Domains,
                             .data$ProteinLength, .data$ProteinPTMLocalization),
                     count = dplyr::n())

  #Get the protein plot if there is domain information####
  linePlot_df <- df_sum[1,] %>%
    tidyr::separate_longer_delim(cols = "Domains", delim = ";") %>%
    dplyr::mutate(Domain = sub("\\(.*", "", .data$Domains),
                  FirstAA = as.integer(sub(".*\\(([^-]+)-.*", "\\1", .data$Domains)),
                  LastAA = as.integer(sub(".*-(.*)\\)", "\\1", .data$Domains))) %>%
    dplyr::select(.data$Domain, .data$FirstAA, .data$LastAA) %>%
    dplyr::bind_rows(data.frame(Domain = "Peptide", FirstAA = 1, LastAA = df_sum$ProteinLength[1])) %>%
    dplyr::mutate(size = ifelse(.data$Domain == "Peptide", 5, 7),
                  centerOfAA = (.data$FirstAA + .data$LastAA) / 2) %>%
    dplyr::arrange(dplyr::desc(.data$Domain == "Peptide"))

  linePlot_df <- linePlot_df[stats::complete.cases(linePlot_df),]

  #Set colors####
  unique_domains <- unique(linePlot_df$Domain)
  colors <- stats::setNames(c("#00394c", "#91b5c5", "#a9a9a9", "#22031F", "#433E0E" )[1:length(unique_domains)], unique_domains)
  linePlot_df <- linePlot_df %>%
    dplyr::mutate(color = colors[.data$Domain])

  #Get the labeling right
  text_df <- linePlot_df[linePlot_df$Domain != "Peptide",] %>%
    dplyr::sample_frac(1) %>%
    dplyr::distinct(.data$Domain, .keep_all=TRUE)
  text_df <- text_df[stats::complete.cases(text_df),]

  #Now the dots####
  #For each site calculate the boxes
  point_df <- df %>%
    dplyr::distinct(.data$ProteinPTMLocalization,
                    .data$ModificationID,
                    .data$TotalGlycanComposition,
                    .data$GlycanType,
                    .data$ProteinPTMLocalization) %>%
    dplyr::summarise(.by = c(.data$ModificationID, .data$GlycanType, .data$ProteinPTMLocalization),
                     count = dplyr::n()) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  colCount = ceiling(sqrt(sum(.data$count))),
                  count = purrr::map_chr(.data$count, ~ paste(1:.x, collapse = ";"))) %>%
    tidyr::separate_longer_delim(cols = "count", delim = ";") %>%
    dplyr::arrange(.data$GlycanType) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  count = as.numeric(.data$count),
                  cumCount = 1:length(.data$count),
                  x = ((.data$cumCount - 1) %% .data$colCount) + 1,
                  y = ceiling(.data$cumCount / .data$colCount))

    #Make sure points are top centered
    maxPointY <- max(point_df$y, na.rm = TRUE)
    point_df <- point_df %>%
      dplyr::mutate(.by = .data$ModificationID,
                    y = (maxPointY - max(.data$y)) + .data$y,
                    y = min(.data$y) + max(.data$y) - .data$y)

  #Value to use to determine location of the points
  pointSpacing <- df_sum$ProteinLength[1] / horizontalPoints
  maxY <- max(point_df$y, na.rm = TRUE)

  #Add the glycans as dots and glycan labels####
  labeldf <- dplyr::distinct(df_sum[c("ProteinPTMLocalization", "ModificationID")])

  #Now calculate the right spacing
  point_df <- point_df %>%
    dplyr::mutate(x_corrected = .data$ProteinPTMLocalization - (pointSpacing * (.data$colCount - 1) * 0.5) +
                    (.data$x - 1) * pointSpacing,
                  y_corrected = (.data$y * yCorrection) + (yNudge * (maxY - max(.data$y * yCorrection))))

  #Check for any clashes in point####
  #Start with the second group Check for any clashes with the previous group
  idVec <- point_df %>%
    dplyr::arrange(.data$ProteinPTMLocalization) %>%
    dplyr::distinct(.data$ModificationID) %>%
    dplyr::pull(.data$ModificationID)

  point_df$Corrected <- "No"

  if(length(idVec) > 1){
    for(i in 2:length(idVec)){
      maxGroup1 <- point_df %>%
        dplyr::filter(.data$ModificationID == idVec[i - 1]) %>%
        dplyr::pull(.data$x_corrected) %>%
        max()

      minGroup2 <- point_df %>%
        dplyr::filter(.data$ModificationID == idVec[i]) %>%
        dplyr::pull(.data$x_corrected) %>%
        min()

      if(minGroup2 - maxGroup1 < 2 * pointSpacing){
        correctionFactor <- (maxGroup1 - minGroup2) + boxSpacing * df_sum$ProteinLength[1]

        point_df <- point_df %>%
          dplyr::mutate(x_corrected = ifelse(.data$ModificationID == idVec[i],
                                             .data$x_corrected + correctionFactor,
                                             .data$x_corrected),
                        Corrected = ifelse(.data$ModificationID == idVec[i],
                                           "Yes", .data$Corrected))
      }
    }
  }

  #Get the point colors based on the glycan composition####
  point_df <- point_df %>%
    dplyr::mutate(glycanColor = dplyr::case_when(.data$GlycanType == "Complex/Hybrid" ~ "#00394a",
                                                 .data$GlycanType == "Sialyl+Fucose" ~ "#ff7f2a",
                                                 .data$GlycanType == "Sialyl" ~ "#2475b5",
                                                 .data$GlycanType == "Fucose" ~ "#aaaaaa",
                                                 .data$GlycanType == "High Mannose" ~ "#28b36d",
                                                 .data$GlycanType == "Truncated" ~ "#D0A5C0",
                                                 .data$GlycanType == "Paucimannose" ~ "#8B1E3F"))

  #Where do we want the boxes####
  rect_df <- point_df %>%
    dplyr::mutate(.by = .data$ModificationID,
                  xmin = min(.data$x_corrected, na.rm = TRUE) - boxSpacing * df_sum$ProteinLength[1] * 0.25,
                  xmax = max(.data$x_corrected, na.rm = TRUE) + boxSpacing * df_sum$ProteinLength[1] * 0.25,
                  ymin = min(.data$y_corrected, na.rm = TRUE) - boxSpacing * maxY,
                  ymax = max(.data$y_corrected, na.rm = TRUE) + boxSpacing * maxY) %>%
    dplyr::distinct(.data$xmin, .data$xmax, .data$ymin, .data$ymax)

  #Lines to the boxes####
  #If box was not adjusted due to overlap use a straight
  line_df <- point_df %>%
    dplyr::summarise(.by = c(.data$ModificationID, .data$ProteinPTMLocalization,
                             .data$Corrected),
                  maximumY = max(.data$y_corrected, na.rm = TRUE),
                  boxX = ifelse(.data$Corrected == "No", .data$ProteinPTMLocalization,
                                mean(c(min(.data$x_corrected), max(.data$x_corrected)))),
                  pointY = maxY + 1) %>%
    dplyr::distinct() %>%
    dplyr::mutate(xCords = paste(.data$ProteinPTMLocalization, .data$boxX, sep = ";"),
                  yCords = paste(.data$pointY, .data$maximumY, sep = ";")) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  xCords = ifelse(.data$Corrected == "No", .data$xCords,
                                  calculateElbowCoords(xVec = .data$xCords, yVec = .data$yCords, return = "x")),
                  yCords = ifelse(.data$Corrected == "No", .data$yCords,
                                  calculateElbowCoords(xVec = .data$xCords, yVec = .data$yCords, return = "y"))) %>%
    tidyr::separate_longer_delim(cols = c("xCords", "yCords"), delim = ";") %>%
    dplyr::mutate(xCords = as.double(.data$xCords),
                  yCords = as.double(.data$yCords))

  #colormapping####
  color_map <- stats::setNames(point_df$glycanColor, point_df$GlycanType)

  #Now plot the actual plot####
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = line_df, ggplot2::aes(x = .data$xCords, y = .data$yCords, group = .data$ModificationID),
                       color = "#27b56e", linewidth = 0.65, show.legend = FALSE) +
    ggplot2::geom_rect(data = rect_df, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                                    ymin = .data$ymin, ymax = .data$ymax),
                       fill = "white", color = "#27b56e", linewidth = 0.65, show.legend = FALSE)+
    ggplot2::geom_segment(data = linePlot_df, mapping = ggplot2::aes(x = .data$FirstAA, xend = .data$LastAA,
                                                   y = maxY + 1, yend = maxY + 1, size = .data$size,
                                                   ), color = linePlot_df$color, size = linePlot_df$size,
                          show.legend = FALSE) +
    ggrepel::geom_text_repel(data = text_df, ggplot2::aes(x = .data$centerOfAA, y = maxY + 1.2, label = .data$Domain),
                             direction = "x", show.legend = FALSE, color = text_df$color,
                             , max.overlaps = Inf, segment.color = NA, verbose = FALSE) +
    ggplot2::geom_point(data = labeldf, ggplot2::aes(x=.data$ProteinPTMLocalization, y = maxY + 1),
                        fill = "#27b56e", color = "white", size = 9, shape = 21, show.legend = FALSE) +
    ggrepel::geom_label_repel(data = labeldf, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                           y = maxY + 1, label = .data$ModificationID),
                              max.overlaps = Inf, nudge_y = 0.4, label.size = NA, fill = NA,
                              color = "#27b56e", show.legend = FALSE, verbose = FALSE) +
    ggplot2::geom_point(data = point_df, ggplot2::aes(x = .data$x_corrected, y = .data$y_corrected,
                                                      color = .data$GlycanType),
                        size = 3, show.legend = TRUE) +
    ggplot2::annotate("text", x = 1 - pointSpacing, y = maxY + 1, label = "1") +
    ggplot2::annotate("text", x = df_sum$ProteinLength[1] + 2 * pointSpacing,
                      y = maxY + 1, label = df_sum$ProteinLength[1]) +
    ggplot2::labs(title = plotTitle) +
    ggplot2::scale_y_continuous(limits =  c(-1, maxY +3)) +
    ggplot2::scale_color_manual(values = color_map, name = "Glycan Type") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    print(p)

}
