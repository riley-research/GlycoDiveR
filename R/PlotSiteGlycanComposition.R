PlotSiteGlycanComposition <- function(input, protein, whichAlias = NULL,
                                      horizontalPoints = 50, yCorrection = 0.5,
                                      yNudge = 1, boxSpacing = 0.05){

  set.seed(1)
  #input = mydata
  # protein = "P12259,A0A0A0MRJ7"
  # protein = "A0A0G2JH66,A0A0G2JPA3,A0A140T8X0,A0A140T9B7,A0A140T9V0,D0EV57"
  # protein = "Q07954"
  #protein = "P02675,D6REL8"
  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(!grepl("C\\(57.0215|M\\(15.9949", .data$AssignedModifications)) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco") %>%
    dplyr::filter(.data$UniprotIDs == protein)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){
    return(fmessage("No glycosylation found on the protein"))
  }

  df_sum <- df %>%
    dplyr::distinct(ModificationID, TotalGlycanComposition, .keep_all = TRUE) %>%
    dplyr::summarise(.by = c(ModificationID, GlycanType, Domains, ProteinLength, ProteinPTMLocalization),
                     count = dplyr::n())

  #Get the protein plot if there is domain information####
  linePlot_df <- df_sum[1,] %>%
    tidyr::separate_longer_delim(cols = "Domains", delim = ";") %>%
    dplyr::mutate(Domain = sub("\\(.*", "", Domains),
                  FirstAA = as.integer(sub(".*\\(([^-]+)-.*", "\\1", Domains)),
                  LastAA = as.integer(sub(".*-(.*)\\)", "\\1", Domains))) %>%
    dplyr::select(Domain, FirstAA, LastAA) %>%
    dplyr::bind_rows(data.frame(Domain = "Peptide", FirstAA = 1, LastAA = df_sum$ProteinLength[1])) %>%
    dplyr::mutate(size = ifelse(Domain == "Peptide", 5, 7),
                  centerOfAA = (FirstAA + LastAA) / 2) %>%
    dplyr::arrange(desc(Domain == "Peptide"))

  #Set colors####
  unique_domains <- unique(linePlot_df$Domain)
  colors <- setNames(c("#00394c", "#91b5c5", "#a9a9a9", "#22031F", "#433E0E" )[1:length(unique_domains)], unique_domains)
  linePlot_df <- linePlot_df %>%
    dplyr::mutate(color = colors[Domain])

  #Get the labeling right
  text_df <- linePlot_df[linePlot_df$Domain != "Peptide",] %>%
    dplyr::sample_frac(1) %>%
    dplyr::distinct(Domain, .keep_all=TRUE)

  #Now the dots####
  #For each site calculate the boxes
  point_df <- df %>%
    dplyr::distinct(.data$ProteinPTMLocalization,
                    .data$ModificationID,
                    .data$TotalGlycanComposition,
                    .data$GlycanType,
                    .data$ProteinPTMLocalization) %>%
    dplyr::summarise(.by = c(ModificationID,.data$GlycanType, .data$ProteinPTMLocalization),
                     count = dplyr::n()) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  colCount = ceiling(sqrt(sum(count))),
                  count = purrr::map_chr(count, ~ paste(1:.x, collapse = ";"))) %>%
    tidyr::separate_longer_delim(cols = "count", delim = ";") %>%
    dplyr::arrange(.data$GlycanType) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  count = as.numeric(count),
                  cumCount = 1:length(.data$count),
                  x = ((cumCount - 1) %% colCount) + 1,
                  y = ceiling(.data$cumCount / .data$colCount))

    #Make sure points are top centered
    maxPointY <- max(point_df$y, na.rm = TRUE)
    point_df <- point_df %>%
      dplyr::mutate(.by = .data$ModificationID,
                    y = (maxPointY - max(.data$y)) + y,
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
                  y_corrected = (y * yCorrection) + yNudge)

  #Check for any clashes in point####
  #Start with the second group Check for any clashes with the previous group
  idVec <- point_df %>%
    dplyr::arrange(ProteinPTMLocalization) %>%
    dplyr::distinct(ModificationID) %>%
    dplyr::pull(ModificationID)

  point_df$Corrected <- "No"

  if(length(idVec) > 1){
    for(i in 2:length(idVec)){
      maxGroup1 <- point_df %>%
        dplyr::filter(ModificationID == idVec[i - 1]) %>%
        dplyr::pull(x_corrected) %>%
        max()

      minGroup2 <- point_df %>%
        dplyr::filter(ModificationID == idVec[i]) %>%
        dplyr::pull(x_corrected) %>%
        min()

      if(maxGroup1 > minGroup2){
        correctionFactor <- (maxGroup1 - minGroup2) + boxSpacing * df_sum$ProteinLength[1]

        point_df <- point_df %>%
          dplyr::mutate(x_corrected = ifelse(ModificationID == idVec[i],
                                             x_corrected + correctionFactor,
                                             x_corrected),
                        Corrected = ifelse(ModificationID == idVec[i],
                                           "Yes", Corrected))
      }
    }
  }

  #Get the point colors based on the glycan composition####
  point_df <- point_df %>%
    dplyr::mutate(glycanColor = dplyr::case_when(.data$GlycanType == "Complex/Hybrid" ~ "#00394a",
                                                 .data$GlycanType == "Sialyl+Fucose" ~ "#ff7f2a",
                                                 .data$GlycanType == "Sialyl" ~ "#2475b5",
                                                 .data$GlycanType == "Fucose" ~ "#aaaaaa",
                                                 .data$GlycanType == "High Mannose" ~ "#28b36d"))

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
                  boxX = ifelse(Corrected == "No", ProteinPTMLocalization, mean(c(min(x_corrected), max(x_corrected)))),
                  pointY = maxY + 1) %>%
    dplyr::distinct() %>%
    dplyr::mutate(xCords = paste(.data$ProteinPTMLocalization, .data$boxX, sep = ";"),
                  yCords = paste(.data$pointY, .data$maximumY, sep = ";")) %>%
    dplyr::mutate(.by = .data$ModificationID,
                  xCords = ifelse(Corrected == "No", xCords,
                                  calculateElbowCoords(xVec = xCords, yVec = yCords, return = "x")),
                  yCords = ifelse(Corrected == "No", yCords,
                                  calculateElbowCoords(xVec = xCords, yVec = yCords, return = "y"))) %>%
    tidyr::separate_longer_delim(cols = c("xCords", "yCords"), delim = ";") %>%
    dplyr::mutate(xCords = as.double(.data$xCords),
                  yCords = as.double(.data$yCords))

  #colormapping####
  color_map <- setNames(point_df$glycanColor, point_df$GlycanType)

  #Now plot the actual plot####
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = line_df, ggplot2::aes(x = .data$xCords, y = .data$yCords, group = .data$ModificationID),
                       color = "#27b56e", linewidth = 1, show.legend = FALSE) +
    ggplot2::geom_rect(data = rect_df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = "white", color = "#27b56e", linewidth = 1, show.legend = FALSE)+
    ggplot2::geom_segment(data = linePlot_df, mapping = aes(x = FirstAA, xend = LastAA,
                                                   y = maxY + 1, yend = maxY + 1, size = .data$size,
                                                   ), color = linePlot_df$color, size = linePlot_df$size, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = text_df, ggplot2::aes(x = centerOfAA, y = maxY + 1.3, label = Domain),
                             direction = "x", show.legend = FALSE, color = text_df$color) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_point(data = labeldf, ggplot2::aes(x=ProteinPTMLocalization, y = maxY + 1),
                        fill = "#27b56e", color = "white", size = 9, shape = 21, show.legend = FALSE) +
    ggrepel::geom_label_repel(data = labeldf, ggplot2::aes(x =.data$ProteinPTMLocalization,
                                                           y = maxY + 1, label = .data$ModificationID),
                              max.overlaps = Inf, nudge_y = 0.4, label.size = NA, fill = NA,
                              color = "#27b56e", show.legend = FALSE) +
    ggplot2::geom_point(data = point_df, ggplot2::aes(x = .data$x_corrected, y = .data$y_corrected,
                                                      color = .data$GlycanType),
                        size = 3, show.legend = TRUE) +
    ggplot2::scale_y_continuous(limits =  c(-1, maxY +3)) +
    # ggplot2::scale_color_identity(guide = "legend") +
    ggplot2::scale_color_manual(values = color_map, name = "Glycan Type") +
    ggplot2::theme_void()

    print(p)

}
