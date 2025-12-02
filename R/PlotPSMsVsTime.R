#' PlotPSMsVsTime
#'
#' Visualize the number of (glyco)PSMs versus the retention time using a lineplot.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param type Choose between "allGlyco", "both", "glyco", "combined", "all", or
#' supply a vector of glycan types, such as c("Multi", "nonGlyco", "Sialyl", "Complex/Hybrid",
#' "Sialyl+Fucose", "Fucose", "Truncated", "High Mannose", "Paucimannose").
#' @param binWidth the bin width used.
#' @param whichAlias  Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param gradientLength Define the LC gradient length. If no length is supplied
#' the maximum PSM retention time is used.
#' @param plotColors Define the colors of the non-glyco types. Default = c("#bdbdbd", "#1b9e77").
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param silent silence printed information (default = FALSE)
#'
#' @returns a lineplot
#' @export
#'
#' @examples \dontrun{
#' PlotPSMsVsTime(mydata)
#'
#' PlotPSMsVsTime(mydata, type = "allGlyco", whichAlias = c("NP0-1-r1"),
#' binWidth = 1, gradientLength = NA)
#' }
PlotPSMsVsTime <- function(input, type = "all", binWidth = 5, gradientLength = NA,
                           plotColors = c("#bdbdbd", "#1b9e77"), whichAlias = NULL,
                           whichPeptide = NULL, whichProtein = NULL,
                           exactProteinMatch = TRUE, silent = FALSE){
  glycoPSMTypes <- .modEnv$GlycanColors$GlycanType

  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(input$PSMTable) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  df <- input$PSMTable %>%
    dplyr::mutate(PSMType = dplyr::case_when(stringr::str_count(.data$GlycanType, paste(glycoPSMTypes, collapse = "|")) > 1 &
                                               stringr::str_count(.data$GlycanType, "Sialyl") == 2 &
                                               grepl("Sialyl+Fucose", .data$GlycanType) == 1 ~ "Multi",
                                             stringr::str_count(.data$GlycanType, paste(glycoPSMTypes, collapse = "|")) == 0 ~ "nonGlyco",
                                             grepl("Sialyl", .data$GlycanType) ~ "Sialyl",
                                             grepl("Complex/Hybrid", .data$GlycanType) ~ "Complex/Hybrid",
                                             grepl("Sialyl+Fucose", .data$GlycanType) ~ "Sialyl+Fucose",
                                             grepl("Fucose", .data$GlycanType) ~ "Fucose",
                                             grepl("Truncated", .data$GlycanType) ~ "Truncated",
                                             grepl("High Mannose", .data$GlycanType) ~ "High Mannose",
                                             grepl("Paucimannose", .data$GlycanType) ~ "Paucimannose",
                                             grepl("OGlycan", .data$GlycanType) ~ "OGlycan",
                                             grepl("NonCanonicalGlyco", .data$GlycanType) ~ "NonCanonicalGlyco",
                                             TRUE ~ "Other"))

  if (length(type) == 1 && identical(type, "allGlyco")) {
    df <- df %>%
      dplyr::filter(.data$PSMType != "nonGlyco")

  } else if (length(type) == 1 && identical(type, "both")) {
    df <- df %>%
      dplyr::mutate(PSMType = ifelse(.data$PSMType == "nonGlyco", "nonGlyco", "Glyco"))

  } else if(length(type) == 1 && identical(type, "glyco")){
    df <- df %>%
      dplyr::mutate(PSMType = ifelse(.data$PSMType == "nonGlyco", "nonGlyco", "glycoPSMs")) %>%
      dplyr::filter(.data$PSMType == "glycoPSMs")

  }else if(length(type) == 1 && identical(type, "combined")){
    df$PSMType <- "PSMs"

  }else if (!(length(type) == 1 && identical(type, "all"))) {
    df <- df %>%
      dplyr::filter(.data$PSMType %in% type)

    if(nrow(df) == 0){
      stop("Tried filtering for: ", type, ", but no values remain.")
    }
  }

  #Now the binning####
  if(!is.na(gradientLength)){
    bins <- seq(0, ceiling(gradientLength / binWidth) * binWidth, binWidth)
  }else{
    bins <- seq(0, ceiling(max(df$RetentionTime) / binWidth) * binWidth, binWidth)
  }

  df <- df %>%
    dplyr::mutate(rawbin = cut(.data$RetentionTime, breaks = bins, include.lowest = TRUE),
                  lowerBound = stringr::str_extract(.data$rawbin, "(?<=\\()[^,]+") |> as.numeric(),
                  upperBound = stringr::str_extract(.data$rawbin, "(?<=,)[^\\]]+") |> as.numeric())

  df$bin <- (df$lowerBound + df$upperBound) / 2

  #Get the summed dataframe####
  df <- df %>%
    dplyr::summarise(.by = c("Alias", "PSMType", "bin"),
              count = dplyr::n())

  #Fill in with 0s for all bins with no IDs
  allLevels <- c()
  for(i in 1:length(bins)){
    if(i == 1){
      next
    }else{
      if(identical(allLevels, c())){
        allLevels <- mean(c(bins[i], bins[i-1]))
      }else{
        allLevels <- c(allLevels, mean(c(bins[i], bins[i-1])))
      }
    }
  }
  df$bin <- factor(df$bin, levels = allLevels)
  df$Alias <- droplevels(df$Alias)
  df <- df %>%
    tidyr::complete(.data$Alias, .data$PSMType, .data$bin, fill = list(count = 0))
  df$bin <- as.numeric(as.character(droplevels(df$bin)))

  #Get the colors right####
  colH <- stats::setNames(c(.modEnv$GlycanColors$color, plotColors[1], rep(plotColors[2], 3)),
                          c(.modEnv$GlycanColors$GlycanType, "nonGlyco", "glycoPSMs", "Glyco", "PSMs"))

  #Plot###
  ggplot2::ggplot(df, ggplot2::aes(x = .data$bin, y = .data$count,
                                   linetype = .data$Alias, color = .data$PSMType)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Time (min)", y = "Number of PSMs") +
    ggplot2::scale_color_manual(values = colH)

}
