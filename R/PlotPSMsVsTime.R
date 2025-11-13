#' PlotPSMsVsTime
#'
#' @param input your formatted data
#' @param type Choose between "allGlyco", "both", "glyco", "combined", "all", or
#' supply a vector of glycan types, such as c("Multi", "nonGlyco", "Sialyl", "Complex/Hybrid",
#' "Sialyl+Fucose", "Fucose", "Truncated", "High Mannose", "Paucimannose")
#' @param bindwidth the bin width used
#' @param whichAlias  provide a vector of Aliases to only select these aliases
#' for plotting
#' @param gradientLength Define the LC gradient length. If no length is supplied
#' the maximum PSM retention time is used
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#'
#' @returns a lineplot
#' @export
#'
#' @examples \dontrun{PlotPSMsVsTime(mydata, type = "allGlyco", whichAlias = c("NP0-1-r1"),
#' bindwidth = 1, gradientLength = NA)
#' }
PlotPSMsVsTime <- function(input, type = "all", bindwidth = 5, whichAlias = NULL,
                           gradientLength = NA, whichPeptide = NA){
  glycoPSMTypes <- c("Sialyl", "Complex/Hybrid", "Sialyl+Fucose",
                "Fucose", "Truncated", "High Mannose", "Paucimannose")

  input <- FilterForCutoffs(input)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
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
                                             TRUE ~ "ERROR"))

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
    bins <- seq(0, ceiling(gradientLength / bindwidth) * bindwidth, bindwidth)
  }else{
    bins <- seq(0, ceiling(max(df$RetentionTime) / bindwidth) * bindwidth, bindwidth)
  }

  df <- df %>%
    dplyr::mutate(rawbin = cut(.data$RetentionTime, breaks = bins, include.lowest = TRUE),
                  lowerBound = stringr::str_extract(.data$rawbin, "(?<=\\()[^,]+") |> as.numeric(),
                  upperBound = stringr::str_extract(.data$rawbin, "(?<=,)[^\\]]+") |> as.numeric())

  df$bin <- (df$lowerBound + df$upperBound) / 2

  #Get the summed dataframe####
  df <- df %>%
    dplyr::summarise(.by = c(.data$Alias, .data$PSMType, .data$bin),
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
  glycan_colors <- c(
    "Complex/Hybrid" = "#1b9e77",
    "Sialyl+Fucose"  = "#d95f02",
    "Sialyl"         = "#7570b3",
    "Fucose"         = "#e7298a",
    "High Mannose"   = "#66a61e",
    "Truncated"      = "#a6761d",
    "Paucimannose"   = "#666666",
    "nonGlyco"       = "#bdbdbd",
    "glycoPSMs" = "#1b9e77",
    "Glyco" = "#1b9e77",
    "PSMs" = "#1b9e77"
  )

  #Plot###
  ggplot2::ggplot(df, ggplot2::aes(x = .data$bin, y = .data$count,
                                   linetype = .data$Alias, color = .data$PSMType)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Time (min)", y = "Number of PSMs") +
    ggplot2::scale_color_manual(values = glycan_colors)

}
