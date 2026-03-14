#' PlotPeptideLength
#'
#' Visualize the identified peptide lengths using a histogram.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param type Choose "both", "glyco", or "nonGlyco".
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
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
#' @param silent silence printed information
#'
#' @returns A ggplot graph.
#' @export
#'
#' @examples \dontrun{
#' PlotPeptideLength(mydata, input = "glyco")
#' }
PlotPeptideLength <- function(input, type = "both", whichAlias = NULL,
                              whichPeptide = NULL, whichProtein = NULL,
                              exactProteinMatch = TRUE, silent = FALSE) {
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  input$PSMTable <- input$PSMTable %>% dplyr::filter(!is.na(.data$Intensity))

  input$PSMTable$Glycan <- sapply(input$PSMTable$TotalGlycanComposition, function(x) ifelse(!is.na(x) & x != "", "Glycosylated", "nonGlycosylated"))

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

  if(type == "both"){
    input$PSMTable %>%
      dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .data$Glycan) %>%
      dplyr::mutate(nakedPep = gsub("[^A-Z]", "", .data$ModifiedPeptide),
                    pepLength = nchar(.data$nakedPep, allowNA = TRUE, keepNA = TRUE),
                    Glycan = factor(.data$Glycan, levels = c("nonGlycosylated", "Glycosylated"))) %>%
      dplyr::summarise(.by = c("Glycan", "pepLength"),
                       count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(.data$count)) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pepLength,
                                   y = .data$count,
                                   fill = .data$Glycan)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::labs(x = "Peptide length", y = "Peptide count",
                    fill = NULL) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5, vjust = 0.5, angle = 0)) +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme)
  }else if (type == "glyco") {
    input$PSMTable %>%
      dplyr::filter(.data$Glycan == "Glycosylated") %>%
      dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
      dplyr::mutate(nakedPep = gsub("[^A-Z]", "", .data$ModifiedPeptide),
                    pepLength = nchar(.data$nakedPep, allowNA = TRUE, keepNA = TRUE)) %>%
      dplyr::summarise(.by ="pepLength",
                       count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(.data$count)) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pepLength,
                                   y = .data$count)) +
      ggplot2::geom_bar(stat = "identity", position = "identity",
                        color = "black", fill = .modEnv$colorScheme[1]) +
      ggplot2::labs(x = "Peptide length", y = "Peptide count") +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5, vjust = 0.5, angle = 0))
  }else if (type == "nonGlyco") {
    input$PSMTable %>%
      dplyr::filter(.data$Glycan == "nonGlycosylated") %>%
      dplyr::distinct(.data$Alias, .data$ModifiedPeptide) %>%
      dplyr::mutate(nakedPep = gsub("[^A-Z]", "", .data$ModifiedPeptide),
                    pepLength = nchar(.data$nakedPep, allowNA = TRUE, keepNA = TRUE)) %>%
      dplyr::summarise(.by ="pepLength",
                       count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(.data$count)) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pepLength,
                                   y = .data$count)) +
      ggplot2::geom_bar(stat = "identity", position = "identity",
                        color = "black", fill = .modEnv$colorScheme[2]) +
      ggplot2::labs(x = "Peptide length", y = "Peptide count") +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5, vjust = 0.5, angle = 0))
  }else {
    stop("Did not understand your type argument")
  }
}

