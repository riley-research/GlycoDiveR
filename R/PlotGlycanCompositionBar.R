#' PlotGlycanCompositionBar
#'
#' Visualize the glycan composition per protein, multiple proteins, or all the
#' proteins.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping Grouping is "technicalReps", "biologicalReps", or "condition".
#' @param scales  Controls plot normalization, choose "fill" or "stack".
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
#' @returns A barplot.
#' @export
#'
#' @examples \dontrun{
#' PlotGlycanCompositionBar(mydatam scales = "stack")
#'
#'
#' }
PlotGlycanCompositionBar <- function(input, grouping = "condition", scales = "fill",
                                     whichAlias = NULL, whichPeptide = NULL,
                                     whichProtein = NULL, exactProteinMatch = TRUE,
                                     silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein,exactProteinMatch)
  input$PTMTable <- input$PTMTable %>% dplyr::filter(!is.na(.data$Intensity))

  if(grouping != "technicalReps"){
    df <- GetMeanTechReps(input$PTMTable) %>%
      dplyr::filter(.data$GlycanType != "NonGlyco")
  }else{
    df <- input$PTMTable %>%
      dplyr::filter(.data$GlycanType != "NonGlyco")
  }

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){return("Nothing left after filtering.")}

  if(scales == "fill"){
    yAxis <- "Number of glycosites\n(Normalized to total)"
  }else{
    yAxis <- "Number of glycosites"
  }

  colH <- stats::setNames(.modEnv$GlycanColors$color, .modEnv$GlycanColors$GlycanType)

  if(grouping == "technicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Run", "Alias", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Run, .data$Alias, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Alias, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                     axis.ticks.y = ggplot2::element_blank())

    return(p)
  }else if(grouping == "biologicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "BioReplicate",
                               "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$BioReplicate, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$x, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                     axis.ticks.y = ggplot2::element_blank())

    return(p)
  }else if(grouping == "condition"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Condition, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                     axis.ticks.y = ggplot2::element_blank())

    return(p)
  }else{
    fmessage("Your grouping argument was not recognized.")
    return(NULL)
  }
}
