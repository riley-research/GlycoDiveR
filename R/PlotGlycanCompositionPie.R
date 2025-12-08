#' PlotGlycanCompositionPie
#'
#' Visualize glycan compositions as a pie chart.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping Grouping is "technicalReps", "biologicalReps", or "condition".
#' @param scales Controls plot normalization, choose "fixed" or "free".
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param  exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param silent silence printed information (default = FALSE)
#'
#' @returns The glycan composition
#' @export
#'
#' @examples \dontrun{
#' PlotGlycanCompositionPie(mydata)
#'
#' PlotGlycanCompositionPie(mydata, grouping = "technicalReps",
#'                          whichProtein = c("P10204", "Q92930"))
#' }
PlotGlycanCompositionPie <- function(input, grouping = "condition", scales = "free",
                                     whichAlias = NULL, whichPeptide = NULL,
                                     whichProtein = NULL, exactProteinMatch = TRUE,
                                     silent = FALSE){

  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

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

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  colH <- stats::setNames(.modEnv$GlycanColors$color, .modEnv$GlycanColors$GlycanType)

  if(grouping == "technicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Run", "Alias", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Run, .data$Alias, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, linewidth=1.25, color= "white") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$Alias, scales = scales) +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())

    return(p)
  }else if(grouping == "biologicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "BioReplicate",
                               "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$BioReplicate, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, linewidth=1.25, color= "white") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$x, scales = scales) +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(p)
  }else if(grouping == "condition"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, linewidth=1.25, color= "white") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$Condition, scales = scales) +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(p)
  }else{
    fmessage("Your grouping argument was not recognized.")
    return(NULL)
  }
}
