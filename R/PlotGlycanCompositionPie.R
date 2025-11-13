#' PlotGlycanCompositionPie
#'
#' @param input Formatted data
#' @param grouping Grouping is "technicalReps", "biologicalReps", or "condition"
#' @param scales Controls plot normalization, choose "fixed" or "free",
#' @param protein Use "all" for all proteins, otherwise provide a vector containing
#' proteins from your data  as in the UniprotIDs column
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#'
#' @returns The glycan composition
#' @export
#'
#' @examples \dontrun{PlotGlycanCompositionPie(mydata, grouping = "technicalReps",
#' protein = c("P10204", "Q92930"))}
PlotGlycanCompositionPie <- function(input, grouping, scales = "free",
                                     protein = "all", whichAlias = NULL,
                                     whichPeptide = NA){
  input <- FilterForCutoffs(input)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

  df <- GetMeanTechReps(input$PTMTable) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(protein != "all"){
    df <- df %>%
      dplyr::filter(.data$UniprotIDs %in% protein)

    if (nrow(df) == 0) {
      return("No glycans found after filtering.")
    }
  }

  if(grouping == "technicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c(.data$Run, .data$Alias, .data$GlycanType),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Run, .data$Alias, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, color= "black") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$Alias, scales = scales) +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())

    return(p)
  }else if(grouping == "biologicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c(.data$Condition, .data$BioReplicate,
                               .data$GlycanType),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$BioReplicate, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, color= "black") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$x, scales = scales) +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(p)
  }else if(grouping == "condition"){
    p <- df %>%
      dplyr::summarise(.by = c(.data$Condition, .data$GlycanType),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x="", y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", width=1, color= "black") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::facet_wrap(~.data$Condition, scales = scales) +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(p)
  }else{
    fmessage("Your grouping argument was not recognized.")
    return(NULL)
  }
}
