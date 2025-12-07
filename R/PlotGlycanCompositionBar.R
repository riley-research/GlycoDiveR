PlotGlycanCompositionBar <- function(input, grouping, position = "fill",
                                     whichAlias = NULL, whichPeptide = NA,
                                     whichProtein = NA, exactProteinMatch = TRUE,
                                     silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein,exactProteinMatch)

  df <- GetMeanTechReps(input$PTMTable) %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){return("Nothing left after filtering.")}

  if(grouping == "technicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Run", "Alias", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Run, .data$Alias, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Condition, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = position, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else if(grouping == "biologicalReps"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "BioReplicate",
                               "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$BioReplicate, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Condition, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = position, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else if(grouping == "condition"){
    p <- df %>%
      dplyr::summarise(.by = c("Condition", "GlycanType"),
                       GlycanCount = dplyr::n()) %>%
      tidyr::complete(.data$Condition, .data$GlycanType, fill = list(GlycanCount = 0)) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Condition, y=.data$GlycanCount, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = position, width=1, color= "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colorScheme, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = NULL, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else{
    fmessage("Your grouping argument was not recognized.")
    return(NULL)
  }
}
