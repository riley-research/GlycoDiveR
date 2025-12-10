#' PlotSiteGlycanCompositionBar
#'
#' Visualize the glycan composition per site.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping Grouping is "none", "technicalReps", "biologicalReps", or "condition".
#' @param intensity Raw log2, or log10 transformed intensity values.
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
#' }
PlotSiteGlycanCompositionBar <- function(input, grouping = "condition", scales = "fill",
                                         intensity = "raw", whichAlias = NULL,
                                         whichPeptide = NULL, whichProtein = NULL,
                                         exactProteinMatch = TRUE, silent = FALSE){
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

  if(CheckForQuantitativeValues(df$Intensity)){
    if(!silent){
      return(fmessage("No quantitative data found."))
    }else{
      return()
    }
  }

  if(nrow(df) == 0){return("Nothing left after filtering.")}

  if(scales == "fill"){
    yAxis <- "Summed intensity\n(Normalized to total)"
  }else{
    yAxis <- "Summed intensity"
  }

  if(intensity == "log2"){
    df$Intensity <- log(df$Intensity, 2)
  }else if(intensity == "log10"){
    df$Intensity <- log(df$Intensity, 10)
  }

  colH <- stats::setNames(c(.modEnv$GlycanColors$color, "#BBBBBB"),
                            c(.modEnv$GlycanColors$GlycanType, "None"))

  if(grouping == "none"){
    p <- df %>%
      dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0) %>%
      dplyr::arrange(.data$ProteinPTMLocalization) %>%
      dplyr::summarise(.by = c("UniprotIDs", "Genes", "ModificationID", "GlycanType"),
                       Intensity = sum(.data$Intensity, na.rm=TRUE)) %>%
      tidyr::complete(.data$UniprotIDs, .data$Genes,
                      .data$ModificationID, fill = list(Intensity = 0.001,
                                                        GlycanType = "None")) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$ModificationID, y=.data$Intensity, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
  }else if(grouping == "technicalReps"){
    p <- df %>%
      dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0) %>%
      dplyr::arrange(.data$ProteinPTMLocalization) %>%
      dplyr::summarise(.by = c("Alias","UniprotIDs", "Genes", "ModificationID", "GlycanType"),
                       Intensity = sum(.data$Intensity, na.rm=TRUE)) %>%
      tidyr::complete(.data$UniprotIDs, .data$Genes,
                      .data$Alias, .data$ModificationID, fill = list(Intensity = 0.001,
                                                                     GlycanType = "None")) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$ModificationID, y=.data$Intensity, fill=.data$GlycanType)) +
      ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
      ggplot2::facet_wrap(~.data$Alias) +
      ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face="bold", size = 12))
    }else if(grouping == "biologicalReps"){
      p <- df %>%
        dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0) %>%
        dplyr::arrange(.data$ProteinPTMLocalization) %>%
        dplyr::mutate(x = paste0(.data$Condition, .data$BioReplicate)) %>%
        dplyr::summarise(.by = c("x","UniprotIDs", "Genes", "ModificationID", "GlycanType"),
                         Intensity = sum(.data$Intensity, na.rm=TRUE)) %>%
        tidyr::complete(.data$UniprotIDs, .data$Genes,
                        .data$x, .data$ModificationID, fill = list(Intensity = 0.001,
                                                                   GlycanType = "None")) %>%
        ggplot2::ggplot(ggplot2::aes(x=.data$ModificationID, y=.data$Intensity, fill=.data$GlycanType)) +
        ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
        ggplot2::facet_wrap(~.data$x) +
        ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
        ggplot2::scale_y_continuous(expand=c(0,0)) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(face="bold", size = 12))
      }else if(grouping == "condition"){
        p <- df %>%
          dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0) %>%
          dplyr::arrange(.data$ProteinPTMLocalization) %>%
          dplyr::summarise(.by = c("Condition","UniprotIDs", "Genes", "ModificationID", "GlycanType"),
                           Intensity = sum(.data$Intensity, na.rm=TRUE)) %>%
          tidyr::complete(.data$UniprotIDs, .data$Genes,
                          .data$Condition, .data$ModificationID, fill = list(Intensity = 0.001,
                                                                             GlycanType = "None")) %>%
          ggplot2::ggplot(ggplot2::aes(x=.data$ModificationID, y=.data$Intensity, fill=.data$GlycanType)) +
          ggplot2::geom_bar(stat="identity", position = scales, width=1, color= "black") +
          ggplot2::facet_wrap(~.data$Condition) +
          ggplot2::scale_fill_manual(values = colH, guide = ggplot2::guide_legend(reverse = TRUE)) +
          ggplot2::labs(fill = NULL, y = yAxis, x = NULL) +
          ggplot2::scale_y_continuous(expand=c(0,0)) +
          ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                         strip.background = ggplot2::element_blank(),
                         strip.text = ggplot2::element_text(face="bold", size = 12))
  }else{
    return(fmessage("Your grouping argument was not recognized."))
  }

  return(p)
}
