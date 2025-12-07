#' Plot Completeness Matrix
#'
#' Plot a matrix showing the identified and missed (glyco)peptides. Each row
#' represents one unique glycopeptide.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param peptideType Choose "glyco" to only view glycopeptide and "other" to  also
#' include non-modified peptides.
#' @param plotColors The colors used for the matrix The default is c("lightgrey", "darkgreen").
#' @param collapseTechReps Do you want to collapse the technical replicates.
#' @param exportDataTo If a system folder is provided, a CSV file will be exported
#' to that folder containing the data presented in the matrix
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting, e.g. whichAlias = c("Alias1", "Alias2").
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
#' @returns A matrix.
#' @export
#'
#' @examples \dontrun{
#' PlotCompletenessMatrix(mydata)
#' PlotCompletenessMatrix(mydata, peptideType = "other", silent = FALSE)
#' }
PlotCompletenessMatrix <- function(input, peptideType  = "glyco",
                                    plotColors = c("grey90", "#44AA99"),
                                    collapseTechReps = FALSE,
                                    whichAlias = NULL, whichPeptide = NULL,
                                    whichProtein = NULL, exactProteinMatch = TRUE,
                                    exportDataTo = FALSE, silent = FALSE){

  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- input$PSMTable %>%
    dplyr::filter(!is.na(.data$Intensity))

  if(collapseTechReps){
    df <- GetMeanTechReps(input$PSMTable)
  }else{
    df <- input$PSMTable
  }

  df <- df %>%
    dplyr::mutate(
      GlycanType = gsub(", NonGlyco|NonGlyco, | NonGlyco|, Unmodified| Unmodified, |Unmodified", "", .data$GlycanType),
      GlycanType = ifelse(grepl(", ", .data$GlycanType), "Multi", .data$GlycanType),
      GlycanType = ifelse(.data$GlycanType == "", "NonGlyco", .data$GlycanType))

  if(peptideType == "glyco"){
    df <- df %>%
      dplyr::filter(.data$GlycanType != "" & .data$GlycanType != "NonGlyco")
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

  df <- df[,c("Alias", "GlycanType", "ModifiedPeptide", "Run")] %>%
    tidyr::pivot_wider(names_from = "Alias", values_from = "Run", values_fn = dplyr::first) %>%
    dplyr::mutate(dplyr::across(-c("ModifiedPeptide", "GlycanType"), ~ ifelse(!is.na(.), 1, 0)))

  #Get the matrix and in the right column order
  mtrx <- data.matrix(df[,3:ncol(df)])
  rownames(mtrx) <- df$ModifiedPeptide

  cols_to_use <- levels(input$PTMTable$Alias)
  cols_to_use <- cols_to_use[cols_to_use %in% colnames(mtrx)]
  mtrx <- mtrx[, cols_to_use, drop = FALSE]

  if(!identical(exportDataTo, FALSE)){
    utils::write.csv(mtrx, paste0(exportDataTo, "/CompletenessHeatmapData.csv"))
  }

  #Get heatmap annotation, color, and legend
  colH <- stats::setNames(c(.modEnv$GlycanColors$color, "#BBBBBB"),
                          c(.modEnv$GlycanColors$GlycanType, "NonGlyco"))

  row_ha = ComplexHeatmap::rowAnnotation(Glycan = df$GlycanType, show_legend = FALSE,
                                         col = list(Glycan = colH))

  col_fun <- c("0" = plotColors[1],
               "1" = plotColors[2])

  lgd <- ComplexHeatmap::Legend(at = c("0", "1"),
    labels = c("Missing", "Identified"),
    legend_gp = grid::gpar(fill = col_fun[c("0", "1")]),
    title = "")

  #Now for clustering the columns
  colCluster <- sapply(colnames(mtrx), function(x) input$annotation$Condition[input$annotation$Alias == x][1])

  ht <- ComplexHeatmap::Heatmap(mtrx, show_row_names = FALSE, cluster_rows = FALSE,
                               cluster_columns = FALSE, left_annotation = row_ha,
                               row_split = df$GlycanType, use_raster = FALSE,
                               cluster_row_slices = TRUE, border = TRUE,
                               col = col_fun, show_heatmap_legend = FALSE,
                               column_split = colCluster, row_title_rot = 0)

  return(ComplexHeatmap::draw(ht, annotation_legend_list = list(lgd)))
}
