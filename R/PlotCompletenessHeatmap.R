#' Plot Completeness Heatmap
#'
#' @param input Formatted data
#' @param grouping peptide or glyco
#' @param peptideType glyco or other
#' @param exportDataTo provide path to a folder to export the heatmap as a csv file
#'
#' @returns Grouped heatmap
#' @export
#'
#' @examples \dontrun{PlotCompletenessHeatmap(mydata)}
PlotCompletenessHeatmap <- function(input, grouping = "peptide", peptideType  = "glyco", exportDataTo = FALSE){
  input <- FilterForCutoffs(input)

  if(grouping == "peptide"){
    df <- GetMeanTechReps(input$PSMTable)
  }

  df$GlycanType <- sapply(df$GlycanType, function(x) gsub(", |NonGlyco|Unmodified", "", x))

  if(peptideType == "glyco"){
    df <- subset(df, GlycanType != "")
  }

  df <- df[,c("Alias", "GlycanType", "ModifiedPeptide", "Run")] %>%
    tidyr::pivot_wider(names_from = Alias, values_from = Run, values_fn = dplyr::first) %>%
    dplyr::mutate(dplyr::across(-c(ModifiedPeptide, GlycanType), ~ ifelse(!is.na(.), 1, 0)))

  #Get the matrix and in the right column order
  mtrx <- data.matrix(df[,3:ncol(df)])
  rownames(mtrx) <- df$ModifiedPeptide

  mtrx <- mtrx[,levels(input$PTMTable$Alias), drop = FALSE]

  if(!identical(exportDataTo, FALSE)){
    write.csv(mtrx, paste0(exportDataTo, "/CompletenessHeatmapData.csv"))
  }

  #Get heatmap annotation, color, and legend
  colH <- setNames(colorScheme[1:length(unique(df$GlycanType))], unique(df$GlycanType))

  row_ha = ComplexHeatmap::rowAnnotation(Glycan = df$GlycanType, show_legend = FALSE,
                                         col = list(Glycan = colH))

  col_fun = circlize::colorRamp2(c(0,1), c("lightgrey", "darkgreen"))
  col_fun(seq(-3, 3))

  lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "", at = c(0, 1),
               labels = c("Missing", "Identified"))

  #Now for clustering the columns
  colCluster <- sapply(colnames(mtrx), function(x) subset(input$annotation, Alias == x)$Condition[1])

  ht <- ComplexHeatmap::Heatmap(mtrx, show_row_names = FALSE, cluster_rows = FALSE,
                               cluster_columns = FALSE, left_annotation = row_ha,
                               row_split = df$GlycanType, use_raster = FALSE,
                               cluster_row_slices = TRUE, border = TRUE,
                               col = col_fun, show_heatmap_legend = FALSE,
                               column_split = colCluster)

  print(ComplexHeatmap::draw(ht, annotation_legend_list = list(lgd)))
}
