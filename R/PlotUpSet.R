#' PlotUpSet
#'
#' Plot UpSet plots to find identification overlaps.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param grouping Choose between "condition", "biologicalReps", and
#' "technicalReps".
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param type Choose between "glyco" or "all" to include only glycopeptides/proteins
#' or all peptides/proteins.
#' @param level Choose "peptide" or "protein".
#' @param plotColor The color of the bars.
#' @param nintersects The maximum number of intersects shown.
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
#' @param silent TRUE if you want info to be printed, FALSE if not.
#'
#' @returns An UpSet plot
#' @export
#'
#' @examples
#' \dontrun{
#'  PlotUpSet(mydata, grouping = "biologicalReps",
#'  type = "glyco", level = "protein", nintersects = 40)
#' }
PlotUpSet <- function(input, grouping = "condition", type = "glyco",
                      level = "peptide", plotColor = "#32006e", whichAlias = NULL,
                      whichProtein = NULL, exactProteinMatch = TRUE,
                      nintersects = 40, whichPeptide = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(type == "glyco"){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "")
  }else if(type == "all"){
  }else{
    warning("Your type argument is not recognized. Including all peptides.")
  }

  if(nrow(input$PSMTable) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if(grouping == "technicalReps"){
    if(level == "peptide"){
      df_list <- split(input$PSMTable$ModifiedPeptide, input$PSMTable$Alias)
    }else if(level == "protein"){
      df_list <- split(input$PSMTable$UniprotIDs, input$PSMTable$Alias)
    }else{
      stop("Check your level argument. Your input is not recognized.")
    }

    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                  nsets = length(names(df_list)), nintersects = nintersects)

    #Extract labels and get color coding right
    colorCode <- as.data.frame(trimws(p$labels)) %>%
      dplyr::left_join(
        input$PSMTable[, c("Alias", "Condition")] %>% dplyr::distinct(.data$Alias, .data$Condition),
        by = c("trimws(p$labels)" = "Alias")
      ) %>%
      dplyr::mutate(
       Condition = factor(.data$Condition, levels = unique(.data$Condition)),
        colorScheme = .modEnv$colorScheme[as.integer(.data$Condition)]
      )

    #replot with correct colors
    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)),
                       sets.bar.color = colorCode$colorScheme,
                       nintersects = nintersects, main.bar.color = plotColor,
                       matrix.color = plotColor, text.scale = 1.5)

    return(p)
  }else if(grouping == "biologicalReps"){
    if(level == "peptide"){
      df <- input$PSMTable %>%
        dplyr::mutate(ID = paste0(.data$Condition, .data$BioReplicate)) %>%
        dplyr::summarise(.by = c("ID", "ModifiedPeptide", "Condition"))

      df_list <- split(df$ModifiedPeptide, df$ID)
    }else if(level == "protein"){
      df <- input$PSMTable %>%
        dplyr::mutate(ID = paste0(.data$Condition, .data$BioReplicate)) %>%
        dplyr::summarise(.by = c("ID", "UniprotIDs", "Condition"))

      df_list <- split(df$UniprotIDs, df$ID)
    }else{
      stop("Check your level argument. Your input is not recognized.")
    }

    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)),
                       nintersects = nintersects)

    #Extract labels and get color coding right
    colorCode <- as.data.frame(trimws(p$labels)) %>%
      dplyr::left_join(
        df[, c("ID", "Condition")] %>% dplyr::distinct(.data$ID, .data$Condition),
        by = c("trimws(p$labels)" = "ID")
      ) %>%
      dplyr::mutate(
        Condition = factor(.data$Condition, levels = unique(.data$Condition)),
        colorScheme = .modEnv$colorScheme[as.integer(.data$Condition)]
      )

    #replot with correct colors
    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)),
                       sets.bar.color = colorCode$colorScheme,
                       nintersects = nintersects, main.bar.color = plotColor,
                       matrix.color = plotColor, text.scale = 1.5)

    return(p)
  }else if(grouping == "condition"){
    if(level == "peptide"){
      df <- input$PSMTable %>%
        dplyr::summarise(.by = c("ModifiedPeptide", "Condition"))

      df_list <- split(df$ModifiedPeptide, df$Condition)
    }else if(level == "protein"){
      df <- input$PSMTable %>%
        dplyr::summarise(.by = c("UniprotIDs", "Condition"))

      df_list <- split(df$UniprotIDs, df$Condition)
    }else{
      stop("Check your level argument. Your input is not recognized.")
    }

    p <- suppressMessages(UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), nintersects = nintersects))

    #Extract labels and get color coding right
    colorCode <- data.frame(Condition = trimws(p$labels)) %>%
      dplyr::mutate(
        Condition = factor(.data$Condition, levels = unique(.data$Condition)),
        colorScheme = .modEnv$colorScheme[as.integer(.data$Condition)]
      )

    #replot with correct colors
    p <- suppressMessages(UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)),
                       sets.bar.color = colorCode$colorScheme,
                       nintersects = nintersects, main.bar.color = plotColor,
                       matrix.color = plotColor, text.scale = 1.5))

    return(p)
  }else{
    stop("Check your grouping argument. Your input is not recognized.")
  }
  return(p)
}
