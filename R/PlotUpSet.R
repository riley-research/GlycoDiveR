#' Plot UpSet plots to find identification overlaps
#'
#' @param input Formatted data
#' @param grouping choose between "condition", "biologicalReps", and
#' "technicalReps"
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param type choose between "glyco" or "all" to include only glycopeptides/proteins
#' or all peptides/proteins
#' @param level choose "peptide" or "protein"
#' @param nintersects The maximum number of intersects shown
#'
#' @returns An UpSet plot
#' @export
#'
#' @examples
#' \dontrun{
#'  PlotUpSet(mydata, grouping = "biologicalReps",
#'  type = "glyco", level = "protein", nintersects = 40)
#' }
PlotUpSet <- function(input, grouping = "condition", whichAlias = NULL,
                      type = "glyco", level = "peptide", nintersects = 40){
  input <- FilterForCutoffs(input)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(Alias %in% whichAlias)
  }

  if(type == "glyco"){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(!is.na(TotalGlycanComposition) & TotalGlycanComposition != "")
  }else if(type == "all"){
  }else{
    warning("Your type argument is not recognized. Including all peptides.")
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
                  nsets = length(names(df_list)), nintersects = 40)

    #Extract labels and get color coding right
    colorCode <- as.data.frame(trimws(p$labels)) %>%
      dplyr::left_join(
        input$PSMTable[, c("Alias", "Condition")] %>% dplyr::distinct(Alias, Condition),
        by = c("trimws(p$labels)" = "Alias")
      ) %>%
      dplyr::mutate(
       Condition = factor(Condition, levels = unique(Condition)),
        colorScheme = colorScheme[as.integer(Condition)]
      )

    #replot with correct colors
    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), cutoff = cutoff,
                       sets.bar.color = colorCode$colorScheme)

    print(p)
  }else if(grouping == "biologicalReps"){
    if(level == "peptide"){
      df <- input$PSMTable %>%
        dplyr::mutate(ID = paste0(Condition, BioReplicate)) %>%
        dplyr::summarise(.by = c(ID, ModifiedPeptide, Condition))

      df_list <- split(df$ModifiedPeptide, df$ID)
    }else if(level == "protein"){
      df <- input$PSMTable %>%
        dplyr::mutate(ID = paste0(Condition, BioReplicate)) %>%
        dplyr::summarise(.by = c(ID, UniprotIDs, Condition))

      df_list <- split(df$UniprotIDs, df$ID)
    }else{
      stop("Check your level argument. Your input is not recognized.")
    }

    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), nintersects = 40)

    #Extract labels and get color coding right
    colorCode <- as.data.frame(trimws(p$labels)) %>%
      dplyr::left_join(
        df[, c("ID", "Condition")] %>% dplyr::distinct(ID, Condition),
        by = c("trimws(p$labels)" = "ID")
      ) %>%
      dplyr::mutate(
        Condition = factor(Condition, levels = unique(Condition)),
        colorScheme = colorScheme[as.integer(Condition)]
      )

    #replot with correct colors
    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), cutoff = cutoff,
                       sets.bar.color = colorCode$colorScheme)

    print(p)
  }else if(grouping == "condition"){
    if(level == "peptide"){
      df <- input$PSMTable %>%
        dplyr::summarise(.by = c(ModifiedPeptide, Condition))

      df_list <- split(df$ModifiedPeptide, df$Condition)
    }else if(level == "protein"){
      df <- input$PSMTable %>%
        dplyr::summarise(.by = c(UniprotIDs, Condition))

      df_list <- split(df$UniprotIDs, df$Condition)
    }else{
      stop("Check your level argument. Your input is not recognized.")
    }

    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), nintersects = 40)

    #Extract labels and get color coding right
    colorCode <- data.frame(Condition = trimws(p$labels)) %>%
      dplyr::mutate(
        Condition = factor(Condition, levels = unique(Condition)),
        colorScheme = colorScheme[as.integer(Condition)]
      )

    #replot with correct colors
    p <- UpSetR::upset(UpSetR::fromList(df_list), order.by = "freq",
                       nsets = length(names(df_list)), cutoff = cutoff,
                       sets.bar.color = colorCode$colorScheme)

    print(p)
  }else{
    stop("Check your grouping argument. Your input is not recognized.")
  }
  return(p)
}
