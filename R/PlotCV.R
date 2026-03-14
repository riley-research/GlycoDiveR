#' PlotCV
#'
#' Plot the protein CVs using violing plots. The number of proteins in each group
#' are labeled. Choose between glyco only, or all, and decide if you want to return
#' a plot or the underlying data. The returned data can directly be used to filter
#' other plots using the whichPeptide argument.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param horizontalLines What horizontal lines to show in the plot.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param type Choose "glyco" or "all"
#' @param returnType Choose "plot" to return the plot, "data" to return the
#' CV data as a dataframe (this dataframe can be used to filter any visualization
#' using the whichPeptide argument), or "GDdata" to return filtered GlycoDiveR
#' data.
#' @param returnThreshold Works in concert with returnType = "data". Defines the
#' the upper CV threshold.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param silent silence printed information
#'
#' @returns A ggplot, dataframe, or GlycoDiveR data
#' @export
#'
#' @examples \dontrun{
#' PlotCV(mydata, type = "glyco", returnType = "plot")
#'
#' PlotCV(mydata, type = "all", returnType = "GDdata", returnThreshold = 20)
#' }
PlotCV <- function(input, horizontalLines = c(0, 20, 50, 100), whichProtein = NULL,
                   type = c("all", "glyco"), returnType = c("plot", "data", "GDdata"),
                   returnThreshold = NULL, whichPeptide = NULL, exactProteinMatch = TRUE,
                   whichAlias = NULL, silent = FALSE){
  type <- match.arg(type)
  returnType <- match.arg(returnType)
  inputRaw <- input
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  df <- FilterForPeptides(input$PSMTable, whichPeptide)
  df <- GetMeanTechReps(df)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  #keep only glyco and remove 0 values, and keep one peptide per Alias
  if (type == "glyco") {
    df <- df %>%
      dplyr::mutate(PSMType = GetPSMGlycanCategory(.data$GlycanType)) %>%
      dplyr::filter(.data$PSMType != "nonGlyco")
  }

  df <- df %>%
    dplyr::filter(.data$Intensity != 0 & !is.na(.data$Intensity)) %>%
    dplyr::arrange(dplyr::desc(.data$Intensity)) %>%
    dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .keep_all = TRUE)

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if(CheckForQuantitativeValues(df$Intensity)){
    if(!silent){
      return(fmessage("No quantitative data found."))
    }else{
      return()
    }
  }

  df_col <- df %>%
    dplyr::mutate(group = .data$Condition) %>%
    dplyr::arrange(.data$Alias) %>%
    dplyr::mutate(group = factor(.data$group, levels = unique(.data$group))) %>%
    dplyr::summarise(.by = "Condition", group = list(.data$group)) %>%
    dplyr::mutate(col = .modEnv$colorScheme[1:dplyr::n()]) %>%
    tidyr::unnest(cols = "group") %>%
    dplyr::distinct(.data$group, .data$col)

  colH <- stats::setNames(df_col$col, df_col$group)

  df <- df %>%
    dplyr::summarise(.by = c("Condition", "ModifiedPeptide"),
                     CV = stats::sd(.data$Intensity, na.rm = TRUE) /
                       mean(.data$Intensity, na.rm = TRUE) * 100) %>%
    dplyr::filter(!is.na(.data$CV))

  if(returnType == "plot"){
    p <- df_label <- df %>%
      dplyr::summarise(.by = "Condition", numberOfPep = dplyr::n())

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Condition, y = .data$CV, fill = .data$Condition)) +
      ggplot2::geom_hline(yintercept = horizontalLines, linetype = "dashed",
                          color = "grey75") +
      ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA) +
      ggplot2::geom_text(data = df_label, ggplot2::aes(x = .data$Condition,
                                                       y = max(df$CV, na.rm = TRUE) * 1.10,
                                                       label = .data$numberOfPep)) +
      ggplot2::labs(x = NULL, y = "CV (%)") +
      ggplot2::scale_fill_manual(values = colH) +
      ggplot2::theme(legend.position = "none")

    return(p)
  }else if (returnType == "data"){
    if(is.null(returnThreshold)) {
      return(df)
    }else if (is.numeric(returnThreshold) &&
              length(returnThreshold) == 1 &&
              is.finite(returnThreshold)){
      df <- df %>%
        dplyr::filter(.data$CV <= returnThreshold)

      return(df)
    }else {
      warning("Did not understand returnThreshold")
    }
  }else if (returnType == "GDdata") {
    inputRaw$PSMTable <- inputRaw$PSMTable %>%
      dplyr::filter(.data$ModifiedPeptide %in% df$ModifiedPeptide)

    inputRaw$PTMTable <- inputRaw$PTMTable %>%
      dplyr::filter(.data$ModifiedPeptide %in% df$ModifiedPeptide)

    return(inputRaw)
  } else {
      warning("Did not understand returnType")
    }
}
