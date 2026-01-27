#' PlotPCA
#'
#' PCA plot based on PSM quantification. Missing values are imputed using
#' kNN.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param quantType "normalized" for normalized intensity values, and "nonNormalized"
#' for non-normalized intensity values.
#' @param dropNoQuant Remove or keep PSMs without a quantitative value or have a quant
#' of 0.
#' @param minPeptideCoverage Works together with the thresholdMode argument.
#' If thresholdMode = "total", provide a single numeric value. Values between 0
#' and 1 are interpreted as proportions (e.g., 0.1 = 10%), while values >= 1 are
#' interpreted as absolute counts (e.g., 2 = at least two peptides identified).
#' If thresholdMode = "group", provide a numeric vector of length two (e.g., c(0.1, 3)).
#' The first value is interpreted in the same way as for "total" (proportion or
#' absolute count). The second value specifies the minimum number of groups in
#' which this threshold must be met.
#' @param thresholdMode Character string specifying the filtering mode: "total"
#' applies a global threshold across all data, while "group" identifies the threshold
#' within groups.
#' @param grouping Grouping is "technicalReps", "biologicalReps", or "condition".
#' @param type Either "glyco" or "all".
#' @param label TRUE for showing Alias labels, FALSE for no labels.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
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
#' @returns A ggplot PCA plot.
#' @export
#'
#' @examples \dontrun{
#'      PlotPCA(mydata)
#'
#'      PlotPCA(mydata, quantType = "normalized", minPeptideCoverage = c(0.4, 2),
#'      thresholdMode = "group", grouping = "biologicalReps", type = "all")
#'
#' }
PlotPCA <- function(input, quantType = "normalized", dropNoQuant = TRUE,
                    minPeptideCoverage = 0.5, thresholdMode = "total",
                    grouping = "technicalReps", type = "glyco", label =TRUE,
                    whichProtein = NULL, whichPeptide = NULL, exactProteinMatch = TRUE,
                    whichAlias = NULL, silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)
  df <- FilterForPeptides(input$PSMTable, whichPeptide)

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(type == "glyco"){
    df <- df %>%
      dplyr::mutate(PSMType = GetPSMGlycanCategory(.data$GlycanType)) %>%
      dplyr::filter(.data$PSMType != "nonGlyco")
  }

  if(grouping != "technicalReps"){
    df <- GetMeanTechReps(df)
  }

  if(dropNoQuant){df <- df %>% dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0)}
  if(!identical(minPeptideCoverage, FALSE)){df <- FilterForMinPeptides(df, minPeptideCoverage, thresholdMode)}

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

  if(quantType == "normalized"){
    qt <- "Intensity"
  }else if(quantType == "nonNormalized"){
    qt <- "RawIntensity"
  }else{
    stop("Did not recognize your quantType. Choose 'normalized', or 'nonNormalized'.")
  }

  if(grouping == "technicalReps"){
    df <- df %>%
      dplyr::select("Alias", "ModifiedPeptide", "Quantification" = dplyr::all_of(qt)) %>%
      dplyr::mutate(Quantification = log(.data$Quantification,2),
                    Quantification = ifelse(is.finite(.data$Quantification),
                                            .data$Quantification, NA_integer_)) %>%
      dplyr::arrange(dplyr::desc(.data$Quantification)) %>%
      dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      tidyr::pivot_wider(names_from = "Alias", values_from = "Quantification")

    mtrx <- as.matrix(df[,2:ncol(df)])
    rownames(mtrx) <- df$ModifiedPeptide

    mtrx <- impute::impute.knn(mtrx, rowmax = 1)$data

    pca <- stats::prcomp(t(mtrx), center = TRUE, scale. = TRUE)
    pca_scores <- as.data.frame(pca$x)

    pca_scores$Alias <- rownames(pca_scores)

    pca_scores <- pca_scores %>%
      dplyr::left_join(input$PSMTable %>% dplyr::select("Alias", "Condition") %>% dplyr::distinct(),
                       by = "Alias")

    if(identical(label, FALSE)){pca_scores$Alias <- NA}

    colH <- stats::setNames(
      .modEnv$colorScheme[1:length(unique(pca_scores$Condition))],
      unique(pca_scores$Condition)
    )

    p <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$Condition)) +
      ggplot2::geom_point(size = 5) +
      ggrepel::geom_label_repel(ggplot2::aes(label = .data$Alias), fill = NA, label.size = NA) +
      ggplot2::labs(
        x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
        y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
        color = NULL
      ) +
      ggplot2::scale_color_manual(values = colH) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else if(grouping == "biologicalReps"){
    df <- df %>%
      dplyr::select("Condition", "BioReplicate", "ModifiedPeptide", "Quantification" = dplyr::all_of(qt)) %>%
      dplyr::mutate(Quantification = log(.data$Quantification,2),
                    Quantification = ifelse(is.finite(.data$Quantification),
                                            .data$Quantification, NA_integer_)) %>%
      dplyr::arrange(dplyr::desc(.data$Quantification)) %>%
      dplyr::distinct(.data$Condition, .data$BioReplicate, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::mutate(Alias = paste(.data$Condition, .data$BioReplicate, sep = "")) %>%
      dplyr::select(-c("Condition", "BioReplicate")) %>%
      tidyr::pivot_wider(names_from = "Alias", values_from = "Quantification")

    mtrx <- as.matrix(df[,2:ncol(df)])
    rownames(mtrx) <- df$ModifiedPeptide

    mtrx <- impute::impute.knn(mtrx, rowmax = 1)$data

    pca <- stats::prcomp(t(mtrx), center = TRUE, scale. = TRUE)
    pca_scores <- as.data.frame(pca$x)

    pca_scores$Alias <- rownames(pca_scores)

    pca_scores <- pca_scores %>%
      dplyr::left_join(input$PSMTable %>% dplyr::select("Alias", "Condition") %>% dplyr::distinct(),
                       by = "Alias")

    if(identical(label, FALSE)){pca_scores$Alias <- NA}

    colH <- stats::setNames(
      .modEnv$colorScheme[1:length(unique(pca_scores$Condition))],
      unique(pca_scores$Condition)
    )

    p <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$Condition)) +
      ggplot2::geom_point(size = 5) +
      ggrepel::geom_label_repel(ggplot2::aes(label = .data$Alias), fill = NA, label.size = NA) +
      ggplot2::labs(
        x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
        y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
        color = NULL
      ) +
      ggplot2::scale_color_manual(values = colH) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else if(grouping == "condition"){
    df <- df %>%
      dplyr::select("Alias", "Condition", "ModifiedPeptide", "Quantification" = dplyr::all_of(qt)) %>%
      dplyr::mutate(Quantification = log(.data$Quantification,2),
                    Quantification = ifelse(is.finite(.data$Quantification),
                                            .data$Quantification, NA_integer_)) %>%
      dplyr::arrange(dplyr::desc(.data$Quantification)) %>%
      dplyr::distinct(.data$Alias, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::summarise(.by = c("Condition", "ModifiedPeptide"),
                       Quantification = stats::median(.data$Quantification, na.rm = TRUE)) %>%
      dplyr::select("Alias" = "Condition", "ModifiedPeptide", "Quantification") %>%
      tidyr::pivot_wider(names_from = "Alias", values_from = "Quantification")

    mtrx <- as.matrix(df[,2:ncol(df)])
    rownames(mtrx) <- df$ModifiedPeptide

    mtrx <- impute::impute.knn(mtrx, rowmax = 1)$data

    pca <- stats::prcomp(t(mtrx), center = TRUE, scale. = TRUE)
    pca_scores <- as.data.frame(pca$x)

    pca_scores$Alias <- rownames(pca_scores)

    if(identical(label, FALSE)){pca_scores$Alias <- NA}

    colH <- stats::setNames(
      .modEnv$colorScheme[1:length(unique(pca_scores$Alias))],
      unique(pca_scores$Alias)
    )

    p <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$Alias)) +
      ggplot2::geom_point(size = 5) +
      ggrepel::geom_label_repel(ggplot2::aes(label = .data$Alias), fill = NA, label.size = NA) +
      ggplot2::labs(
        x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
        y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
        color = NULL
      ) +
      ggplot2::scale_color_manual(values = colH) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    return(p)
  }else{
    stop("Your grouping argument was not recognized")
  }
}
