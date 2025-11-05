#' GetComparison
#'
#' Computes pvalues, BH-corrected pvalues, and log2 fold changes for the specified
#' conditions. Comparisons are done on the glycopeptide level using the
#' max intensity value for duplicate PSMs.
#'
#' @param input formatted data
#' @param comparisons the comparisons you would like to include
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
#' @param type Choose "glyco" to use only glycopeptides. Otherwise
#' choose "all".
#'
#' @returns a comparison dataframe
#' @export
#'
#' @examples \dontrun{GetComparison(mydata, list(c("PBS", "Lis")))}
GetComparison <- function(input, comparisons, type = "glyco", whichAlias = NULL){
  df <- GetMeanTechReps(input$PSMTable)

  if(type == "glyco"){
    df <- df %>%
      dplyr::filter(.data$GlycanType != "NonGlyco" & .data$GlycanType != "Unmodified")
  }

  df <- df %>%
    dplyr::select(.data$UniprotIDs, "Proteins" = .data$Genes,
                  .data$ModifiedPeptide, .data$Condition,
                  .data$Intensity, .data$BioReplicate) %>%
    dplyr::mutate(sampleName = paste(.data$Condition,
                                     .data$BioReplicate, sep = ";"),
                  Intensity = ifelse(.data$Intensity == 0, NA, .data$Intensity),
                  Intensity = log(.data$Intensity, 2)) %>%
    dplyr::select(-c(.data$Condition, .data$BioReplicate)) %>%
    tidyr::pivot_wider(names_from = .data$sampleName, values_from = "Intensity")

  tempColumns <- names(df)[4:length(names(df))]

  for(i in 1:length(comparisons)){
    compi <- comparisons[i][[1]]
    df[paste(compi[1], compi[2], sep="-")] <- NA
    index1 <- which(grepl(paste0(compi[1], ";"), names(df)))
    index2 <- which(grepl(paste0(compi[2], ";"), names(df)))

    if(length(index1) < 2 || length(index2) <2){
      warning("Tried to find this comparison: ", paste(compi, collapse = "-"),
              ", but at least one group does not have more than 1 measurements")
    }

    df[paste(compi[1], compi[2], sep="-")] <- apply(df, 1, function(row) {
      TTest_log2FC(as.numeric(row[index1]), as.numeric(row[index2]))
    })
  }

  df <- df %>%
    dplyr::select(-tidyselect::all_of(tempColumns)) %>%
    tidyr::pivot_longer(cols = -tidyselect::all_of(c("UniprotIDs", "Proteins", "ModifiedPeptide")),
                        names_to = "Label", values_to = "LabelLog2Fc") %>%
    tidyr::separate_wider_delim(cols = "LabelLog2Fc",  names = c("log2FC", "pvalue"), delim = ";") %>%
    dplyr::mutate(.by = "Label",
                  pvalue = ifelse(.data$pvalue == "NA", NA, .data$pvalue),
                  pvalue = as.double(.data$pvalue),
                  adjpvalue = stats::p.adjust(.data$pvalue, method = "BH"))
  return(df)
}
