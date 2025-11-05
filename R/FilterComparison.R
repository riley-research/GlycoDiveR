#' FilterComparison
#'
#' This functions helps filter your comparison dataframe, which can subsequently
#' be used for other plots and functions. The dataframe needs to have a column
#' supplied by the whichStatistic argument (typically, "pvalue", or "adjpvalue").
#' Alternatively, you can filter the comparison dataframe yourself, for example,
#' using subset or dplyr::filter.
#'
#' @param comparison A comparison dataframe.
#' @param statisticalCutoff What statistical cutoff to use (default = 0.05)
#' @param whichStatistic Typically, "pvalue", or "adjpvalue". Must correspond to
#' an identically named column in your comparison dataframe.
#' @param whichLabel What label(s) to use. Corresponds to the Label column of
#' the input comparison dataframe.
#' @param minf TRUE/FALSE if -inf should be retained (default = FALSE).
#' @param inf TRUE/FALSE if inf should be retained (default = FALSE)
#'
#' @returns A subset of your comparison dataframe
#' @export
#'
#' @examples \dontrun{FilterComparison(comparison, whichLabel =
#' c("SampleA-SampleB", "SampleB-SampleC")) }
FilterComparison <- function(comparison, statisticalCutoff = 0.05, whichStatistic = "pvalue",
                             whichLabel = NA, minf = FALSE, inf = FALSE){

  if(is.data.frame(comparison) && whichStatistic %in% names(comparison)){
    if("log2FC" %in% names(comparison) && !minf){
      comparison <- comparison %>%
        dplyr::filter(.data$log2FC != -Inf)
    }
    if("log2FC" %in% names(comparison) && !inf){
      comparison <- comparison %>%
        dplyr::filter(.data$log2FC != Inf)
    }
    if(!is.na(whichLabel)){
      comparison <- comparison %>%
        dplyr::filter(.data$Label %in% whichLabel)
    }
    comparison <- comparison %>%
      dplyr::filter(!!rlang::sym(whichStatistic) <= statisticalCutoff)
  }else{
    stop("Input is either not a dataframe, or does not have the columns ", whichStatistic)
  }
  return(comparison)
}
