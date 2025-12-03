#' ManageGlycanDatabase
#'
#' Manage the glycan database.
#'
#' @param edit Supply a path the a csv file containing the headers "FullName",
#' "ShortName" and "GMass", or supply a data frame with those headers.
#' @param action Choose between the following:
#' View: to view the current modification database.
#' Append: append to the currrent modification database with the supplied edit.
#' Replace: replace the modification database with the supplied edit.
#' Reset: reset to the default values.
#'
#' @returns The updated database.
#' @export
#'
#' @examples \dontrun{
#' ManageGlycanDatabase(edit = "C:/locationOfDataframe.csv", action="Append")
#' }
ManageGlycanDatabase <- function(edit = NULL, action = c("View", "Append", "Replace", "Reset")) {
  action <- match.arg(action)

  prepare_edit <- function(edit) {
    if(!is.data.frame(edit)) {
      edit <- utils::read.csv(edit)
    }
    edit %>%
      dplyr::select("FullName", "ShortName", "GMass") %>%
      dplyr::mutate(
        FullName = as.character(.data$FullName),
        ShortName = as.character(.data$ShortName),
        GMass = as.numeric(.data$GMass)
      )
  }

  if(action == "View") {
    print(.modEnv$GlycanDatabase)
  } else if(action == "Append") {
    edit <- prepare_edit(edit)

    if(any(edit$FullName %in% .modEnv$GlycanDatabase$FullName) ||
       any(edit$ShortName %in% .modEnv$GlycanDatabase$ShortName)){
      stop("Duplicated found in either FullName or ShortName column. Please check.")
    }

    .modEnv$GlycanDatabase <- dplyr::bind_rows(.modEnv$GlycanDatabase, edit)
  } else if(action == "Replace") {
    edit <- prepare_edit(edit)
    .modEnv$GlycanDatabase <- edit
  } else if(action == "Reset") {
    if(exists("GlycanDatabase")) {
      .modEnv$GlycanDatabase <- GlycanDatabase
    } else {
      stop("No original GlycanDatabase found to reset.")
    }
  }

  invisible(.modEnv$GlycanDatabase)
}
