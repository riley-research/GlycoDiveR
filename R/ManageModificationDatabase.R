#' ManageModificationDatabase
#'
#' Manage the modification database.
#'
#' @param edit Supply a path the a csv file containing the headers "FullName" and
#' "ModificationMass", or supply a data frame with the headers "FullName" and
#' "ModificationMass".
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
#' ManageModificationDatabase(edit = mydataframe, action="Append")
#' }
ManageModificationDatabase <- function(edit = NULL, action = c("View", "Append", "Replace", "Reset")) {
  action <- match.arg(action)

  prepare_edit <- function(edit) {
    if(!is.data.frame(edit)) {
      edit <- utils::read.csv(edit)
    }
    edit %>%
      dplyr::select("FullName", "ModificationMass") %>%
      dplyr::mutate(
        FullName = as.character(.data$FullName),
        ModificationMass = as.character(.data$ModificationMass)
      )
  }

  if(action == "View") {
    print(.modEnv$ModificationDatabase)
  } else if(action == "Append") {
    edit <- prepare_edit(edit)
    .modEnv$ModificationDatabase <- dplyr::bind_rows(.modEnv$ModificationDatabase, edit)
  } else if(action == "Replace") {
    edit <- prepare_edit(edit)
    .modEnv$ModificationDatabase <- edit
  } else if(action == "Reset") {
    if(exists("ModificationDatabase")) {
      .modEnv$ModificationDatabase <- ModificationDatabase
    } else {
      stop("No original ModificationDatabase found to reset.")
    }
  }

  invisible(.modEnv$ModificationDatabase)
}
