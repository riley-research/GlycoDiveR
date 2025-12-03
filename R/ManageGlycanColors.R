#' ManageGlycanColors
#'
#' Manage the color palette used for glycans.
#'
#' @param edit Supply a path the a csv file containing the headers "GlycanType",
#' and "color", or supply a data frame with those headers.
#' @param action Choose between the following:
#' View: to view the current modification database.
#' Append: append to the currrent modification database with the supplied edit. If
#' duplicate GlycanType values are found, it will overwrite the existing values.
#' Replace: replace the modification database with the supplied edit.
#' Reset: reset to the default values.
#'
#' @returns The updated database.
#' @export
#'
#' @examples \dontrun{
#' ManageGlycanColors(edit = mydataframe, action="Append")
#' }
ManageGlycanColors <- function(edit = NULL, action = c("View", "Append", "Replace", "Reset")) {
  action <- match.arg(action)

  prepare_edit <- function(edit) {
    if(!is.data.frame(edit)) {
      edit <- utils::read.csv(edit)
    }
    edit %>%
      dplyr::select("GlycanType", "color") %>%
      dplyr::mutate(
        GlycanType = as.character(.data$GlycanType),
        color = as.character(.data$color)
      )
  }

  if(action == "View") {
    print(.modEnv$GlycanColors)
  } else if(action == "Append") {
    edit <- prepare_edit(edit)

    .modEnv$GlycanColors <- .modEnv$GlycanColors %>%
      dplyr::select(!(.data$GlycanType %in% edit$GlycanType))

    .modEnv$GlycanColors <- dplyr::bind_rows(.modEnv$GlycanColors, edit)
  } else if(action == "Replace") {
    edit <- prepare_edit(edit)
    .modEnv$GlycanColors <- edit
  } else if(action == "Reset") {
    if(exists("GlycanColors")) {
      .modEnv$GlycanColors <- GlycanColors
    } else {
      stop("No original GlycanColors found to reset.")
    }
  }

  invisible(.modEnv$GlycanColors)
}
