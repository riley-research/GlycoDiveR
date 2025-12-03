#' ManageColorScheme
#'
#' Manage the color palette used for plotting
#'
#' @param edit Supply a vector containing R supported colors. These can for example
#' be "white" or "#E1D89F".
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
#' ManageColorScheme(edit = c("white", "black", "#E1D89F"), action="Append")
#' }
ManageColorScheme <- function(edit = NULL, action = c("View", "Append", "Replace", "Reset")) {
  action <- match.arg(action)

  if(action == "View") {
    print(.modEnv$colorScheme)
  } else if(action == "Append") {
    .modEnv$colorScheme <- c(.modEnv$colorScheme, edit)
  } else if(action == "Replace") {
    .modEnv$colorScheme <- edit
  } else if(action == "Reset") {
    if(exists("colorScheme")) {
      .modEnv$colorScheme <- colorScheme
    } else {
      stop("No original colorScheme found to reset.")
    }
  }

  invisible(.modEnv$colorScheme)
}
