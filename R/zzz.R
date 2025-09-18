#' @importFrom rlang .data
NULL
.onLoad <- function(libname, pkgname) {
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                               panel.grid = ggplot2::element_blank()))
}
