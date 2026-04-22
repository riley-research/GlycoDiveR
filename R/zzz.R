#' @importFrom rlang .data
NULL
.onLoad <- function(libname, pkgname) {
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                               panel.grid = ggplot2::element_blank()))

  .modEnv$GlycanDatabase <- GlycanDatabase
  .modEnv$ModificationDatabase <- ModificationDatabase
  .modEnv$colorScheme <- colorScheme
  .modEnv$GlycanColors <- GlycanColors
  .modEnv$useExtendedOGlycanCategories <- FALSE
  .modEnv$GlycanCategories <- data.frame(
    Priority = 1:11,
    Residue  = c("N", "N", "N", "N", "N", "N", "N", "S|T", "S|T", "S|T", "S|T"),
    Pattern  = c(
      "(?=.*(A|G))(?=.*F)", "A|G", "F", "Phospho",
      ".*", ".*", ".*", # N-linked placeholders
      "(?=.*(A|G))(?=.*F)", "A|G", "F", ".*" # O-linked placeholders
    ),
    Min_H    = c(0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0),
    Max_H    = c(Inf, Inf, Inf, Inf, 3, Inf, Inf, Inf, Inf, Inf, Inf),
    Min_N    = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    Max_N    = c(Inf, Inf, Inf, Inf, 2, 2, Inf, Inf, Inf, Inf, Inf),
    Category = c(
      "Sialofucosylated",
      "Sialylated",
      "Fucosylated",
      "Phosphomannose",
      "Truncated",      # Priority 5: N <= 2 AND H <= 3
      "Oligomannose",   # Priority 6: N <= 2 AND H >= 4
      "Complex/Hybrid", # Priority 7: N > 2
      "Sialofucosylated",
      "Sialylated",
      "Fucosylated",
      "OGlycan"
    ),
    stringsAsFactors = FALSE
  )
}
