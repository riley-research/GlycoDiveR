#' Manage glycan categories
#'
#' Manage the glycan categories used by GlycoDiveR. Run `ManageGlycanCategories(action = "View")`
#' to see what this database looks like. To use custom glycan categories, please run this function before
#' you import your data. The colors for new glycan categories should be registered using `ManageGlycanColors()`.
#'
#' @param edit A data frame used to replace the existing glycan categories.
#'   You can generate a template for this data frame using `action = "showCode"`.
#'   The columns for the dataframe are:
#'
#'   \describe{
#'     \item{Priority}{Numeric priority used to resolve matches. The category with the highest priority is assigned.}
#'     \item{Residue}{Amino acid at the glycosylation site.}
#'     \item{Pattern}{Pattern matched against the `TotalGlycanComposition` column (see `mydata$PSMTable$TotalGlycanComposition`).
#'                    The pattern is a Regular Expressions (Regex). AI is super helpful getting the Regex
#'                    that you need.}
#'     \item{Min_H / Max_H}{Minimum and maximum number of hexoses.}
#'     \item{Min_N / Max_N}{Minimum and maximum number of HexNAc residues.}
#'     \item{Category}{Assigned glycan category for the highest-priority match.}
#'   }
#'
#'   If you introduce new glycan categories, use `ManageGlycanColors()` to assign
#'   colors to them.
#'
#' @param action Character string specifying the action to perform. One of:
#'   \describe{
#'     \item{"View"}{Return the current glycan categories.}
#'     \item{"Replace"}{Replace the current glycan categories with `edit`.}
#'     \item{"Reset"}{Reset glycan categories to the default values.}
#'     \item{"showCode"}{Return example code for creating a new `edit` data frame.}
#'   }
#'
#' @return A data frame containing the (updated) glycan categories.
#' @export
#'
#' @examples
#' \dontrun{
#' # View current categories
#' ManageGlycanCategories(action = "View")
#'
#' # Replace categories using a template
#' ManageGlycanCategories(
#'   edit = df_from_action_equals_showCode,
#'   action = "Replace"
#' )
#' }
ManageGlycanCategories <- function(edit = NULL, action = c("View", "Replace", "Reset", "showCode")) {
  action <- match.arg(action)

  if(action == "View") {
    print(.modEnv$GlycanCategories)
  } else if(action == "Replace") {
    .modEnv$GlycanCategories <- edit
    warning(
      paste0("Important: If adding new categories, ensure you use ",
        "ManageGlycanColors() to register colors for them.")
    )
  } else if(action == "Reset") {
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
  } else if(action == "showCode"){
    code <- 'edit <- tibble::tribble(
      ~Priority, ~Residue, ~Pattern,             ~Min_H, ~Max_H, ~Min_N, ~Max_N, ~Category,
      #----------|---------|---------------------|-------|-------|-------|-------|------------------
      1,         "N",      "(?=.*(A|G))(?=.*F)", 0,      Inf,    0,      Inf,    "Sialofucosylated",
      2,         "N",      "A|G",                0,      Inf,    0,      Inf,    "Sialylated",
      3,         "N",      "F",                  0,      Inf,    0,      Inf,    "Fucosylated",
      4,         "N",      "Phospho",            0,      Inf,    0,      Inf,    "Phosphomannose",
      5,         "N",      ".*",                 0,      3,      0,      2,      "Truncated",
      6,         "N",      ".*",                 4,      Inf,    0,      2,      "Oligomannose",
      7,         "N",      ".*",                 0,      Inf,    0,      Inf,    "Complex/Hybrid",

      # O-linked logic
      8,         "S|T",    "(?=.*(A|G))(?=.*F)", 0,      Inf,    0,      Inf,    "Sialofucosylated",
      9,         "S|T",    "A|G",                0,      Inf,    0,      Inf,    "Sialylated",
      10,        "S|T",    "F",                  0,      Inf,    0,      Inf,    "Fucosylated",
      11,        "S|T",    ".*",                 0,      Inf,    0,      Inf,    "OGlycan"
    )'

    return(cat(code, "\n"))
  }

  invisible(.modEnv$GlycanDatabase)
}
