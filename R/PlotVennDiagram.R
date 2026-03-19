#' PlotVennDiagram
#'
#' Creates a Venn diagram based for findings overlap in peptides, proteins, or
#' glycans, either for glyco only or all. Only supports comparison between two
#' or three groups. Two-group comparisons will be scaled proportionally, while
#' three-way comparisons are not as this is mathemetically only possible in a few
#' cases.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param groups Supply a vector of two or three groups of choice. This works in
#' tandem with the grouping argument. If "technicalReps" is selected, grouping
#' corresponds to the values in the Alias column. For "biologicalReps" this is
#' the value in the Condition column and the BioReplicate column connected with a
#' _. For example, Condition_1. For "condition" this is the value in the Condition
#' column.
#' @param grouping Choose "condition", "biologicalReps", or "technicalReps".
#' @param type Choose either "glyco" to filter for only glycosylated PSMs, or
#' select "all".
#' @param level Compare either on the levels of "peptide", "protein", or
#' "glycan".
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param silent silence printed information (default = FALSE).
#'
#' @returns A ggplot visualization.
#' @export
#'
#' @examples \dontrun{
#' PlotVennDiagram(mydata, groups = c("Treated", "Untreated"), level = "glycan")
#'
#' PlotVennDiagram(mydata, groups = c("Treated_1", "Treated_2"),
#' grouping = "biologicalReps", level = "protein")
#' }
PlotVennDiagram <- function(input, groups, grouping = "condition", type = "glyco",
                      level = c("peptide", "protein", "glycan"),
                      whichAlias = NULL, whichProtein = NULL, exactProteinMatch = TRUE,
                      whichPeptide = NULL, silent = FALSE){
  level <- match.arg(level)
  input <- FilterForCutoffs(input, silent)
  input$PSMTable <- FilterForPeptides(input$PSMTable, whichPeptide)
  input$PSMTable <- FilterForProteins(input$PSMTable, whichProtein, exactProteinMatch)

  if(!is.null(whichAlias)){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(type == "glyco"){
    input$PSMTable <- input$PSMTable %>%
      dplyr::filter(!is.na(.data$TotalGlycanComposition) & .data$TotalGlycanComposition != "")
  }

  if(nrow(input$PSMTable) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  if (grouping == "technicalReps" & level != "glycan") {
    df <- input$PSMTable %>%
      dplyr::mutate(id = .data$Alias)
  }else if (grouping == "biologicalReps" & level != "glycan") {
    df <- GetMeanTechReps(input$PSMTable) %>%
      dplyr::mutate(id = paste0(.data$Condition, "-", .data$BioReplicate))
  }else if (grouping == "condition" & level != "glycan") {
    df <- GetMeanTechReps(input$PSMTable) %>%
      dplyr::mutate(id = .data$Condition)
  }else if (grouping == "technicalReps" & level == "glycan") {
    df <- input$PTMTable %>%
      tidyr::separate_longer_delim(cols = "TotalGlycanComposition", delim = ",") %>%
      dplyr::mutate(id = .data$Alias)
  }else if (grouping == "biologicalReps" & level == "glycan") {
    df <- GetMeanTechReps(input$PTMTable) %>%
      tidyr::separate_longer_delim(cols = "TotalGlycanComposition", delim = ",") %>%
      dplyr::mutate(id = paste0(.data$Condition, "-", .data$BioReplicate))
  }else if (grouping == "condition" & level == "glycan") {
    df <- GetMeanTechReps(input$PTMTable) %>%
      tidyr::separate_longer_delim(cols = "TotalGlycanComposition", delim = ",") %>%
      dplyr::mutate(id = .data$Condition)
  }else {
    warning("No match.")
  }

  if (level == "glycan") {
    df <- df %>%
      dplyr::filter(.data$id %in% groups) %>%
      dplyr::distinct(.data$id, .data$TotalGlycanComposition, .keep_all = TRUE) %>%
      dplyr::rename("var" = "TotalGlycanComposition")
  } else if (level == "peptide") {
    df <- df %>%
      dplyr::filter(.data$id %in% groups) %>%
      dplyr::distinct(.data$id, .data$ModifiedPeptide, .keep_all = TRUE) %>%
      dplyr::rename("var" = "ModifiedPeptide")
  } else {
    df <- df %>%
      dplyr::filter(.data$id %in% groups) %>%
      tidyr::separate_longer_delim(cols = "UniprotIDs", delim = ",") %>%
      dplyr::distinct(.data$id, .data$UniprotIDs, .keep_all = TRUE) %>%
      dplyr::rename("var" = "UniprotIDs")
  }

  group_result <- unique(as.character(df$id))

  if (length(group_result) != 2 & length(group_result) != 3) {
    if(!silent) {
      return(fmessage(paste0("Need 2 or 3 groups, you supplied ", length(group_result))))
    }else {
      return()
    }
  }

  if (length(group_result) == 3) {
    #Calculate overlap and difference
    group_list <- split(df$var, as.character(df$id))

    A_set <- group_list[[group_result[1]]]
    B_set <- group_list[[group_result[2]]]
    C_set <- group_list[[group_result[3]]]

    ABC_set <- Reduce(intersect, list(A_set, B_set, C_set))
    vals_excl <- c(
      A = length(setdiff(A_set, union(B_set, C_set))),
      B = length(setdiff(B_set, union(A_set, C_set))),
      C = length(setdiff(C_set, union(A_set, B_set)))
    )
    over_excl <- c(
      AB = length(setdiff(intersect(A_set, B_set), ABC_set)),
      BC = length(setdiff(intersect(B_set, C_set), ABC_set)),
      AC = length(setdiff(intersect(A_set, C_set), ABC_set)),
      ABC = length(ABC_set)
    )

    scale_factor <- 1
    radii <- sqrt(c(length(A_set), length(B_set), length(C_set)) / pi) * scale_factor
    r1 <- radii[1]; r2 <- radii[2]; r3 <- radii[3]

    # Advanced distance function
    get_dist <- function(count1, count2, shared, r1, r2) {
      if (shared <= 0) return((r1 + r2) * 1.1) # Slight gap if no overlap

      # Calculate overlap ratio relative to the smaller set
      overlap_ratio <- shared / min(count1, count2)

      # We use a non-linear power (0.5) so that even moderate overlap
      # pulls the circles significantly toward each other.
      # The 0.1 multiplier ensures they never perfectly stack unless ratio is 1.0
      d <- (r1 + r2) * (1 - (overlap_ratio^0.5 * 0.9))
      return(d)
    }

    #Calculate Distances
    d12 <- get_dist(length(A_set), length(B_set), over_excl[1], r1, r2)
    d23 <- get_dist(length(B_set), length(C_set), over_excl[2], r2, r3)
    d13 <- get_dist(length(A_set), length(C_set), over_excl[3], r1, r3)

    # Calculate Center Coordinates (Law of Cosines)
    if(d12 + d23 < d13) d13 <- d12 + d23
    if(d12 + d13 < d23) d23 <- d12 + d13

    cX <- (d12^2 + d13^2 - d23^2) / (2 * d12)
    cY <- sqrt(max(0, d13^2 - cX^2))

    centers <- data.frame(
      group = group_result,
      x = c(0, d12, cX),
      y = c(0, 0, cY),
      r = radii
    )

    #Create Polygons
    create_poly <- function(h, k, r, g) {
      theta <- seq(0, 2*pi, length.out = 100)
      data.frame(x = h + r*cos(theta), y = k + r*sin(theta), group = g)
    }

    poly_df <- do.call(rbind, lapply(1:3, function(i) {
      create_poly(centers$x[i], centers$y[i], centers$r[i], centers$group[i])
    }))

    # Label positions
    labels_df <- data.frame(
      count = c(vals_excl, over_excl),
      x = c(
        centers$x[1] - r1*0.4, centers$x[2] + r2*0.4, centers$x[3], # Unique A, B, C
        (centers$x[1] + centers$x[2])/2,                           # AB
        (centers$x[2] + centers$x[3])/2,                           # BC
        (centers$x[1] + centers$x[3])/2,                           # AC
        mean(centers$x)                                            # ABC
      ),
      y = c(
        centers$y[1], centers$y[2], centers$y[3] + r3*0.4,
        (centers$y[1] + centers$y[2])/2,
        (centers$y[2] + centers$y[3])/2,
        (centers$y[1] + centers$y[3])/2,
        mean(centers$y)
      )
    )

    # Plot
    ggplot2::ggplot() +
      ggplot2::geom_polygon(data = poly_df, ggplot2::aes(x = .data$x,
                                                         y = .data$y,
                                                         fill = .data$group),
                            alpha = 0.5, color = "white", linewidth = 0.5) +
      # Main Labels (Set Names)
      ggplot2::geom_text(data = centers, ggplot2::aes(x = .data$x,
                                                      y = .data$y + .data$r + (max(radii)*0.1),
                                                      label = .data$group),
                         fontface = "bold", size = 5) +
      ggplot2::geom_text(data = labels_df, ggplot2::aes(x = .data$x,
                                                        y = .data$y,
                                                        label = .data$count),
                         size = 4, color = "black") +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")
  }else if (length(group_result) == 2) {
    # Calculate overlap and difference
    group_list <- split(df$var, as.character(df$id))

    A_set <- group_list[[group_result[1]]]
    B_set <- group_list[[group_result[2]]]

    # Exclusive and Overlap sets
    AB_set <- intersect(A_set, B_set)

    vals_excl <- c(
      A = length(setdiff(A_set, B_set)),
      B = length(setdiff(B_set, A_set))
    )
    over_excl <- c(
      AB = length(AB_set)
    )

    # Radii based on total set sizes (Area = pi * r^2)
    r1 <- sqrt(length(A_set) / pi)
    r2 <- sqrt(length(B_set) / pi)
    target_overlap <- over_excl

    # Function to calculate area of overlap given distance d
    calc_overlap_area <- function(d, r1, r2) {
      if (d >= r1 + r2) return(0)         # No overlap
      if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2) # One inside the other

      part1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
      part2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
      part3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))

      return(part1 + part2 - part3)
    }

    # Find the exact distance 'd' using a root finder (uniroot)
    if (target_overlap == 0) {
      d_final <- r1 + r2
    } else if (target_overlap >= min(length(A_set), length(B_set))) {
      d_final <- abs(r1 - r2)
    } else {
      d_final <- stats::uniroot(function(x) calc_overlap_area(x, r1, r2) - target_overlap,
                         lower = abs(r1 - r2), upper = r1 + r2)$root
    }

    #Define Centers and Polygons
    centers <- data.frame(
      group = group_result,
      x = c(0, d_final),
      y = c(0, 0),
      r = c(r1, r2)
    )

    create_poly <- function(h, k, r, g) {
      theta <- seq(0, 2*pi, length.out = 100)
      data.frame(x = h + r*cos(theta), y = k + r*sin(theta), group = g)
    }

    poly_df <- do.call(rbind, lapply(1:2, function(i) {
      create_poly(centers$x[i], centers$y[i], centers$r[i], centers$group[i])
    }))

    # Dynamic Label Positioning
    labels_df <- data.frame(
      count = c(vals_excl[1], vals_excl[2], over_excl),
      x = c(centers$x[1] - (r1 * 0.4), # Push A label left
            centers$x[2] + (r2 * 0.4), # Push B label right
            (centers$x[1] + centers$x[2]) / 2), # Overlap in center
      y = 0
    )

    # Plot
    ggplot2::ggplot() +
      ggplot2::geom_polygon(data = poly_df,
                            ggplot2::aes(x = .data$x,
                                         y = .data$y,
                                         fill = .data$group),
                            alpha = 0.5, color = "white") +
      ggplot2::geom_text(data = centers,
                         ggplot2::aes(x = .data$x,
                                      y = .data$y + .data$r + 0.5,
                                      label = .data$group),
                         fontface = "bold") +
      ggplot2::geom_text(data = labels_df,
                         ggplot2::aes(x = .data$x,
                                      y = .data$y,
                                      label = .data$count)) +
      ggplot2::scale_fill_manual(values = .modEnv$colorScheme) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void()
  }
}
