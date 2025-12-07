#' PlotGlycoproteinGlycanNetwork
#'
#' Creates a network where each node in the circle represents one glycan
#' composition (e.g., H1N1, or H1N1A), connecting to a protein represented as
#' a bar plot and sorted by number of identified glycan sites.
#'
#' @param input Formatted data imported through a GlycoDiveR importer.
#' @param edgeWidth Defines the linewidth of the edges.
#' @param verticeSize Defines the size of the glycan and protein nodes,
#' verticeSize = c(3,5) means a glycan node size of 3 and a protein node size of 5.
#' @param whichAlias Provide a vector of Aliases to only select these aliases
#' for plotting.
#' @param highlight Specify what glycan category to highlight, or use "all" to highlight all.
#' @param whichPeptide Filter what peptides to plot. This can either be a dataframe
#' with a ModifiedPeptide peptide column, or a vector with the ModifiedPeptide sequences
#' that you want to keep. Inputted data with the comparison importer functions is
#' directly usable, also after filtering using the FilterComparison function.
#' @param whichProtein Filter what proteins to plot. These are the IDs as presented
#' in the UniprotIDs column in your GlycoDiveR data. This can either be a dataframe
#' with a UniprotIDs column, or a vector with the UniprotIDs you want to keep.
#' @param exactProteinMatch This is only relevant if you select for proteins using
#' the whichProtein argument. When set to TRUE (default), your supplied UniprotIDs
#' must be an exact match to the UniprotIDs in the dataframe. When set to FALSE,
#' it will select non-exact matches. For example, "P61224" will only match to
#' "P61224,P62834" when set to FALSE.
#' @param silent silence printed information (default = FALSE).
#'
#' @returns A glycoprotein to glycan network
#' @export
#'
#' @examples
#' \dontrun{
#' PlotGlycoProteinGlycanNetwork(mydata)
#'
#' PlotGlycoProteinGlycanNetwork(
#'   mydata,
#'   highlight = c("Sialyl", "Truncated"),
#'   edgeWidth = 3,
#'   verticeSize = c(3, 4)
#' )
#' }
PlotGlycoProteinGlycanNetwork <- function(input, edgeWidth = 1.5,
                                          verticeSize = c(5,3), whichAlias = NULL,
                                          highlight = "all", whichPeptide = NULL,
                                          whichProtein = NULL, exactProteinMatch = TRUE,
                                          silent = FALSE){
  input <- FilterForCutoffs(input, silent)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)
  input$PTMTable <- FilterForProteins(input$PTMTable, whichProtein, exactProteinMatch)

  df <- input$PTMTable %>%
    dplyr::filter(.data$GlycanType != "NonGlyco")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
  }

  if(nrow(df) == 0){
    if(!silent){
      return(fmessage("No data is left after filtering."))
    }else{
      return()
    }
  }

  df <- df[,c("TotalGlycanComposition", "GlycanType", "UniprotIDs", "NumberOfNSites")] %>%
    dplyr::distinct()

  #Get the proteins, colors, and coordinates
  dfprot <- df[,c("UniprotIDs", "NumberOfNSites")] %>%
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::desc(.data$NumberOfNSites)) %>%
    dplyr::mutate(x = 1)
  dfprot$y <- seq(1,nrow(dfprot))

  dfprot$col <- sapply(dfprot$NumberOfNSites, function(x) dplyr::case_when(x == 1 ~ "grey90",
                                                                           x == 2 ~ "grey70",
                                                                           x == 3 ~ "grey50",
                                                                           x == 4 ~ "grey30",
                                                                           x == 5 ~ "grey15",
                                                                           TRUE ~ "black"))

  dfleg <- dfprot %>%
    dplyr::arrange(.data$NumberOfNSites) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(NumberOfNSites = ifelse(.data$NumberOfNSites > 5, ">5", as.character(.data$NumberOfNSites))) %>%
    dplyr::ungroup()

  #Get the glycans and colors
  dfgly <- df[,c("TotalGlycanComposition", "GlycanType")] %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$GlycanType)

  dfcolmatch <- .modEnv$GlycanColors %>%
    dplyr::rename("col" = "color") %>%
    dplyr::filter(.data$GlycanType %in% unique(dfgly$GlycanType))

  dfgly <- dfgly %>%
    dplyr::left_join(by = "GlycanType", dfcolmatch)

  #Get the circular coordinates
  numv <- nrow(dfgly)
  radius <- nrow(dfprot) * 0.5

  angles <- seq(pi / 2, 2 * pi + pi / 2, length.out = numv + 1)[-numv - 1]

  x_coords <- radius * cos(angles)
  x_coords <- sapply(x_coords, function(x) ifelse(x <= 0.01,
                                                  x - (0.2 * radius),
                                                  x + (0.2 * radius)))

  y_coords <- radius * sin(angles)
  y_coords <- sapply(y_coords, function(x) x + radius)

  circle_coords <- data.matrix(data.frame(x = x_coords, y = y_coords))

  # --- Protein coordinates (centered) ---
  center_y <- radius
  nprot <- nrow(dfprot)

  if (nprot == 1) {
    # Only one protein → put it exactly in the center
    prot_y <- center_y
  } else {
    # Multiple proteins → spread them evenly with the middle at center_y
    prot_y <- seq(center_y - (nprot - 1) / 2,
                  center_y + (nprot - 1) / 2,
                  length.out = nprot)
  }

  negX <- circle_coords[,"x"][circle_coords[,"x"] < 0]
  posX <- circle_coords[,"x"][circle_coords[,"x"] > 0]

  negX <- negX[which.max(negX)]
  posX <- posX[which.min(posX)]

  dfprot$x <- mean(c(negX, posX))
  dfprot$y <- prot_y

  circle_coords <- data.matrix(data.frame(x = x_coords, y = y_coords))

  #Add both the a graph
  g <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g <- igraph::add_vertices(g, nv = nrow(dfprot), attr = list(name = dfprot$UniprotIDs, frame.color = "transparent"))
  g <- igraph::add_vertices(g, nv = nrow(dfgly), attr = list(name = dfgly$TotalGlycanComposition, frame.color = "black"))

  igraph::V(g)$shape <- c(rep("square", nrow(dfprot)), rep("circle", nrow(dfgly)))
  igraph::V(g)$size <- c(rep(verticeSize[1], nrow(dfprot)), rep(verticeSize[2], nrow(dfgly)))
  igraph::V(g)$color <- c(dfprot$col, dfgly$col)

  #Combine the coords
  coords <- rbind(data.matrix(dfprot[,c("x", "y")]), circle_coords)

  #Now add color to the edge df and add the edges
  df <- df %>%
    dplyr::left_join(by = "GlycanType", dfcolmatch)

  if(highlight == "all"){
    df$colv <- sapply(df$col, function(x) paste0(x, as.character(80)))
  }else{
    df <- df %>%
      dplyr::mutate(colv = dplyr::if_else(.data$GlycanType %in% highlight,
                                          paste0(.data$col),
                                          NA))
  }

  edge_list <- as.vector(t(df[,c("TotalGlycanComposition", "UniprotIDs")]))

  g <- igraph::add_edges(g, edges = edge_list)

  plot(g, layout = coords, vertex.label = NA, vertex.size = igraph::V(g)$size,
       vertex.color = igraph::V(g)$color, edge.color = df$colv,
       edge.width = edgeWidth)

  graphics::legend(#"topright", inset = c(-0.18, 0),
    x=1.05,
    y=1,
    legend = dfcolmatch$GlycanType,
    fill = dfcolmatch$col,
    title = "Glycan Type",
    cex = 0.8,
    pt.cex = 0.6,
    bty = "n",
    xpd = NA)

  graphics::legend(#"bottomright", inset = c(0, 0),
    x = 1.05,
    y= 0,
         legend = unique(dfleg$NumberOfNSites),
         fill = unique(dfleg$col),
         title = "Sites Per\nProtein",
         cex = 0.8,
         pt.cex = 0.6,
         bty = "n",
         xpd = NA)

  p <- grDevices::recordPlot()
  return(p)
}
