#' PlotGlycoproteinGlycanNetwork
#'
#' @param input Formatted data
#' @param condition What condition to plot as listed in the annotation dataframe
#' @param type Glycan type "N"
#' @param edgeWidth defines the linewidth of the edges
#' @param verticeSize defines the size of the glycan and protein
#' nodes c(3,5) means a glycan node size of 3 and a protein node size of 5
#' @param whichAlias provide a vector of Aliases to only select these aliases
#' for plotting
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
#'   condition = "DB_treated",
#'   edgeWidth = 3,
#'   verticeSize = c(3, 4)
#' )
#' }
PlotGlycoProteinGlycanNetwork <- function(input, condition = NA, type = "N",
                                          edgeWidth = 1.5, verticeSize = c(5,3),
                                          whichAlias = NULL){
  if(!is.na(condition)){
    input$PTMTable <- input$PTMTable[input$PTMTable[["Condition"]] %in% condition, ]

    if(nrow(input$PTMTable) == 0){
      stop("Just tried selecting the following condition: ", condition, ".\nNo rows are left after filtering.")
    }
  }

  input <- FilterForCutoffs(input)

  df <- input$PTMTable %>%
    dplyr::filter(.data$GlycanType != "NonGlyco" & .data$GlycanType != "OGlycan")

  if(!is.null(whichAlias)){
    df <- df %>%
      dplyr::filter(.data$Alias %in% whichAlias)
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

  dfcolmatch <- data.frame(GlycanType = unique(df$GlycanType),
                           col = NA)

  for(i in 1:nrow(dfcolmatch)){
    dfcolmatch$col[i] <- colorScheme[i]
  }

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

  dfprot$x <- 0  # all proteins in a vertical column
  dfprot$y <- prot_y

  circle_coords <- data.matrix(data.frame(x = x_coords, y = y_coords))

  #Add both the a graph
  g <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g <- igraph::add_vertices(g, nv = nrow(dfprot), attr = list(name = dfprot$UniprotIDs))
  g <- igraph::add_vertices(g, nv = nrow(dfgly), attr = list(name = dfgly$TotalGlycanComposition))

  igraph::V(g)$shape <- c(rep("square", nrow(dfprot)), rep("circle", nrow(dfgly)))
  igraph::V(g)$size <- c(rep(verticeSize[1], nrow(dfprot)), rep(verticeSize[2], nrow(dfgly)))
  igraph::V(g)$color <- c(dfprot$col, dfgly$col)

  #Combine the coords
  coords <- rbind(data.matrix(dfprot[,c("x", "y")]), circle_coords)

  #Now add color to the edge df and add the edges
  df <- df %>%
    dplyr::left_join(by = "GlycanType", dfcolmatch)

  df$colv <- sapply(df$col, function(x) paste0(x, as.character(80)))

  edge_list <- as.vector(t(df[,c("TotalGlycanComposition", "UniprotIDs")]))

  g <- igraph::add_edges(g, edges = edge_list)

  plot(g, layout = coords, vertex.label = NA, vertex.size = igraph::V(g)$size,
       vertex.color = igraph::V(g)$color, edge.color = df$colv,
       edge.width = edgeWidth)

  graphics::legend("topright", inset = c(-0.18, 0),
         legend = dfcolmatch$GlycanType,
         fill = dfcolmatch$col,
         title = "Glycan Type",
         cex = 0.8,
         pt.cex = 0.6,
         bty = "n")

  graphics::legend("bottomright", inset = c(0, 0),
         legend = unique(dfleg$NumberOfNSites),
         fill = unique(dfleg$col),
         title = "Sites Per\nProtein",
         cex = 0.8,
         pt.cex = 0.6,
         bty = "n")
}
