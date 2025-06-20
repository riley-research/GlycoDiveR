#' PlotGlycoproteinGlycanNetwork
#'
#' @param input Formatted data
#' @param condition What condition to plot as listed in the annotation dataframe
#' @param type Glycan type "N"
#'
#' @returns A glycoprotein to glycan network
#' @export
#'
#' @examples PlotGlycoProteinGlycanNetwork(mydata, condition = "DB_treated")
PlotGlycoProteinGlycanNetwork <- function(input, condition = NA, type = "N"){
  if(!is.na(condition)){
    input$PTMTable <- subset(input$PTMTable, Condition %in% condition)
  }

  input <- FilterForCutoffs(input)

  df <- subset(input$PTMTable, GlycanType != "NonGlyco" & GlycanType != "OGlycan")
  df <- df[,c("TotalGlycanComposition", "GlycanType", "UniprotIDs", "NumberOfNSites")] %>%
    dplyr::distinct()

  #Get the proteins, colors, and coordinates
  dfprot <- df[,c("UniprotIDs", "NumberOfNSites")] %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(NumberOfNSites)) %>%
    dplyr::mutate(x = 1)
  dfprot$y <- seq(1,nrow(dfprot))

  dfprot$col <- sapply(dfprot$NumberOfNSites, function(x) dplyr::case_when(x == 1 ~ "#FAFAFA",
                                                                           x == 2 ~ "#D9D9D9",
                                                                           x == 3 ~ "#BDBDBD",
                                                                           x == 4 ~ "#737574",
                                                                           x == 5 ~ "#535554",
                                                                           TRUE ~ "#1E1E1E"))

  dfleg <- dfprot %>%
    dplyr::arrange(NumberOfNSites) %>%
    dplyr::rowwise() %>%
    mutate(NumberOfNSites = ifelse(NumberOfNSites > 5, ">5", as.character(NumberOfNSites)))


  #Get the glycans and colors
  dfgly <- df[,c("TotalGlycanComposition", "GlycanType")] %>%
    dplyr::distinct() %>%
    dplyr::arrange(GlycanType)

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
  x_coords <- sapply(x_coords, function(x) ifelse(x <= 0.01, x - (0.2 * radius), x + (0.2 * radius)))

  y_coords <- radius * sin(angles)
  y_coords <- sapply(y_coords, function(x) x + radius)

  circle_coords <- data.matrix(data.frame(x = x_coords, y = y_coords))

  #Add both the a graph
  g <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g <- igraph::add_vertices(g, nv = nrow(dfprot), attr = list(name = dfprot$UniprotIDs))
  g <- igraph::add_vertices(g, nv = nrow(dfgly), attr = list(name = dfgly$TotalGlycanComposition))

  igraph::V(g)$shape <- c(rep("square", nrow(dfprot)), rep("circle", nrow(dfgly)))
  igraph::V(g)$size <- c(rep(5, nrow(dfprot)), rep(2, nrow(dfgly)))
  igraph::V(g)$color <- c(dfprot$col, dfgly$col)

  #Combine the coords
  coords <- rbind(data.matrix(dfprot[,c("x", "y")]), circle_coords)

  #Now add color to the edge df and add the edges
  df <- df %>%
    dplyr::left_join(by = "GlycanType", dfcolmatch)

  df$colv <- sapply(df$col, function(x) paste0(x, as.character(50)))

  edge_list <- as.vector(t(df[,c("TotalGlycanComposition", "UniprotIDs")]))

  g <- igraph::add_edges(g, edges = edge_list)

  plot(g, layout = coords, vertex.label = NA, vertex.size = igraph::V(g)$size,
       vertex.color = igraph::V(g)$color, edge.color = df$colv)

  legend('topright', legend = dfcolmatch$GlycanType, fill = dfcolmatch$col, title="Glycan Type")
  legend('bottomright', legend = unique(dfleg$NumberOfNSites), fill = unique(dfleg$col), title = "Sites Per\nProtein")
}
