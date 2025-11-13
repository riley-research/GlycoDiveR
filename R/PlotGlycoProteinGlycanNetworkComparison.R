PlotGlycoProteinGlycanNetworkComparison <- function(input, left, right, type = "N",
                                          edgeWidth = 1.5, verticeSize = c(5,3),
                                          whichAlias = NULL, highlight = "all",
                                          whichPeptide = NA){
  if(!is.na(condition)){
    input$PTMTable <- input$PTMTable[input$PTMTable[["Alias"]] %in% c(condition1,condition2) ]

    if(nrow(input$PTMTable) == 0){
      stop("Just tried selecting the following condition: ", condition1, " and ",
           condition2, ".\nNo rows are left after filtering.")
    }
  }

  input <- FilterForCutoffs(input)
  input$PTMTable <- FilterForPeptides(input$PTMTable, whichPeptide)

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

  View(dfglt)



}
