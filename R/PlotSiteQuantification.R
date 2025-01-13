PlotSiteQuantification <- function(input, protein, site, cutoff = NA){
  input <- FilterForCutoffs(input)

  df <- GetMeanTechReps(input$PTMTable)

  df <- subset(df, UniprotIDs == protein & ModificationID == site) %>%
    dplyr::ungroup()

  if(!is.na(cutoff)){
    if(substr(cutoff, nchar(cutoff), nchar(cutoff)) == "%"){
      ctfp <- as.numeric(substr(cutoff, 1, nchar(cutoff) -1))
      dftemp <- subset(df, Intensity > (ctfp/100) * max(Intensity, na.rm = TRUE))
      df <- subset(df, TotalGlycanComposition %in% dftemp$TotalGlycanComposition)
    }else{
      dftemp <- subset(df, Intensity > cutoff)
      df <- subset(df, TotalGlycanComposition %in% dftemp$TotalGlycanComposition)
    }
  }

  dfsum <- df %>%
    dplyr::group_by(Condition, TotalGlycanComposition) %>%
    dplyr::reframe(Condition = Condition, TotalGlycanComposition = TotalGlycanComposition,
                   mean = mean(Intensity, na.rm = TRUE), sd = sd(Intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  df$TotalGlycanComposition <- factor(df$TotalGlycanComposition, levels = unique(df$TotalGlycanComposition))

  p <- ggplot2::ggplot(data = dfsum) +
    ggplot2::geom_bar(data = dfsum, aes(x = TotalGlycanComposition, y = mean, fill = Condition), stat = "identity", color = "black") +
    ggplot2::geom_errorbar(data = dfsum, aes(x = TotalGlycanComposition, ymin = mean, ymax = mean+sd)) +
    ggplot2::geom_point(data = df, aes(x = TotalGlycanComposition, y = Intensity)) +
    ggplot2::labs(x = "", y = "Intensity (a.u.)") +
    ggplot2::facet_wrap(~Condition) +
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, max(df$Intensity)*1.05)) +
    ggplot2::scale_fill_manual(values = colorScheme)

  print(p)
}
