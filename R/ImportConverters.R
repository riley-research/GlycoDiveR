MSFraggerConverter <- function(unfiltereddf, annotationdf, fastaPath, quantdf, scrape,
                               normalization, convertFPModCodeToMass, OPairLevelConversion){
  fmessage("Now starting import.")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  #Run####
  if("Spectrum.File" %in% existingCols){
    unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][length(strsplit(x, "\\", fixed = T)[[1]])-1])
    filtereddf <- cbind(filtereddf, Run = as.character(unfiltereddf$Run))
    fmessage("Successfully imported Run column.")
  }
  else{stop("The column Run was not found in the input dataframe.")}

  #ModifiedPeptide####
  if("Modified.Peptide" %in% existingCols){
    filtereddf$ModifiedPeptide <- apply(unfiltereddf[,c("Peptide", "Modified.Peptide")], 1, function(x) GetPeptide(pep = x[1], modpep = x[2]))
    fmessage("Successfully imported Modified Peptide column.")}
  else{stop("The column Modified.Peptide was not found in the input dataframe.")}

  #HasNSequon####
  if("Has.N.Glyc.Sequon" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(HasNSequon = unfiltereddf$Has.N.Glyc.Sequon)
    fmessage("Successfully imported HasNSequon column.")}
  else{fmessage("No NSequon information(OPair only) in data, skipping")}

  #PSMScore####
  if ("Hyperscore" %in% existingCols) {
    filtereddf <- cbind(filtereddf, PSMScore = as.double(unfiltereddf$Hyperscore))
    fmessage("Successfully imported Hyperscore column and renamed to PSMScore.")}
  else{stop("The column Hyperscore was not found in the input dataframe.")}

  #AssignedModifications####
  if ("Assigned.Modifications" %in% existingCols) {
    filtereddf <- cbind(filtereddf, AssignedModifications = as.character(unfiltereddf$Assigned.Modifications))
    fmessage("Successfully imported Assigned Modifications column.")}
  else {stop("The column Assigned.Modifications was not found in the input dataframe.")}

  #ModifiedPeptide clean####
  if(convertFPModCodeToMass){
    filtereddf$ModifiedPeptide <- FPModCodeToModMass(filtereddf$ModifiedPeptide, filtereddf$AssignedModifications)
    fmessage("Successfully cleaned ModifiedPeptide column.")
  }

  #TotalGlycanComposition####
  if ("Total.Glycan.Composition" %in% existingCols) {
    filtereddf <- cbind(filtereddf, TotalGlycanComposition = as.character(unfiltereddf$Total.Glycan.Composition))
    filtereddf$TotalGlycanComposition <- CleanGlycanNames(filtereddf$TotalGlycanComposition)
    fmessage("Successfully imported Total Glycan Composition column.")}
  else {stop("The column Total.Glycan.Composition was not found in the input dataframe.")}

  #GlycanQValue####
  if ("Glycan.q.value" %in% existingCols & !("Confidence.Level" %in% existingCols)) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(GlycanQValue = as.double(unfiltereddf$Glycan.q.value))
    fmessage("Successfully imported Glycan q Value column.")
  }else if(!("Glycan.q.value" %in% existingCols) & "Confidence.Level" %in% existingCols){
    filtereddf <- filtereddf %>%
      dplyr::mutate(GlycanQValue = unfiltereddf$Confidence.Level)

    filtereddf$GlycanQValue <- dplyr::case_when(filtereddf$GlycanQValue == "Level1" ~ as.character(OPairLevelConversion[1]),
                                                filtereddf$GlycanQValue == "Level1b" ~ as.character(OPairLevelConversion[2]),
                                                filtereddf$GlycanQValue == "Level2" ~ as.character(OPairLevelConversion[3]),
                                                filtereddf$GlycanQValue == "Level3" ~ as.character(OPairLevelConversion[4]),
                                                TRUE ~ filtereddf$GlycanQValue)

    filtereddf$GlycanQValue <- as.numeric(filtereddf$GlycanQValue)
    fmessage("Successfully imported Glycan q Value column.")
    }else if("Glycan.q.value" %in% existingCols & "Confidence.Level" %in% existingCols){
      filtereddf <- filtereddf %>%
        dplyr::mutate(GlycanQValue = as.double(unfiltereddf$Glycan.q.value),
                      ConfidenceLevel = unfiltereddf$Confidence.Level)

      filtereddf$GlycanQValue <- as.double(dplyr::case_when(is.na(filtereddf$GlycanQValue) & filtereddf$ConfidenceLevel == "Level1" ~ OPairLevelConversion[1],
                                                  is.na(filtereddf$GlycanQValue) & filtereddf$ConfidenceLevel == "Level1b" ~ OPairLevelConversion[2],
                                                  is.na(filtereddf$GlycanQValue) & filtereddf$ConfidenceLevel == "Level2" ~ OPairLevelConversion[3],
                                                  is.na(filtereddf$GlycanQValue) & filtereddf$ConfidenceLevel == "Level3" ~ OPairLevelConversion[4],
                                                  TRUE ~ filtereddf$GlycanQValue))

      filtereddf <- filtereddf %>%
        dplyr::select(-"ConfidenceLevel")
      fmessage("Successfully imported Glycan q Value column.")
    }
  else {stop("The column Glycan.q.value was not found in the input dataframe.")}

  #IsUnique####
  if ("Is.Unique" %in% existingCols) {
    filtereddf <- cbind(filtereddf, IsUnique = as.logical(unfiltereddf$Is.Unique))
    fmessage("Successfully imported Is Unique column.")}
  else {stop("The column Is.Unique was not found in the input dataframe.")}

  #UniprotIDs####
  if ("Protein.ID" %in% existingCols) {
    if ("Mapped.Proteins" %in% existingCols) {
      filtereddf$UniprotIDs <- apply(unfiltereddf[, c("Protein.ID", "Mapped.Proteins")], 1, function(x) CombineProtCols(x))
      fmessage("Successfully imported and merged Protein ID and Mapped protein columns.")
    } else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Proteins column found, only taking from Protein ID column.")
    }
  } else {
    stop("The column Protein.ID was not found in the input dataframe.")}

  #Genes####
  if ("Gene" %in% existingCols) {
    if("Mapped.Genes" %in% existingCols){
      filtereddf$Genes <- apply(unfiltereddf[, c("Gene", "Mapped.Genes")], 1, function(x) CombineTwoCols(x))
      fmessage("Successfully imported and merged Gene and Mapped Genes columns.")
    }else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Genes columns found, only taking from Gene column")}
  } else {stop("The column Gene was not found in the input dataframe.")}

  #ProteinLength####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = .data$UniprotIDs,
                      ProteinLength = GetProteinLength(IDVec = .data$UniprotIDs,
                                                       fastaFile = fastaFile))
    }else{warning("Fasta path does not exist.")}
  }

  #NumberOfNSites/NumberOfOSites####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = (.data$UniprotIDs),
                      NumberOfSites = GetGlycoSitesPerProtein(IDVec = .data$UniprotIDs,
                                                              fastaFile = fastaFile))

      filtereddf <- filtereddf %>%
        tidyr::separate_wider_delim(.data$NumberOfSites, delim = ";", names = c("NumberOfNSites", "NumberOfOSites")) %>%
        dplyr::mutate(NumberOfNSites = as.numeric(.data$NumberOfNSites),
                      NumberOfOSites = as.numeric(.data$NumberOfOSites))
      fmessage("Successfully mapped number of N and O glycosites per protein.")
    }else{warning("Fasta path does not exist.")}
  }

  #ProteinStart####
  if ("Protein.Start" %in% existingCols) {
    filtereddf <- cbind(filtereddf, ProteinStart = as.numeric(unfiltereddf$Protein.Start))}
  else {filtereddf$ProteinStart = NA
    warning("The column Protein Start was not found in the input dataframe.")}

  #RetentionTime####
  if ("Retention" %in% existingCols) {
    filtereddf <- cbind(filtereddf, RetentionTime = as.numeric(unfiltereddf$Retention)/60)}
  else {filtereddf$RetentionTime = NA
  warning("The column Retention was not found in the input dataframe.")}

  #GlycanType####
  filtereddf$GlycanType <- apply(filtereddf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))
  filtereddf <- filtereddf %>%
    dplyr::mutate(GlycanType = sapply(.data$GlycanType, toString))
  fmessage("Successfully added GlycanType column.")

  if(scrape){
    fmessage("Now scraping Uniprot. Set 'scrape = FALSE' to the importer function to skip this step.
             Grab a coffee. This might take a while.")

    #Get subcellular localization and domain information####
    filtereddf <- filtereddf %>%
      dplyr::bind_cols(GetUniprotSubcellularInfo(filtereddf$UniprotIDs))
    fmessage("Successfully added subcellular localization and domain information.")
  }else{
    filtereddf$SubcellularLocalization <- NA
    filtereddf$Domains <- NA
  }

  #RawIntensity####
  if("Intensity" %in% existingCols) {
    filtereddf <- cbind(filtereddf, RawIntensity = as.numeric(unfiltereddf$Intensity))
    fmessage("Successfully imported RawIntensity column.")}
  else{filtereddf <- cbind(filtereddf, RawIntensity = as.numeric(NA))
  warning("Intensity column not found. Filled with NA.")}

  #Intensity####
  filtereddf$Intensity <- filtereddf$RawIntensity
  if(!all(is.na(filtereddf$Intensity)) & sum(filtereddf$Intensity != 0)){
    if(normalization == "none"){
      fmessage("Successfully imported Intensity column without normalization.")
    }else if(normalization == "median"){
      globalMedian = stats::median(filtereddf$Intensity[filtereddf$Intensity != 0], na.rm = TRUE)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = .data$Run,
                      Intensity = medianNormalization(intensityVec = .data$Intensity,
                                                      globalMedian = globalMedian))
      fmessage("Successfully median normalized the intensities.")
    }else if(normalization %in% c("FP_Normalized", "FP_MaxLFQ")){
      filtereddf <- UpdateFPIntensities(filtereddf, quantdf, normalization)
      fmessage(paste0("Successfully imported Intensity column using: ", normalization))
    }
  }else{
    fmessage("Successfully imported Intensity column. Note: No quantitative values found.")
  }

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf, by = "Run")

  filtereddf$Alias <- factor(filtereddf$Alias, levels = unique(annotationdf$Alias))

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}

ByonicConverter <- function(unfiltereddf, annotationdf, fastaPath,
                            modification_df, scrape){
  fmessage("Now starting import.")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  #Run####
  if("Run" %in% existingCols){
    filtereddf <- filtereddf %>%
      dplyr::mutate(Run = as.character(unfiltereddf$Run))
    fmessage("Successfully imported Run column.")
  }
  else{stop("The column Run was not found in the input dataframe.")}

  #ModifiedPeptide####
  if("Peptide\r\n< ProteinMetrics Confidential >" %in% existingCols){
    filtereddf <- filtereddf %>%
      dplyr::mutate(ModifiedPeptide = as.character(unfiltereddf$`Peptide\r\n< ProteinMetrics Confidential >`),
                    ModifiedPeptide = stringr::str_extract(.data$ModifiedPeptide, "(?<=\\.).*(?=\\.[^\\.]*$)"),
                    ModifiedPeptide = gsub("\\+", "", .data$ModifiedPeptide))

    fmessage("Successfully imported Peptide ProteinMetrics Confidential  column.")}
  else{stop("The column Peptide ProteinMetrics Confidential was not found in the input dataframe.")}

  #RawIntensity & Intensity####
  filtereddf <- filtereddf %>%
    dplyr::mutate(RawIntensity = 0,
                  Intensity = 0)
  fmessage("No intensity values, so proceeding with all NA.")

  #PSMScore####
  if ("Score" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(PSMScore = as.double(unfiltereddf$Score))
    fmessage("Successfully imported Score column and renamed to PSMScore.")}
  else{stop("The column Score was not found in the input dataframe.")}

  #AssignedModifications####
  if ("Modification Type(s)" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(AssignedModifications = as.character(unfiltereddf$`Modification Type(s)`))

    filtereddf$AssignedModifications <- ConvertByonicAssignedModifications(filtereddf$ModifiedPeptide, filtereddf$AssignedModifications)
    fmessage("Successfully imported Modification Type(s) column.")}
  else {stop("The column Modification Type(s) was not found in the input dataframe.")}

  #TotalGlycanComposition####
  if("Glyco" %in% modification_df$ModificationType){
    tempMod <- modification_df %>%
      dplyr::filter(.data$ModificationType == "Glyco") %>%
      dplyr::mutate(Mass = round(as.double(.data$Mass), 3))

    filtereddf <- filtereddf %>%
      dplyr::mutate(TotalGlycanComposition = AssignedModsToGlycanComp_Byonic(.data$AssignedModifications,tempMod))
    fmessage("Successfully imported Total Glycan Composition column.")
  }else{
    filtereddf$TotalGlycanComposition <- NA
    fmessage("No glyco found in the glyco table. All glyco are listed as NA")
  }

  #GlycanQValue####
  if ("|Log Prob|" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(GlycanQValue = 10^(-as.double(unfiltereddf$`|Log Prob|`)))
    fmessage("Successfully imported |Log Prob| (as Glycan Q Value) Value column.")
  }else{
    filtereddf$GlycanQValue <- NA
    fmessage("No |Log Prob| found.")
  }

  #UniprotIDs####
  if ("Protein Name" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(UniprotIDs = stringr::str_extract(unfiltereddf$`Protein Name`, "(?<=\\|).*?(?=\\|)"))
    fmessage("Successfully imported Protein Name column.")
  } else {
    stop("The column Protein.ID was not found in the input dataframe.")}

  #Genes####
  if ("Protein Name" %in% existingCols) {
    filtereddf$Genes <- sapply(unfiltereddf$`Protein Name`, function(x) {
      temp <- strsplit(x, "\\|")[[1]][3]
      temp <-  strsplit(temp, "_")[[1]][1]
      return(temp)
    })
    fmessage("Successfully imported and merged Gene and Mapped Genes columns.")
  } else {stop("The column Protein Name was not found in the input dataframe.")}

  #ProteinLength####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = "UniprotIDs",
                      ProteinLength = GetProteinLength(IDVec = .data$UniprotIDs,
                                                       fastaFile = fastaFile))
    }else{warning("Fasta path does not exist.")}
  }

  #NumberOfNSites/NumberOfOSites####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = "UniprotIDs",
                      NumberOfSites = GetGlycoSitesPerProtein(IDVec = .data$UniprotIDs,
                                                              fastaFile = fastaFile))

      filtereddf <- filtereddf %>%
        tidyr::separate_wider_delim(.data$NumberOfSites, delim = ";", names = c("NumberOfNSites", "NumberOfOSites")) %>%
        dplyr::mutate(NumberOfNSites = as.numeric(.data$NumberOfNSites),
                      NumberOfOSites = as.numeric(.data$NumberOfOSites))
      fmessage("Successfully mapped number of N and O glycosites per protein.")
    }else{warning("Fasta path does not exist.")}
  }

  #ProteinStart####
  if ("Starting\r\nposition" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(ProteinStart = as.numeric(unfiltereddf$`Starting\r\nposition`))
    }else {filtereddf$ProteinStart = NA
  warning("The column Protein Start was not found in the input dataframe.")}

  #RetentionTime####
  if ("Scan Time" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(RetentionTime = as.double(unfiltereddf$`Scan Time`))
  }else {filtereddf$RetentionTime = NA
  warning("The column Scan Time was not found in the input dataframe.")}

  #GlycanType####
  filtereddf$GlycanType <- apply(filtereddf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))
  filtereddf <- filtereddf %>%
    dplyr::mutate(GlycanType = sapply(.data$GlycanType, toString))
  fmessage("Successfully added GlycanType column.")

  if(scrape){
    fmessage("Now scraping Uniprot. Set 'scrape = FALSE' to the importer function to skip this step.
             Grab a coffee. This might take a while.")

    #Get subcellular localization and domain information####
    filtereddf <- filtereddf %>%
      dplyr::bind_cols(GetUniprotSubcellularInfo(filtereddf$UniprotIDs))
    fmessage("Successfully added subcellular localization and domain information.")
  }else{
    filtereddf$SubcellularLocalization <- NA
    filtereddf$Domains <- NA
  }

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf, by = "Run")

  filtereddf$Alias <- factor(filtereddf$Alias, levels = unique(annotationdf$Alias))

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}

pGlycoConverter <- function(unfiltereddf, annotationdf, fastaPath,
                            modification_df, normalization, scrape){
  fmessage("Now starting import.")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  #Run####
  if("RawName" %in% existingCols){
    filtereddf <- filtereddf %>%
      dplyr::mutate(Run = as.character(unfiltereddf$RawName))
    fmessage("Successfully imported Run column.")
  }
  else{stop("The column Run was not found in the input dataframe.")}

  #ModifiedPeptide####
  if("Peptide" %in% existingCols & "Mod" %in% existingCols){
    filtereddf$ModifiedPeptide <- GetModifiedPeppGlyco(unfiltereddf$Peptide, unfiltereddf$Mod, modification_df)
    fmessage("Successfully imported Peptide ProteinMetrics Confidential  column.")}
  else{stop("The column Peptide ProteinMetrics Confidential was not found in the input dataframe.")}

  #RawIntensity####
  if("IsotopeArea" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(RawIntensity = as.numeric(unfiltereddf$IsotopeArea))
    fmessage("Successfully imported RawIntensity column.")}
  else{filtereddf <- cbind(filtereddf, RawIntensity = as.numeric(0))
  warning("Intensity column not found. Filled with NA.")}

  #Intensity####
  filtereddf$Intensity <- filtereddf$RawIntensity
  if(!all(is.na(filtereddf$Intensity)) & sum(filtereddf$Intensity != 0)){
    if(normalization == "none"){
      fmessage("Successfully imported Intensity column without normalization.")
    }else if(normalization == "median"){
      globalMedian = stats::median(filtereddf$Intensity[filtereddf$Intensity != 0], na.rm = TRUE)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = .data$Run,
                      Intensity = medianNormalization(intensityVec = .data$Intensity,
                                                      globalMedian = globalMedian))
      fmessage("Successfully median normalized the intensities.")
    }
  }else{
    fmessage("Successfully imported Intensity column. Note: No quantitative values found.")
  }

  #PSMScore####
  if ("PeptideFDR" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(PSMScore = as.double(unfiltereddf$PeptideFDR))
    fmessage("Successfully imported Score column and renamed to PSMScore.")}
  else{stop("The column Score was not found in the input dataframe.")}

  #AssignedModifications####
  if (all(c("Mod", "GlyMass", "GlySite", "Peptide") %in% existingCols)) {
    filtereddf$AssignedModifications <- AssignmedModspGlyco(unfiltereddf$Mod,
                                                            unfiltereddf$GlyMass,
                                                            unfiltereddf$GlySite,
                                                            unfiltereddf$Peptide,
                                                            modification_df)
    fmessage("Successfully imported AssignedModifications column.")}
  else {stop("One or more of these columns were not found: Mod, GlyMass, GlySite, Peptide.")}

  #TotalGlycanComposition####
  if("GlycanComposition" %in% existingCols){
    filtereddf <- filtereddf %>%
      dplyr::mutate(TotalGlycanComposition = gsub(";", ",", unfiltereddf$GlycanComposition))
    fmessage("Successfully imported Total Glycan Composition column.")
  }else{
    filtereddf$TotalGlycanComposition <- NA
    fmessage("No glyco found in the glyco table. All glyco are listed as NA")
  }

  #GlycanQValue####
  if ("GlycanFDR" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(GlycanQValue = as.double(unfiltereddf$GlycanFDR))
    fmessage("Successfully imported GlycanFDR (as Glycan Q Value) Value column.")
  }else{
    filtereddf$GlycanQValue <- NA
    fmessage("No GlycanFDR found.")
  }

  #UniprotIDs####
  if ("Proteins" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(UniprotIDs = stringr::str_extract(unfiltereddf$Proteins, "(?<=\\|).*?(?=\\|)"))
    fmessage("Successfully imported UniprotIDs column.")
  } else {
    stop("The column Proteins was not found in the input dataframe.")}

  #Genes####
  if ("Genes" %in% existingCols) {
    filtereddf$Genes <- apply(unfiltereddf[c("Proteins", "Genes")], 1, function(x)
                        ifelse(x[2] == "" | is.na(x[2]),
                               strsplit(strsplit(x[1], "\\|")[[1]][3], "_")[[1]][1],
                               toupper(x[2])))
    fmessage("Successfully imported genes column.")
  } else {stop("The column Genes was not found in the input dataframe.")}

  #ProteinLength####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = "UniprotIDs",
                      ProteinLength = GetProteinLength(IDVec = .data$UniprotIDs,
                                                       fastaFile = fastaFile))
    }else{stop("Fasta path does not exist.")}
  }

  #NumberOfNSites/NumberOfOSites####
  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = "UniprotIDs",
                      NumberOfSites = GetGlycoSitesPerProtein(IDVec = .data$UniprotIDs,
                                                              fastaFile = fastaFile))

      filtereddf <- filtereddf %>%
        tidyr::separate_wider_delim(.data$NumberOfSites, delim = ";", names = c("NumberOfNSites", "NumberOfOSites")) %>%
        dplyr::mutate(NumberOfNSites = as.numeric(.data$NumberOfNSites),
                      NumberOfOSites = as.numeric(.data$NumberOfOSites))
      fmessage("Successfully mapped number of N and O glycosites per protein.")
    }else{warning("Fasta path does not exist.")}
  }

  #ProteinStart####
  if ("ModifiedPeptide" %in% names(filtereddf)) {
    fastaFile <- seqinr::read.fasta(file = fastaPath)
    filtereddf <- filtereddf %>%
      dplyr::mutate(.by = c("ModifiedPeptide", "UniprotIDs"),
        ProteinStart = GetPeptideLocInProtein(.data$UniprotIDs, .data$ModifiedPeptide, fastaFile))
  }else {filtereddf$ProteinStart = NA
  warning("The column Protein Start was not found in the input dataframe.")}

  #RetentionTime####
  if ("RT" %in% existingCols) {
    filtereddf <- filtereddf %>%
      dplyr::mutate(RetentionTime = as.double(unfiltereddf$RT / 60))
  }else {filtereddf$RetentionTime = NA
  warning("The column Scan Time was not found in the input dataframe.")}

  #GlycanType####
  filtereddf$GlycanType <- apply(filtereddf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))
  filtereddf <- filtereddf %>%
    dplyr::mutate(GlycanType = sapply(.data$GlycanType, toString))
  fmessage("Successfully added GlycanType column.")

  if(scrape){
    fmessage("Now scraping Uniprot. Set 'scrape = FALSE' to the importer function to skip this step.
             Grab a coffee. This might take a while.")

    #Get subcellular localization and domain information####
    filtereddf <- filtereddf %>%
      dplyr::bind_cols(GetUniprotSubcellularInfo(filtereddf$UniprotIDs))
    fmessage("Successfully added subcellular localization and domain information.")
  }else{
    filtereddf$SubcellularLocalization <- NA
    filtereddf$Domains <- NA
  }

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf, by = "Run")

  filtereddf$Alias <- factor(filtereddf$Alias, levels = unique(annotationdf$Alias))

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}
