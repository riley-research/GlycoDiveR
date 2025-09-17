MSFraggerConverter <- function(unfiltereddf, annotationdf, fastaPath){
  fmessage("Now starting import.")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  if("Run" %in% existingCols){
    filtereddf <- cbind(filtereddf, Run = as.character(unfiltereddf$Run))
    fmessage("Successfully imported Run column.")
  }
  else{stop("The column Run was not found in the input dataframe.")}

  if("Modified.Peptide" %in% existingCols){
    filtereddf$ModifiedPeptide <- apply(unfiltereddf[,c("Peptide", "Modified.Peptide")], 1, function(x) GetPeptide(pep = x[1], modpep = x[2]))
    fmessage("Successfully imported Modified Peptide column.")}
  else{stop("The column Modified.Peptide was not found in the input dataframe.")}

  if("Intensity" %in% existingCols) {
    filtereddf <- cbind(filtereddf, Intensity = as.numeric(unfiltereddf$Intensity))
    fmessage("Successfully imported Intensity column.")}
  else{filtereddf <- cbind(filtereddf, Intensity = as.numeric(NA))
       warning("Intensity column not found. Filled with NA.")}

  if ("Hyperscore" %in% existingCols) {
    filtereddf <- cbind(filtereddf, PSMScore = as.double(unfiltereddf$Hyperscore))
    fmessage("Successfully imported Hyperscore column and renamed to PSMScore.")}
  else{stop("The column Hyperscore was not found in the input dataframe.")}

  if ("Assigned.Modifications" %in% existingCols) {
    filtereddf <- cbind(filtereddf, AssignedModifications = as.character(unfiltereddf$Assigned.Modifications))
    fmessage("Successfully imported Assigned Modifications column.")}
  else {stop("The column Assigned.Modifications was not found in the input dataframe.")}

  if ("Total.Glycan.Composition" %in% existingCols) {
    filtereddf <- cbind(filtereddf, TotalGlycanComposition = as.character(unfiltereddf$Total.Glycan.Composition))
    filtereddf$TotalGlycanComposition <- sapply(filtereddf$TotalGlycanComposition, function(x) CleanGlycanNames(x))
    fmessage("Successfully imported Total Glycan Composition column.")}
  else {stop("The column Total.Glycan.Composition was not found in the input dataframe.")}

  if ("Glycan.q.value" %in% existingCols) {
    filtereddf <- cbind(filtereddf, GlycanQValue = as.double(unfiltereddf$Glycan.q.value))
    fmessage("Successfully imported Glycan q Value column.")}
  else if("Confidence.Level" %in% existingCols){
    filtereddf <- cbind(filtereddf, GlycanQValue = unfiltereddf$Confidence.Level)

    filtereddf$GlycanQValue <- dplyr::case_when(filtereddf$GlycanQValue %in% c("Level1", "Level1b") ~ "0",
                                                    filtereddf$GlycanQValue %in% c("Level2") ~ "0.05",
                                                    filtereddf$GlycanQValue %in% c("Level3") ~ "0.1",
                                                    TRUE ~ filtereddf$GlycanQValue)
    filtereddf$GlycanQValue <- as.numeric(filtereddf$GlycanQValue)
    fmessage("Successfully imported Glycan q Value column.")}
  else {stop("The column Glycan.q.value was not found in the input dataframe.")}

  if ("Is.Unique" %in% existingCols) {
    filtereddf <- cbind(filtereddf, IsUnique = as.logical(unfiltereddf$Is.Unique))
    fmessage("Successfully imported Is Unique column.")}
  else {stop("The column Is.Unique was not found in the input dataframe.")}

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

  if ("Gene" %in% existingCols) {
    if("Mapped.Genes" %in% existingCols){
      filtereddf$Genes <- apply(unfiltereddf[, c("Gene", "Mapped.Genes")], 1, function(x) CombineTwoCols(x))
      fmessage("Successfully imported and merged Gene and Mapped Genes columns.")
    }else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Genes columns found, only taking from Gene column")}
  } else {stop("The column Gene was not found in the input dataframe.")}

  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::mutate(.by = .data$UniprotIDs,
                      ProteinLength = GetProteinLength(IDVec = .data$UniprotIDs,
                                                       fastaFile = fastaFile))
    }else{warning("Fasta path does not exist.")}
  }

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

  if ("Protein.Start" %in% existingCols) {
    filtereddf <- cbind(filtereddf, ProteinStart = as.numeric(unfiltereddf$Protein.Start))}
  else {filtereddf$ProteinStart = NA
    warning("The column Is.Unique was not found in the input dataframe.")}

  filtereddf$GlycanType <- apply(filtereddf[,c("AssignedModifications", "TotalGlycanComposition")], 1, function(x) GlycanComptToGlycanType(mod = x[1], glycanComp = x[2]))
  filtereddf <- filtereddf %>%
    dplyr::mutate(GlycanType = sapply(.data$GlycanType, toString))
  fmessage("Successfully added GlycanType column.")

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf, by = "Run")

  filtereddf$Alias <- factor(filtereddf$Alias, levels = unique(annotationdf$Alias))

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}
