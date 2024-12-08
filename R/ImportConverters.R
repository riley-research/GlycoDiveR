source(here::here("R/GlycoDiveRUtils.R"))

MSFraggerConverter <- function(unfiltereddf, annotationdf, fastaPath){
  message("INFO : NOW STARTING IMPORT")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  if("Run" %in% existingCols){filtereddf <- cbind(filtereddf, Run = as.character(unfiltereddf$Run))
  }else{stop("The column Run was not found in the input dataframe.")}

  if("Modified.Peptide" %in% existingCols){filtereddf <- cbind(filtereddf, ModifiedPeptide = as.character(unfiltereddf$Run))}
  else(stop("The column Modified.Peptide was not found in the input dataframe."))

  if("Intensity" %in% existingCols) {filtereddf <- cbind(filtereddf, Intensity = as.character(unfiltereddf$Run))}
  else{filtereddf <- cbind(filtereddf, Intensity = as.numeric(NA))
       warning("Intensity column not found. Filled with NA.")}

  #if("Charge" %in% existingCols) {filtereddf <- cbind(filtereddf, Charge = as.numeric(unfiltereddf$Charge))}
  #else{stop("The column Charge was not found in the input dataframe.")}

  #if("Retention" %in% existingCols) {filtereddf <- cbind(filtereddf, RetentionTime = as.numeric(unfiltereddf$Retention))}
  #else{stop("The column Retention was not found in the input dataframe.")}

  #if("Observed.M.Z" %in% existingCols) {filtereddf <- cbind(filtereddf, ObservedMZ = as.double(unfiltereddf$Observed.M.Z))}
  #else{stop("The column Observed.M.Z. was not found in the input dataframe.")}

  #if ("Delta.Mass" %in% existingCols) {filtereddf <- cbind(filtereddf, DeltaMass = as.double(unfiltereddf$Delta.Mass))}
  #else {stop("The column Delta.Mass was not found in the input dataframe.")}

  if ("Hyperscore" %in% existingCols) {filtereddf <- cbind(filtereddf, PSMScore = as.double(unfiltereddf$Hyperscore))}
  else{stop("The column Hyperscore was not found in the input dataframe.")}

  #if ("Nextscore" %in% existingCols) {filtereddf <- cbind(filtereddf, Nextscore = as.double(unfiltereddf$Nextscore))}
  #else {stop("The column Nextscore was not found in the input dataframe.")}

  if ("Assigned.Modifications" %in% existingCols) {filtereddf <- cbind(filtereddf, AssignedModifications = as.character(unfiltereddf$Assigned.Modifications))}
  else {stop("The column Assigned.Modifications was not found in the input dataframe.")}

  if ("Total.Glycan.Composition" %in% existingCols) {filtereddf <- cbind(filtereddf, TotalGlycanComposition = as.character(unfiltereddf$Total.Glycan.Composition))}
  else {stop("The column Total.Glycan.Composition was not found in the input dataframe.")}

  #if ("Glycan.Score" %in% existingCols) {filtereddf <- cbind(filtereddf, GlycanScore = as.numeric(unfiltereddf$Glycan.Score))}
  #else {stop("The column Glycan.Score was not found in the input dataframe.")}

  if ("Glycan.q.value" %in% existingCols) {filtereddf <- cbind(filtereddf, GlycanQValue = as.double(unfiltereddf$Glycan.q.value))}
  else {stop("The column Glycan.q.value was not found in the input dataframe.")}

  if ("Is.Unique" %in% existingCols) {filtereddf <- cbind(filtereddf, IsUnique = as.logical(unfiltereddf$Is.Unique))}
  else {stop("The column Is.Unique was not found in the input dataframe.")}

  if ("Protein.ID" %in% existingCols) {
    if("Mapped.Proteins" %in% existingCols){
      filtereddf$UniprotIDs <- apply(unfiltereddf[, c("Protein.ID", "Mapped.Proteins")], 1, function(x) CombineProtCols(x))
    }else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Proteins columns found, only taking from Protein ID column")}
    } else {stop("The column Protein.ID was not found in the input dataframe.")}

  if ("Gene" %in% existingCols) {
    if("Mapped.Genes" %in% existingCols){
      filtereddf$Genes <- apply(unfiltereddf[, c("Gene", "Mapped.Genes")], 1, function(x) CombineTwoCols(x))
    }else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Genes columns found, only taking from Gene column")}
  } else {stop("The column Gene was not found in the input dataframe.")}

  #if ("Protein.Description" %in% existingCols) {filtereddf <- cbind(filtereddf, ProteinDescription = as.character(unfiltereddf$Protein.Description))}
  #else {stop("The column Protein.Description was not found in the input dataframe.")}

  #if ("Mapped.Genes" %in% existingCols) {filtereddf <- cbind(filtereddf, MappedGenes = as.character(unfiltereddf$Mapped.Genes))}
  #else {stop("The column Mapped.Genes was not found in the input dataframe.")}

  if("UniprotIDs" %in% names(filtereddf)){
    if(file.exists(fastaPath)){
      fastaFile <- seqinr::read.fasta(file = fastaPath)
      filtereddf <- filtereddf %>%
        dplyr::ungroup() %>%
        dplyr::group_by(UniprotIDs) %>%
        dplyr::mutate(ProteinLength = GetProteinLength(IDVec = UniprotIDs, fastaFile = fastaFile)) %>%
        dplyr::ungroup()
    }else{warning("Fasta path did not exist.")}
  }

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf)

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}
