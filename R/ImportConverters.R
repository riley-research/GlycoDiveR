source(here::here("R/GlycoDiveRUtils.R"))

MSFraggerConverter <- function(unfiltereddf, annotationdf, fastaPath){
  message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Now starting import.\033[0m")
  filtereddf <- data.frame(ID = seq(1:nrow(unfiltereddf)))
  existingCols <- unique(names(unfiltereddf))

  if("Run" %in% existingCols){
    filtereddf <- cbind(filtereddf, Run = as.character(unfiltereddf$Run))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Run column.\033[0m")
  }
  else{stop("The column Run was not found in the input dataframe.")}

  if("Modified.Peptide" %in% existingCols){
    filtereddf <- cbind(filtereddf, ModifiedPeptide = as.character(unfiltereddf$Modified.Peptide))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Modified Peptide column.\033[0m")}
  else(stop("The column Modified.Peptide was not found in the input dataframe."))

  if("Intensity" %in% existingCols) {
    filtereddf <- cbind(filtereddf, Intensity = as.character(unfiltereddf$Run))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Intensity column.\033[0m")}
  else{filtereddf <- cbind(filtereddf, Intensity = as.numeric(NA))
       warning("Intensity column not found. Filled with NA.")}

  if ("Hyperscore" %in% existingCols) {
    filtereddf <- cbind(filtereddf, PSMScore = as.double(unfiltereddf$Hyperscore))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Hyperscore column and renamed to PSMScore.\033[0m")}
  else{stop("The column Hyperscore was not found in the input dataframe.")}

  if ("Assigned.Modifications" %in% existingCols) {
    filtereddf <- cbind(filtereddf, AssignedModifications = as.character(unfiltereddf$Assigned.Modifications))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Assigned Modifications column.\033[0m")}
  else {stop("The column Assigned.Modifications was not found in the input dataframe.")}

  if ("Total.Glycan.Composition" %in% existingCols) {
    filtereddf <- cbind(filtereddf, TotalGlycanComposition = as.character(unfiltereddf$Total.Glycan.Composition))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Total Glycan Composition column.\033[0m")}
  else {stop("The column Total.Glycan.Composition was not found in the input dataframe.")}

  if ("Glycan.q.value" %in% existingCols) {
    filtereddf <- cbind(filtereddf, GlycanQValue = as.double(unfiltereddf$Glycan.q.value))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Glycan q Value column.\033[0m")}
  else {stop("The column Glycan.q.value was not found in the input dataframe.")}

  if ("Is.Unique" %in% existingCols) {
    filtereddf <- cbind(filtereddf, IsUnique = as.logical(unfiltereddf$Is.Unique))
    message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported Is Unique column.\033[0m")}
  else {stop("The column Is.Unique was not found in the input dataframe.")}

  if ("Protein.ID" %in% existingCols) {
    if ("Mapped.Proteins" %in% existingCols) {
      filtereddf$UniprotIDs <- apply(unfiltereddf[, c("Protein.ID", "Mapped.Proteins")], 1, function(x) CombineProtCols(x))
      message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported and merged Protein ID and Mapped protein columns.\033[0m")
    } else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Proteins column found, only taking from Protein ID column.")
    }
  } else {
    stop("The column Protein.ID was not found in the input dataframe.")}

  if ("Gene" %in% existingCols) {
    if("Mapped.Genes" %in% existingCols){
      filtereddf$Genes <- apply(unfiltereddf[, c("Gene", "Mapped.Genes")], 1, function(x) CombineTwoCols(x))
      message("\033[30m[", base::substr(Sys.time(), 1, 16), "] INFO: Successfully imported and merged Gene and Mapped Genes columns.\033[0m")
    }else {
      filtereddf <- cbind(filtereddf, UniprotIDs = as.character(unfiltereddf$Protein.ID))
      warning("No Mapped.Genes columns found, only taking from Gene column")}
  } else {stop("The column Gene was not found in the input dataframe.")}

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

  if ("Protein.Start" %in% existingCols) {
    filtereddf <- cbind(filtereddf, ProteinStart = as.numeric(unfiltereddf$Protein.Start))}
  else {filtereddf$ProteinStart = NA
    warning("The column Is.Unique was not found in the input dataframe.")}

  filtereddf <- filtereddf %>%
    dplyr::left_join(annotationdf, by = "Run")

  if(anyNA(filtereddf[c("Condition", "BioReplicate", "TechReplicate", "Alias")])){
    warning("NA detected. Please verify annotation dataframe!")
  }

  return(filtereddf)
}
