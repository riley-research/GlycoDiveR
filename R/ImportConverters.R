MSFraggerConverter <- function(unfiltereddf){
  filtereddf <- unfiltereddf[,c("Run", "Modified.Peptide", "Intensity", "Charge",
                                "Retention", "Observed.M.Z", "Delta.Mass", "Hyperscore",
                                "Nextscore", "Assigned.Modifications", "Total.Glycan.Composition",
                                "Glycan.Score", "Glycan.q.value", "Is.Unique",
                                "Protein.ID", "Gene", "Protein.Description",
                                "Mapped.Genes")]

  colnames(filtereddf) <- c("Run", "ModifiedPeptide", "Intensity", "Charge",
                            "RetentionTime", "ObservedMZ", "DeltaMass", "Score",
                            "DeltaScore", "Modifications", "GlycanComposition",
                            "GlycanScore", "GlycanQValue", "IsUnique", "UniprotID",
                            "Gene", "ProteinDescription", "Proteins")
  return(filtereddf)
}
