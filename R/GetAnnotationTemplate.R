GetAnnotationTemplate <- function(path, tool){
  if(tool == "MSFragger"){
    fileList <- list.files(path, recursive = TRUE)
    fileList <- fileList[grepl("psm.tsv", fileList)]
    unfiltereddf <- data.frame()

    for(file in fileList){
      temptable <- utils::read.table(paste0(path, "/", file), sep = "\t", header = T)
      unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
    }

    lengthList <- length(strsplit(unfiltereddf$Spectrum.File[1], "\\", fixed = T)[[1]])
    unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][lengthList-1])

    tempdf <- data.frame(Run = unique(unfiltereddf$Run),
                         Condition = NA,
                         Alias = NA,
                         BioReplicate = NA,
                         TechReplicate = NA)

    if(grepl("/$", path)){
      utils::write.csv(tempdf, paste0(path,"annotation.csv"), row.names=FALSE)
      message("INFO: ", "Annotation dataframe exported to: ", paste0(path,"annotation.csv"))
    }else{
      utils::write.csv(tempdf, paste0(path,"/annotation.csv"), row.names=FALSE)
      message("INFO: ", "Annotation dataframe exported to: ", paste0(path,"/annotation.csv"))}

  }
  else{
    warning("tool does not match any recognized tools.")
  }
}
