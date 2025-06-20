GetAnnotationTemplate <- function(path, tool){
  if(tool == "MSFragger"){
    fileList <- list.files(path, recursive = TRUE)
    fileList <- fileList[grepl("psm.tsv", fileList)]
    unfiltereddf <- data.frame()

    if(length(fileList) == 0){stop("No files found.")}
    for(file in fileList){
      print(paste0(path, "/", file))
      temptable <- data.table::fread(paste0(path, "/", file), sep = "\t", check.names = TRUE, fill = TRUE)
      #temptable <- utils::read.table(paste0(path, "/", file), sep = "\t", col.names = T, quote = "")
      unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
    }
    View(unfiltereddf)
    lengthList <- length(strsplit(unfiltereddf$Spectrum.File[1], "\\", fixed = T)[[1]])
    unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][lengthList-1])

    tempdf <- data.frame(Run = unique(unfiltereddf$Run),
                         Condition = NA,
                         Alias = NA,
                         BioReplicate = NA,
                         TechReplicate = NA)

    if(grepl("/$", path)){
      utils::write.csv(tempdf, paste0(path,"annotation.csv"), row.names=FALSE)
      fmessage(paste0("Annotation dataframe exported to: ", paste0(path,"annotation.csv")))
    }else{
      utils::write.csv(tempdf, paste0(path,"/annotation.csv"), row.names=FALSE)
      fmessage(paste0("Annotation dataframe exported to: ", paste0(path,"/annotation.csv")))}

  }
  else{
    warning("tool does not match any recognized tools.")
  }
}
