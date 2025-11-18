#' GetAnnotationTemplate
#'
#' @param path The path the the folder that holds the subforlder with the PSM files
#' @param tool The search engine used. Options: MSFragger
#'
#' @returns Exported annotation file
#' @export
#'
#' @examples \dontrun{GetAnnotationTemplate("Z:/Subfolder", tool = "MSFragger")}
#'
GetAnnotationTemplate <- function(path, tool){
  if(file.exists(paste0(path,"annotation.csv")) || file.exists(paste0(path,"/annotation.csv"))){
    stop("Annotation file already exists. Please remove or edit existing file.")
  }

  if(tool == "MSFragger"){
    fileList <- list.files(path, recursive = TRUE)
    fileList <- fileList[grepl("psm.tsv", fileList)]
    unfiltereddf <- data.frame()

    if(length(fileList) == 0){stop("No files found.")}
    for(file in fileList){
      print(paste0(path, "/", file))
      temptable <- data.table::fread(paste0(path, "/", file), sep = "\t", check.names = TRUE, fill = TRUE)
      unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)
    }
    unfiltereddf$Run <- sapply(unfiltereddf$Spectrum.File, function(x) strsplit(x, "\\", fixed = T)[[1]][length(strsplit(x, "\\", fixed = T)[[1]])-1])

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
  }else if(tool == "Byonic"){
    fileList <- list.files(path, recursive = TRUE)
    fileList <- fileList[grepl("Byonic.xlsx", fileList)]
    unfilteredvec <- c()

    if(length(fileList) == 0){stop("No files found.")}
    for(file in fileList){
      print(paste0(path, "/", file))
      unfilteredvec <- c(unfilteredvec,
                         gsub("3) Parameter file:  |\\\\objs\\\\params.prf", "", as.character(readxl::read_xlsx(paste0(path, "/", file), sheet = 1, range = "B3:B4"))))
    }
    tempdf <- data.frame(Run = unique(unfilteredvec),
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
