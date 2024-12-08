CombineTwoCols <- function(row){
  if (is.na(row[2]) | row[2] == ""){return(as.character(row[1]))}
  else{return(as.character(paste(row[1], row[2], sep = ",")))}}

CombineProtCols <- function(row){
  if (is.na(row[2]) | row[2] == ""){
    return(row[1])
  }else{
    addcol <- unlist(regmatches(row[2], gregexpr("(?<=\\|)[^|]+(?=\\|)", row[2], perl = TRUE)))
    addcol <- addcol[!grepl(",|_", addcol)]
    return(as.character(paste(row[1], paste(addcol, collapse=","), sep = ",")))
  }
}

GetProteinLength <- function(IDVec, fastaFile){
  inputVector <- strsplit(IDVec[1], ",")[[1]]
  uniprotID <- inputVector[1]
  IDhit <- fastaFile[grepl(uniprotID, names(fastaFile)) & !grepl("rev", names(fastaFile))]
  protLength <- nrow(as.data.frame(IDhit[[1]]))
  return(as.integer(protLength))
}
