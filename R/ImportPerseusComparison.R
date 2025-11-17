#' ImportPerseusComparison
#'
#' This imports a Perseus comparison and returns it as a dataframe. The Conditions
#' specified in Perseus must match the Condition values in the GlycoDiveR annotation
#' file. Perseus formats column names as such: Condition1_Condition2. GlycoDiveR
#' creates identical mapping.
#'
#' @param input The GlycoDiveR inputted data.
#' @param path The path to the Perseus .txt file.
#' @param cleanCCarbamidomethylation Removes \\[57.0215\\] from the MSstats
#' ModifiedPeptide column (default = TRUE)
#'
#' @returns A comparison dataframe.
#' @export
#'
#' @examples \dontrun{ImportPerseusComparison(mydata, C:/PerseusOutput.txt)}
ImportPerseusComparison <- function(input, path, cleanCCarbamidomethylation = TRUE){
  # Read the data using read.delim from text
  columnNames <- names(utils::read.delim(path))
  Perseus_raw <- utils::read.delim(path, skip = 2)
  names(Perseus_raw) <- columnNames

  existingCols <- names(Perseus_raw)

  #Find the comparison
  AllPossibleComparisons <- expand.grid(unique(input$PSMTable$Condition), unique(input$PSMTable$Condition))
  AllPossibleComparisons <- AllPossibleComparisons %>%
    dplyr::mutate(Comparison = paste(.data$Var1, .data$Var2, sep = "_"),
                  FormattedComparison = paste(.data$Var1, .data$Var2, sep = "-"),
                  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ExistsInData = dplyr::case_when(any(grepl(.data$Comparison, existingCols)) ~ "yes",
                                                  TRUE ~ "no")) %>%
    dplyr::ungroup()

  cleanCombinations <- AllPossibleComparisons %>%
    dplyr::filter(.data$ExistsInData == "yes")

  if(nrow(cleanCombinations) > 0){
    fmessage(paste0("Identified the following comparisons: ", paste(cleanCombinations$Comparison, collapse = "; ")))
  }else{
    stop(paste0(
      "Did not identify any comparisons.\nLooked for the following substrings in the column names:\n",
      paste(AllPossibleComparisons$Comparison, collapse = "\n")
    ))
  }

  formattedCombinations <- data.frame(rawColumnName = c(paste0("X.Log.Student.s.T.test.p.value.", cleanCombinations$Comparison),
                                                                         paste0("Student.s.T.test.q.value.", cleanCombinations$Comparison),
                                                                         paste0("Student.s.T.test.Difference.", cleanCombinations$Comparison)),
                                                       CleanColumnName = c(paste0("pvalue;", cleanCombinations$FormattedComparison),
                                                                           paste0("adjpvalue;", cleanCombinations$FormattedComparison),
                                                                           paste0("log2FC;", cleanCombinations$FormattedComparison)))

  #Change colnames
  for(i in seq_len(nrow(formattedCombinations))) {
    Perseus_raw <- Perseus_raw %>%
      dplyr::rename(!! as.character(formattedCombinations[i, 2]) := !!as.character(formattedCombinations[i, 1]))
  }

  #Select the columns of interest and make longer
  Perseus_raw <- Perseus_raw %>%
    dplyr::select("Protein", "Protein.ID", "Modified.Sequence",
                  formattedCombinations$CleanColumnName) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with(c("log2FC", "pvalue", "adjpvalue")),
                        names_to = c(".value", "Label"),
                        names_sep = ";") %>%
    dplyr::mutate(pvalue = 10^(-.data$pvalue))

  fmessage("Finished identifying and extracting the columns of interest.")

  #UniprotIDs####
  if("Protein.ID" %in% existingCols){
    Perseus_raw <- Perseus_raw %>%
      dplyr::rename("UniprotIDs" = "Protein.ID")
    fmessage("Successfully converted UniprotIDs column.")
  }else{stop("The column Protein ID was not found in the input dataframe.")}

  #Proteins####
  if("Protein" %in% existingCols){
    Perseus_raw <- Perseus_raw %>%
      dplyr::mutate(Protein = stringr::str_extract(.data$Protein, "(?<=\\|)[^|_]+(?=_)"))
    fmessage("Successfully converted UniprotIDs column.")
  }else{warning("The column Protein was not found in the input dataframe.")}

  #ModificationID####
  if("Modified.Sequence" %in% existingCols){
    if(cleanCCarbamidomethylation){
      Perseus_raw$Modified.Sequence <- gsub("\\[57\\.0215\\]", "", Perseus_raw$Modified.Sequence)
    }
    ModPepPerseus <- unique(Perseus_raw$Modified.Sequence)
    ModPepData <- unique(input$PSMTable$ModifiedPeptide)

    if(!all(ModPepPerseus %in% ModPepData)){
      warning(paste0(
        "Did not identify all Modified Peptides.\nThe following were not present in the GlycoDiveR data:\n",
        paste(ModPepPerseus[!ModPepPerseus %in% ModPepData], collapse = "\n")
      ))
    }

    nrowBefore <- (nrow(Perseus_raw))
    Perseus_raw <- Perseus_raw %>%
      dplyr::rename("ModifiedPeptide" = "Modified.Sequence")
    Perseus_raw <- input$PTMTable %>%
      dplyr::filter(.data$GlycanType != "NonGlyco") %>%
      dplyr::select("ModifiedPeptide", "ModificationID") %>%
      dplyr::distinct() %>%
      dplyr::right_join(Perseus_raw, by = "ModifiedPeptide")
    nrowAfter <- (nrow(Perseus_raw))
    if(nrowBefore != nrowAfter){
      warning(paste0("Columns before ModificationID joining: ", nrowBefore, ". and after joining: ", nrowAfter))
    }else{
      fmessage("Successfully added ModificationID column.")
    }
  }else{warning("The column Modified Sequence was not found in the input dataframe.")}

  Perseus_raw <- Perseus_raw %>%
    dplyr::select("UniprotIDs", "Proteins" = "Protein", "ModificationID", "ModifiedPeptide",
                  "Label", "log2FC", "pvalue", "adjpvalue")
  return(Perseus_raw)
}

