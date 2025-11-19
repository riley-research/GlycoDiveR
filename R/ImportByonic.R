#' ImportByonic
#'
#' This is the importer for Byonic data. It needs to untouched Byonic output, which
#' is a Excel file with at least a Summary and a Spectra tab. This script will look
#' for all files in the folder and all subfolders using the name tag "Byonic.xslx".
#'
#' @param path The folder path, which will be used to find all "Byonic.xlsl" files
#' in this folder and all subfolders.
#' @param annotation The annotation dataframe. Please generate a template using
#' GetAnnotationTemplate("path", tool = "Byonic")
#' @param fastaPath The path to the fasta file that was used.
#' @param peptideScoreCutoff The score cutoff from the "Score" column. The cutoff
#' is the lower limit.
#' @param glycanScoreCutoff The score cutoff in the Log Prob column. The score is the
#' upper limit.
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param removeReverse Remove reverse hits as noted with the ">Reverse " tag.
#'
#' @returns Formatted GlycoDiveR data file.
#' @export
#'
#' @examples \dontrun{ImportByonic(path = "C:/Byonic",
#' annotation = "C:/annotation.csv",
#' fastaPath = "C:/fastafile.fasta",
#' peptideScoreCutoff = 0,
#' glycanScoreCutoff = 1)}
ImportByonic <- function(path, annotation, fastaPath, peptideScoreCutoff, glycanScoreCutoff,
                              scrape = FALSE, removeReverse = TRUE){
  unfiltereddf <- data.frame()
  modification_df <- data.frame()
  annotationdf <- utils::read.csv(annotation)
  CheckAnnotation(annotationdf)

  fileList <- list.files(path, recursive = TRUE)
  fileList <- fileList[grepl("Byonic.xlsx", fileList)]

  #Bind all data together
  if(length(fileList) == 0){stop("No files found.")}
  for(file in fileList){
    Run <- gsub("3) Parameter file:  |\\\\objs\\\\params.prf", "", as.character(readxl::read_xlsx(paste0(path, "/", file), sheet = 1, range = "B3:B4")))
    temptable <- readxl::read_xlsx(paste0(path, "/", file), 2)
    temptable$Run <- Run
    unfiltereddf <- plyr::rbind.fill(unfiltereddf, temptable)

    temptable2 <- readxl::read_xlsx(paste0(path, "/", file), 1, skip = 14) %>%
      dplyr::filter(grepl("@", .data$Rule) ) %>% #& !grepl("Glycan", .data$Rule)
      dplyr::select("Rule")

    modification_df <- plyr::rbind.fill(modification_df, temptable2)
  }

  if(removeReverse){
    unfiltereddf <- unfiltereddf %>%
      dplyr::filter(!grepl(">Reverse ", .data$`Protein Name`))
    fmessage("Removed >Reverse proteins")
  }

  #Clean the modification df
  modification_df <- modification_df %>%
    dplyr::mutate(
      ModificationType = ifelse(grepl("Glycan", .data$Rule), "Glyco", "NonGlyco"),
      Rule = ifelse(
        grepl("Glycan", .data$Rule),
        stringr::str_extract(.data$Rule, "^[^@]+"),  # everything before the first " @ "
        .data$Rule
      ),
      Rule = ifelse(.data$ModificationType == "Glyco",
                    CleanGlycanNames(.data$Rule),
                    .data$Rule))

  modification_df <- modification_df %>%
    dplyr::distinct() %>%
    dplyr::mutate(Mass = stringr::str_extract(.data$Rule, "[+-]?\\d+\\.\\d+"),
      Mass = sub("^\\+", "", .data$Mass),
      Mass = purrr::pmap_chr(
        list(.data$ModificationType,
             .data$Rule,
             .data$Mass),
        ~ if (..1 == "Glyco") ComputeGlycanMass(..2) else ..3))

  #Add variable mods to ModificationDatabase
  modToAdd <- modification_df %>%
    dplyr::filter(.data$ModificationType != "Glyco") %>%
    dplyr::select("FullName" = "Rule",
                  "ModificationMass" = "Mass") %>%
    dplyr::mutate(ModificationMass = sprintf("%.3f", as.numeric(.data$ModificationMass)))

  .modEnv$ModificationDatabase <- dplyr::bind_rows(.modEnv$ModificationDatabase, modToAdd)

  filtereddf <- ByonicConverter(unfiltereddf, annotationdf, fastaPath,
                  modification_df, scrape)

  PTMdf <- PSMToPTMTable(filtereddf)

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "Byonic",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = FALSE)

  return(data)
}
