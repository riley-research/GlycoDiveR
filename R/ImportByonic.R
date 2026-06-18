#' ImportByonic
#'
#' This is the importer for Byonic data. It needs to untouched Byonic output, which
#' is a Excel file with at least a Summary and a Spectra tab. This script will look
#' for all files in the folder and all subfolders using the name tag "Byonic.xslx".
#'
#' @param path The folder path, which will be used to find all "Byonic.xsls" files
#' in this folder and all subfolders.
#' @param annotation The annotation dataframe. Please generate a template using
#' GetAnnotationTemplate("path", tool = "Byonic")
#' @param fastaPath The path to the fasta file that was used.
#' @param peptideScoreCutoff The score cutoff from the "Score" column. The cutoff
#' is the lower limit.
#' @param glycanScoreCutoff The score cutoff for the Log Prob column. The Log Prob
#' column is converted using 10^(-|Log Prob|). The cutoff is the
#' upper limit.
#' @param deltaModCutoff The score cutoff for the Delta Mod score. The score is the
#' lower limit.
#' @param scrape set TRUE/FALSE to use scraping of Uniprot data.
#' @param removeReverse Remove reverse hits as noted with the ">Reverse " tag.
#' @param cutoffFilter Filter for the glycan score, peptide score, N Sequon, and
#' confidence levels.
#' @param dropNoQuant Remove or keep PSMs without a quantitative value or have a quant
#' of 0.
#' @param minPeptideCoverage Works together with the thresholdMode argument.
#' If thresholdMode = "total", provide a single numeric value. Values between 0
#' and 1 are interpreted as proportions (e.g., 0.1 = 10%), while values >= 1 are
#' interpreted as absolute counts (e.g., 2 = at least two peptides identified).
#' If thresholdMode = "group", provide a numeric vector of length two (e.g., c(0.1, 3)).
#' The first value is interpreted in the same way as for "total" (proportion or
#' absolute count). The second value specifies the minimum number of groups in
#' which this threshold must be met.
#' @param thresholdMode Character string specifying the filtering mode: "total"
#' applies a global threshold across all data, while "group" identifies the threshold
#' within groups.
#' @param GlyToucan set to TRUE will connect to the GlyCosmos API to retrieve
#' GlyToucan identifiers.
#'
#' @returns Formatted GlycoDiveR data file.
#' @export
#'
#' @examples \dontrun{ImportByonic(path = "C:/Byonic",
#' annotation = "C:/annotation.csv",
#' fastaPath = "C:/fastafile.fasta",
#' peptideScoreCutoff = 0,
#' glycanScoreCutoff = 1)}
ImportByonic <- function(path, annotation, fastaPath, peptideScoreCutoff = 0, glycanScoreCutoff = 0.01,
                         deltaModCutoff = 1, scrape = TRUE, removeReverse = TRUE, cutoffFilter = TRUE,
                         dropNoQuant = FALSE, minPeptideCoverage = FALSE,
                         thresholdMode = c("group", "total"), GlyToucan = TRUE){
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

  unfiltereddf$ID <-   seq(1:nrow(unfiltereddf))

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
        ~ if (..1 == "Glyco") as.character(ComputeGlycanMass(..2)) else ..3))

  #Add variable mods to ModificationDatabase
  modToAdd <- modification_df %>%
    dplyr::filter(.data$ModificationType != "Glyco") %>%
    dplyr::select("FullName" = "Rule",
                  "ModificationMass" = "Mass") %>%
    dplyr::mutate(ModificationMass = sprintf("%.3f", as.numeric(.data$ModificationMass)))

  .modEnv$ModificationDatabase <- dplyr::bind_rows(.modEnv$ModificationDatabase, modToAdd)

  filtereddf <- ByonicConverter(unfiltereddf, annotationdf, fastaPath,
                                modification_df, scrape, GlyToucan)

  if(cutoffFilter){
    filtereddf <- FilterPSMTable(filtereddf,
                                 peptideScoreCutoff,
                                 glycanScoreCutoff,
                                 filterForNoNSequon = FALSE,
                                 confidenceLevels = FALSE,
                                 deltaModCutoff,
                                 searchEngine = "Byonic",
                                 SScoreCutoff = FALSE,
                                 structureConfidenceCutoff = FALSE,
                                 AScoreCutoff = FALSE)}

  if(dropNoQuant){filtereddf <- filtereddf %>% dplyr::filter(!is.na(.data$Intensity) & .data$Intensity != 0)}

  if(!identical(minPeptideCoverage, FALSE)){filtereddf <- FilterForMinPeptides(filtereddf, minPeptideCoverage, thresholdMode)}

  PTMdf <- PSMToPTMTable(filtereddf)

  data <- list(PSMTable = filtereddf,
               rawPSMTable = unfiltereddf,
               PTMTable = PTMdf,
               annotation = annotationdf,
               searchEngine = "Byonic",
               peptideScoreCutoff = peptideScoreCutoff,
               glycanScoreCutoff = glycanScoreCutoff,
               filterForNoNSequon = FALSE,
               confidenceLevels = FALSE,
               deltaModCutoff = deltaModCutoff,
               SScoreCutoff = FALSE,
               structureConfidenceCutoff = FALSE,
               AScoreCutoff = FALSE)

  return(data)
}
