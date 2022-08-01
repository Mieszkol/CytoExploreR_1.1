#' Import, pregate and extract population of interest from your samples to make the files smaller.
#' By default it imports all of the ".fcs" and ".FCS" files from specified directory.
#' The function skips the compensation and transformation - pregating can be only done using FSC and SSC channels.
#'
#' @param path path points to the location of the .fcs files to read in. Preferably
#'   the name of folder in current working directory
#'
#' @param gatingTemplate name of a gatingTemplate csv file to be used for gate saving.
#' Set to "gatingTemplate_pregating.csv" by default
#'
#' @param output name of the folder to which the written FCS files should be
#'   saved, set to NULL by default to save the files to the current working
#'   directory
#'
#' @param gates sets gates to use for exporting the population of interest. Set to c(
#' cyto_gate_draw(gs, parent = "root", alias = "Cells", channels = c("FSC-A", "SSC-A")),
#' cyto_gate_draw(gs, parent = "Cells", alias = "Single Cells", channels = c("FSC-A", "FSC-H"))
#' ) by default
#'
#' @param stats defines the statistics to be exported to control csv sheet. Set to stats = c(
#'cyto_stats_compute(gs, parent = "root", alias = "Cells", stat = "freq"),
#'cyto_stats_compute(gs, parent = "Cells", alias = "Single Cells", stat = "freq")
#') by default
#'
#' @param control name of the file compiling the frequency of each gate for all of the samples. Set to NULL by default - no file created
#'
#'
#' @return new fcs files extracted from subpopulation of interest
#'
#' @importFrom utils write.csv
#' @importFrom dplyr bind_rows bind_cols select %>%
#' @importFrom tools file_ext
#'
#' @examples
#'
#' #Load CytoExploreR
#' library(CytoExploreR)
#'
#' # Get path to Activation .fcs files in CytoExploreRData
#' datadir <- system.file("extdata", package = "CytoExploreRData")
#' path <- paste0(datadir, "/Activation")
#'
#' # Load and preprocess the fcs files
#' cyto_fcs_pregating(
#' path,
#' output = "pregated",
#' control = "pregated"
#' )
#'
#' @author Mieszko Lachota
#'
#' @export

cyto_fcs_pregating <- function(path = ".",
                               gatingTemplate = "gatingTemplate_pregating.csv",
                               output = NULL,
                               gates = c(
                                 cyto_gate_draw(
                                   gs,
                                   parent = "root",
                                   alias = "Cells",
                                   channels = c("FSC-A", "SSC-A")
                                 ),
                                 cyto_gate_draw(
                                   gs,
                                   parent = "Cells",
                                   alias = "Single Cells",
                                   channels = c("FSC-A", "FSC-H")
                                 )
                               ),
                               stats = c(
                                 cyto_stats_compute(
                                   gs,
                                   parent = "root",
                                   alias = "Cells",
                                   stat = "freq"
                                 ),
                                 cyto_stats_compute(
                                   gs,
                                   parent = "Cells",
                                   alias = "Single Cells",
                                   stat = "freq"
                                 )
                               ),
                               control = NULL,
                               output_gate = "Single Cells") {
  # INITIAL CHECKS
  wd <- getwd()
  setwd(path)
  if (is.null(output)) {
    stop("Please specify output folder")
  }
  if (file.exists(output)) {
    stop("There are files in the specified output folder")
  }
  if (file.exists(gatingTemplate)) {
    stop("There is a gatingTemplate with confiltcting name")
  }
  if (file.exists(paste(control, ".csv", sep = ""))) {
    stop("There is a control datasheet with confilcting name")
  }
  # CONVERT FUNCTIONS TO CHAR STRINGS
  gates_char <- as.character(substitute(gates))
  gates_char <- gates_char[2:length(gates_char)]
  stats_char <- as.character(substitute(stats))
  stats_char <- stats_char[2:length(stats_char)]
  
  # CREATES OUTPUT DIRECTORY
  dir.create(output, recursive = TRUE)
  # FILE PATHS
  file_paths <- list.files(full.names = FALSE)
  # FILE NAMES
  file_names <- basename(file_paths)
  # FILE EXTENSIONS
  file_ext <- file_ext(file_names)
  # FCS FILES
  fcs_files <- file_paths[file_ext %in% c("fcs", "FCS")]
  # CHECKS
  if (is.null(fcs_files)) {
    stop("There are no fcs files in the specified folder")
  }
  # FILE NAMES WITHOUT EXTENSIONS
  no_ext_names <- file_path_sans_ext(fcs_files)
  # NUMBER OF FCS FILES
  fcs_number <- NULL
  dir.names <- NULL
  preprocessing.stats <- NULL
  for (i in 1:length(no_ext_names)) {
    fcs_number <- c(fcs_number, i)
    # CHECKS
    if (file.exists(paste(fcs_number[i], no_ext_names[i], sep = "_"))) {
      stop("There are folders/files with confilcting name")
    }
    # MAKES A SEPARATE FOLDER FOR EVERY FILE
    dir.create(paste(fcs_number[i], no_ext_names[i], sep = "_"), recursive = TRUE)
    # COPIES EVERY FILE INTO ITS FOLDER
    file.copy(fcs_files[i], paste(fcs_number[i], no_ext_names[i], sep = "_"))
    # LOADS THE OTHER FCS FILES ONE-BY-ONE
    gs <-
      cyto_setup(
        paste(fcs_number[i], no_ext_names[i], sep = "_"),
        gatingTemplate = gatingTemplate,
        clean = FALSE,
        markers = FALSE,
        details = FALSE
      )
    # SETS UP GATES
    if (!(length(cyto_gatingTemplate_read(gatingTemplate)$alias) > 0)) {
      for (l in 1:length(gates_char)) {
        eval(parse(text = gates_char[l]))
      }
    }
    # APPLIES GATINGTEMPLATE FROM THE FIRST FILE
    if (length(cyto_gatingTemplate_read(gatingTemplate)$alias) > 0) {
      cyto_gatingTemplate_apply(gs, gatingTemplate)
    }
    # SAVES THE EXTRACTED POPULATION
    cyto_save(gs, parent = output_gate, save_as = output)
    # SAVES THE FREQ STATS
    preprocessing.stats <- rbind(for (j in 1:length(stats_char)) {
      eval(parse(text = stats_char[j]))
    })
    # DELETES CREATED DIRECTORY AND FILE
    unlink(paste(fcs_number[i], no_ext_names[i], sep = "_"), recursive = TRUE)
  }
  # EXPORTS FREQ STATS IF THE NAME IS PROVIDED
  if (!is.null(control)) {
    write.csv(preprocessing.stats,
              paste(control, ".csv", sep = ""),
              row.names = FALSE)
  }
 
  # SETS WORKING DIRECTORY BACK TO THE ORIGINAL
  setwd(wd)
}