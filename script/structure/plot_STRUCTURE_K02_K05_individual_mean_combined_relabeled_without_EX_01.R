args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}

repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source_script <- file.path(
  repo_root,
  "script",
  "structure",
  "plot_STRUCTURE_allK_individual.R"
)

script_text <- paste(readLines(source_script, warn = FALSE), collapse = "\n")

# Retain the original functions and K02-K05 plotting settings, skip unrelated
# single-K outputs, and remove EX_01 only from the combined plot data.
script_text <- sub(
  "repo_root <- normalizePath(file.path(script_dir, \"..\"))",
  "repo_root <- get(\"repo_root\", envir = .GlobalEnv)",
  script_text,
  fixed = TRUE
)
script_text <- sub("for (k in 1:5)", "for (k in integer(0))", script_text, fixed = TRUE)
script_text <- sub(
  "  return(out)\n\n  \n})",
  "  return(out %>% filter(population_label != \"EX_01\"))\n\n  \n})",
  script_text,
  fixed = TRUE
)
script_text <- strsplit(script_text, "\n########\n", fixed = TRUE)[[1]][1]

Sys.setenv(
  PHYSARIA_WEN_LABEL_FILE = file.path(
    repo_root,
    "data",
    "metadata",
    "Physaria_globosa_population_counts_from_individuals_relabeled.csv"
  ),
  PHYSARIA_OUTPUT_SUFFIX = "_relabeled_without_EX_01"
)

eval(parse(text = script_text), envir = new.env(parent = globalenv()))
