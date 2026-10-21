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
  "maps",
  "plot_P_globosa_water_overlay_scatterpie_K02_K08_hfp.R"
)

script_lines <- readLines(source_script, warn = FALSE)
script_text <- paste(script_lines, collapse = "\n")

# Reuse the original map workflow in memory while limiting this copy to K = 3,
# using the relabeled STRUCTURE input, and excluding EX_01 from every point layer.
script_text <- sub(
  "repo_root <- normalizePath(file.path(script_dir, \"..\"))",
  "repo_root <- get(\"repo_root\", envir = .GlobalEnv)",
  script_text,
  fixed = TRUE
)
script_text <- sub("k_values <- 2:8", "k_values <- 3L", script_text, fixed = TRUE)
script_text <- sub(
  "line_tbl <- point_xy_base %>%",
  "line_tbl <- point_xy_base %>%\n  filter(point_id != \"EX_01\") %>%",
  script_text,
  fixed = TRUE
)
script_text <- sub(
  "      )\n    )\n\n  if (any(!stats::complete.cases(point_xy[, cluster_cols, drop = FALSE])))",
  "      )\n    ) %>%\n    filter(point_id != \"EX_01\")\n\n  if (any(!stats::complete.cases(point_xy[, cluster_cols, drop = FALSE])))",
  script_text,
  fixed = TRUE
)
script_text <- sub(
  "structure_file <- file.path(structure_dir, suffix_filename(sprintf(\"STRUCTURE_K%02d_individual_mean.csv\", k)))",
  "structure_file <- file.path(structure_dir, sprintf(\"STRUCTURE_K%02d_individual_mean_relabeled.csv\", k))",
  script_text,
  fixed = TRUE
)
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
