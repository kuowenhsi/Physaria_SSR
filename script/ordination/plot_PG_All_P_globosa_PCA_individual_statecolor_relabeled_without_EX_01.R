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
  "ordination",
  "plot_PG_All_pca_dapc.R"
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

# Load the original functions and settings without running its output block.
script_text <- paste(readLines(source_script, warn = FALSE), collapse = "\n")
script_text <- strsplit(script_text, "\nparsed <- read_pg_all", fixed = TRUE)[[1]][1]
script_text <- sub(
  "repo_root <- normalizePath(file.path(script_dir, \"..\", \"..\"))",
  "repo_root <- get(\"repo_root\", envir = .GlobalEnv)",
  script_text,
  fixed = TRUE
)

plot_env <- new.env(parent = globalenv())
eval(parse(text = script_text), envir = plot_env)

with(plot_env, {
  parsed <- read_pg_all(input_file)
  genind_globosa <- build_genind(parsed)
  pca_globosa <- run_pca(genind_globosa, parsed$metadata)

  pca_axis_labels <- c(
    sprintf("PC1 (%.1f%%)", 100 * pca_globosa$explained[1]),
    sprintf("PC2 (%.1f%%)", 100 * pca_globosa$explained[2])
  )

  plot_scores <- pca_globosa$ind_scores %>%
    filter(population_label != "EX_01")
  plot_pop_scores <- pca_globosa$pop_scores %>%
    filter(population_label != "EX_01")

  p <- plot_individual_ordination_with_population_overlay_statecolor(
    plot_scores,
    plot_pop_scores,
    method = "PCA",
    axis_labels = pca_axis_labels,
    subtitle = NULL
  ) +
    # Invisible full-data layer preserves the original x/y scaling.
    geom_blank(
      data = pca_globosa$ind_scores,
      aes(x = Axis1, y = Axis2),
      inherit.aes = FALSE
    )

  save_plot(p, "PG_All_P_globosa_PCA_individual_statecolor.png")
})
