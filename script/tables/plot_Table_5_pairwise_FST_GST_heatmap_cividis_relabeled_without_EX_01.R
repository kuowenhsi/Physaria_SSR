args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}

repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source_script <- file.path(repo_root, "script", "tables", "plot_Table_5_heatmap.R")

Sys.setenv(
  PHYSARIA_WEN_LABEL_FILE = file.path(
    repo_root,
    "data",
    "metadata",
    "Physaria_globosa_population_counts_from_individuals_relabeled.csv"
  )
)

# Load the original calculation and plotting functions without running its
# output block.
script_text <- paste(readLines(source_script, warn = FALSE), collapse = "\n")
script_text <- strsplit(script_text, "\ntable_5 <- build_pairwise_table", fixed = TRUE)[[1]][1]
script_text <- sub(
  "repo_root <- normalizePath(file.path(script_dir, \"..\", \"..\"))",
  "repo_root <- get(\"repo_root\", envir = .GlobalEnv)",
  script_text,
  fixed = TRUE
)

plot_env <- new.env(parent = globalenv())
eval(parse(text = script_text), envir = plot_env)

with(plot_env, {
  table_5_full <- build_pairwise_table(read_pg_all(input_file))
  table_5_long_full <- tidy_pairwise_table(table_5_full)

  table_5_without_ex_01 <- table_5_full %>%
    filter(population != "EX_01") %>%
    select(-any_of("EX_01"))
  table_5_long <- tidy_pairwise_table(table_5_without_ex_01)

  heatmap_plot <- build_heatmap(table_5_long, fill_option = "cividis") +
    # Preserve the original full-data color mapping after removing the row
    # and column from the rendered heatmap.
    scale_fill_viridis_c(
      option = "cividis",
      limits = range(table_5_long_full$value, na.rm = TRUE),
      name = NULL
    )

  ggsave(
    file.path(
      repo_root,
      "figure",
      "Table_5_pairwise_FST_GST_heatmap_cividis_relabeled_without_EX_01.png"
    ),
    plot = heatmap_plot,
    width = 3.7,
    height = 3.5,
    dpi = 600
  )
})
