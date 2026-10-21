args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}

repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source_script <- file.path(repo_root, "script", "tables", "plot_Table_3_summary.R")

Sys.setenv(
  PHYSARIA_WEN_LABEL_FILE = file.path(
    repo_root,
    "data",
    "metadata",
    "Physaria_globosa_population_counts_from_individuals_relabeled.csv"
  ),
  PHYSARIA_TABLE_4_FILE = file.path(repo_root, "data", "Table_4_relabeled.csv"),
  PHYSARIA_OUTPUT_SUFFIX = "_relabeled_without_EX_01"
)

# Load the original data preparation, palettes, and saving function without
# running its output block.
script_text <- paste(readLines(source_script, warn = FALSE), collapse = "\n")
script_text <- strsplit(script_text, "\n# 1. Diversity metrics", fixed = TRUE)[[1]][1]
script_text <- sub(
  "repo_root <- normalizePath(file.path(script_dir, \"..\", \"..\"))",
  "repo_root <- get(\"repo_root\", envir = .GlobalEnv)",
  script_text,
  fixed = TRUE
)

plot_env <- new.env(parent = globalenv())
eval(parse(text = script_text), envir = plot_env)

with(plot_env, {
  all_species_scatter_tbl_full <- population_tbl %>%
    mutate(Population_display = dplyr::coalesce(Wen_label, Population_plot))

  all_species_scatter_tbl <- all_species_scatter_tbl_full %>%
    filter(Population_display != "EX_01")

  p_all_species_ar_ap <- ggplot(
    all_species_scatter_tbl,
    aes(x = AR, y = AP, color = species)
  ) +
    geom_point(size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = all_species_scatter_tbl %>% filter(species != "P. globosa"),
      aes(label = Population_display),
      size = 2.3,
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    geom_blank(
      data = all_species_scatter_tbl_full,
      aes(x = AR, y = AP),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = species_colors, drop = FALSE) +
    labs(
      x = "Allelic richness (AR)",
      y = "Private alleles (AP)",
      color = "Species"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA)
    )

  p_all_species_he_ho <- ggplot(
    all_species_scatter_tbl,
    aes(x = He, y = Ho, color = species)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray65", linewidth = 0.4) +
    geom_point(size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = all_species_scatter_tbl %>% filter(species != "P. globosa"),
      aes(label = Population_display),
      size = 2.3,
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    geom_blank(
      data = all_species_scatter_tbl_full,
      aes(x = He, y = Ho),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = species_colors, drop = FALSE) +
    labs(
      x = "Expected heterozygosity (He)",
      y = "Observed heterozygosity (Ho)",
      color = "Species"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  p_all_species_scatter_combined <- cowplot::plot_grid(
    p_all_species_ar_ap,
    p_all_species_he_ho,
    nrow = 1,
    align = "h",
    axis = "tb",
    labels = c("A", "B")
  )

  save_plot(
    p_all_species_scatter_combined,
    "Table_3_all_species_scatter_combined.png",
    width = 7,
    height = 3.5
  )
})
