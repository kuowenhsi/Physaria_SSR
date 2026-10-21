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
  globosa_scatter_tbl_full <- population_tbl %>%
    filter(species == "P. globosa") %>%
    mutate(
      `Geographic area` = factor(`Geographic area`, levels = names(state_colors)),
      highlight_label = Population_plot %in% c("KY_02", "KY_05", "TR_02", "TNE_01", "TNW_02", "TNW_01")
    )

  globosa_scatter_tbl <- globosa_scatter_tbl_full %>%
    filter(Population_plot != "EX_01")

  p_globosa_ar_ap <- ggplot(
    globosa_scatter_tbl,
    aes(x = AR, y = AP, color = `Geographic area`)
  ) +
    geom_point(size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(
      aes(label = Population_plot),
      size = 2.4,
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    geom_blank(
      data = globosa_scatter_tbl_full,
      aes(x = AR, y = AP),
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(
        floor(min(globosa_scatter_tbl_full$AP, na.rm = TRUE)),
        ceiling(max(globosa_scatter_tbl_full$AP, na.rm = TRUE)),
        by = 1
      )
    ) +
    labs(
      x = "Allelic richness (AR)",
      y = "Private alleles (AP)",
      color = "State"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  p_globosa_he_ho_fb <- ggplot(
    globosa_scatter_tbl,
    aes(x = He, y = Ho, fill = FB)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray65", linewidth = 0.4) +
    geom_point(shape = 21, color = "black", size = 2.8, alpha = 0.9, stroke = 0.35) +
    ggrepel::geom_text_repel(
      data = globosa_scatter_tbl %>% filter(!highlight_label),
      aes(label = Population_plot),
      size = 2.4,
      color = "black",
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    geom_text(
      data = globosa_scatter_tbl %>% filter(highlight_label),
      aes(label = Population_plot),
      size = 2.4,
      color = "#c00000",
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_blank(
      data = globosa_scatter_tbl_full,
      aes(x = He, y = Ho),
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#fddbc7",
      high = "#b2182b",
      midpoint = 0,
      limits = range(globosa_scatter_tbl_full$FB, na.rm = TRUE)
    ) +
    labs(
      x = "Expected heterozygosity (He)",
      y = "Observed heterozygosity (Ho)",
      fill = "FB"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA)
    )

  p_globosa_scatter_combined <- cowplot::plot_grid(
    p_globosa_ar_ap,
    p_globosa_he_ho_fb,
    nrow = 1,
    align = "h",
    axis = "tb",
    labels = c("C", "D")
  )

  save_plot(
    p_globosa_scatter_combined,
    "Table_3_P_globosa_scatter_combined.png",
    width = 7,
    height = 3.5
  )
})
