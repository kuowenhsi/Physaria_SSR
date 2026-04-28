library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
  dirname(normalizePath(script_path))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(repo_root, "script", "utils", "utils_wen_labels.R"))

input_file <- file.path(repo_root, "data", "Table_3.csv")
table_4_file <- file.path(repo_root, "data", "Table_4.csv")
output_dir <- file.path(repo_root, "figure", "Table_3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
wen_lookup <- read_wen_label_lookup(repo_root)

read_table_3 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  tbl <- readr::read_csv(I(data_lines), na = c("", "NA"), show_col_types = FALSE)
  numeric_cols <- c("n", "A", "AR", "AP", "Ho", "He", "F", "% null", "FB")

  tbl %>%
    mutate(across(all_of(numeric_cols), as.numeric)) %>%
    mutate(
      row_type = case_when(
        str_detect(Population, "^Average across") ~ "average",
        if_all(all_of(c("n", "A", "AR", "AP", "Ho", "He", "F", "% null", "FB")), is.na) ~ "section",
        TRUE ~ "population"
      ),
      species_header = if_else(row_type == "section", Population, NA_character_)
    ) %>%
    tidyr::fill(species_header, .direction = "down") %>%
    mutate(
      species = species_header,
      species = case_when(
        str_detect(Population, "^Average across P\\. globosa") ~ "P. globosa",
        str_detect(Population, "^Average across P\\. ouachitensis") ~ "P. ouachitensis",
        str_detect(Population, "^Average across P\\. filiformis") ~ "P. filiformis",
        str_detect(Population, "^Average across P\\. gracilis") ~ "P. gracilis",
        TRUE ~ species
      )
    ) %>%
    select(-species_header)
}

extract_significance <- function(x) {
  case_when(
    str_detect(x, fixed("***")) ~ "P < 0.0001",
    str_detect(x, fixed("**")) ~ "P < 0.001",
    str_detect(x, fixed("*")) ~ "P < 0.01",
    TRUE ~ "Not significant"
  )
}

parse_stat_value <- function(x) {
  readr::parse_number(x)
}

read_table_4 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
    transmute(
      Population,
      `Geographic area`,
      het_excess_z_raw = `Excess in heterozygosity Z-test`,
      het_excess_w_raw = `Excess in heterozygosity Wilcoxon signed-rank test`,
      mratio_z_raw = `Deficiency in M-Ratio Z-test`,
      mratio_w_raw = `Deficiency in M-Ratio Wilcoxon signed-rank test`,
      het_excess_z = parse_stat_value(het_excess_z_raw),
      het_excess_w = parse_stat_value(het_excess_w_raw),
      mratio_z = parse_stat_value(mratio_z_raw),
      mratio_w = parse_stat_value(mratio_w_raw),
      het_excess_z_sig = extract_significance(het_excess_z_raw),
      het_excess_w_sig = extract_significance(het_excess_w_raw),
      mratio_z_sig = extract_significance(mratio_z_raw),
      mratio_w_sig = extract_significance(mratio_w_raw)
    ) %>%
    mutate(
      het_bottleneck_sig = case_when(
        het_excess_z_sig == "P < 0.0001" | het_excess_w_sig == "P < 0.0001" ~ "P < 0.0001",
        het_excess_z_sig == "P < 0.001" | het_excess_w_sig == "P < 0.001" ~ "P < 0.001",
        het_excess_z_sig == "P < 0.01" | het_excess_w_sig == "P < 0.01" ~ "P < 0.01",
        TRUE ~ "Not significant"
      ),
      mratio_bottleneck_sig = case_when(
        mratio_z_sig == "P < 0.0001" | mratio_w_sig == "P < 0.0001" ~ "P < 0.0001",
        mratio_z_sig == "P < 0.001" | mratio_w_sig == "P < 0.001" ~ "P < 0.001",
        mratio_z_sig == "P < 0.01" | mratio_w_sig == "P < 0.01" ~ "P < 0.01",
        TRUE ~ "Not significant"
      )
    )
}

table_3 <- read_table_3(input_file)
table_4 <- read_table_4(table_4_file)

population_tbl <- table_3 %>%
  filter(row_type == "population") %>%
  left_join(table_4, by = c("Population", "Geographic area")) %>%
  mutate(
    species = factor(species, levels = c("P. globosa", "P. ouachitensis", "P. filiformis", "P. gracilis")),
    `Geographic area` = factor(`Geographic area`),
    Population_plot = apply_wen_label_from_canonical(Population, wen_lookup),
    het_bottleneck_sig = replace_na(het_bottleneck_sig, "Not tested"),
    mratio_bottleneck_sig = replace_na(mratio_bottleneck_sig, "Not tested")
  )

average_tbl <- table_3 %>%
  filter(row_type == "average") %>%
  mutate(
    species = factor(species, levels = c("P. globosa", "P. ouachitensis", "P. filiformis", "P. gracilis"))
  )

species_colors <- c(
  "P. globosa" = "#1f78b4",
  "P. ouachitensis" = "#e31a1c",
  "P. filiformis" = "#33a02c",
  "P. gracilis" = "#ff7f00"
)

state_colors <- c(
  "IN" = "#1b9e77",
  "KY" = "#d95f02",
  "EX-S" = "#7f7f7f",
  "TNE" = "#7570b3",
  "TNW" = "#e7298a"
)

sig_shapes <- c(
  "Not tested" = 21,
  "Not significant" = 21,
  "P < 0.01" = 24,
  "P < 0.001" = 22,
  "P < 0.0001" = 23
)

sig_fills <- c(
  "Not tested" = "white",
  "Not significant" = "white",
  "P < 0.01" = "#fdd835",
  "P < 0.001" = "#fb8c00",
  "P < 0.0001" = "#d7301f"
)

metric_labels <- c(
  A = "Alleles per locus (A)",
  AR = "Allelic richness (AR)",
  AP = "Private alleles (AP)",
  Ho = "Observed heterozygosity (Ho)",
  He = "Expected heterozygosity (He)",
  F = "Fixation index (F)",
  `% null` = "Null alleles (%)",
  FB = "Adjusted fixation index (FB)"
)

save_plot <- function(plot_obj, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

# 1. Diversity metrics across populations
diversity_long <- population_tbl %>%
  select(Population, Population_plot, species, `Geographic area`, A, AR, AP) %>%
  pivot_longer(cols = c(A, AR, AP), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("A", "AR", "AP"), labels = metric_labels[c("A", "AR", "AP")]),
    Population_plot = fct_reorder(Population_plot, value, .desc = FALSE)
  )

p_diversity <- ggplot(diversity_long, aes(x = value, y = Population_plot, color = species)) +
  geom_segment(aes(x = 0, xend = value, y = Population_plot, yend = Population_plot), linewidth = 0.35, alpha = 0.5) +
  geom_point(size = 2.2) +
  facet_wrap(~metric, scales = "free_x", ncol = 1) +
  scale_color_manual(values = species_colors) +
  labs(
    title = "Microsatellite diversity metrics by population",
    x = NULL,
    y = NULL,
    color = "Species"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

save_plot(p_diversity, "Table_3_diversity_metrics.png", width = 10, height = 11)

# 2. Observed vs expected heterozygosity per population
het_long <- population_tbl %>%
  mutate(Population_plot = fct_reorder(Population_plot, He)) %>%
  select(Population_plot, species, Ho, He, het_bottleneck_sig) %>%
  pivot_longer(cols = c(Ho, He), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, Ho = "Observed (Ho)", He = "Expected (He)")
  )

p_het_profile <- ggplot(het_long, aes(x = value, y = Population_plot, color = metric)) +
  geom_line(aes(group = Population_plot), color = "gray70", linewidth = 0.4) +
  geom_point(aes(shape = het_bottleneck_sig, fill = het_bottleneck_sig), size = 2.6, stroke = 0.35) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Observed (Ho)" = "#1b9e77", "Expected (He)" = "#d95f02")) +
  scale_shape_manual(values = sig_shapes, drop = FALSE) +
  scale_fill_manual(values = sig_fills, drop = FALSE) +
  labs(
    title = "Observed versus expected heterozygosity",
    subtitle = "Point symbols show heterozygosity-excess bottleneck significance from Table 4",
    x = "Heterozygosity",
    y = NULL,
    color = NULL,
    shape = "Het. excess test",
    fill = "Het. excess test"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

save_plot(p_het_profile, "Table_3_heterozygosity_profile.png", width = 11, height = 8.5)

# 3. Ho vs He scatter
p_het_scatter <- ggplot(population_tbl, aes(x = He, y = Ho, color = species, size = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55") +
  geom_point(aes(shape = het_bottleneck_sig, fill = het_bottleneck_sig), alpha = 0.95, stroke = 0.4) +
  ggrepel::geom_text_repel(aes(label = Population_plot), size = 2.7, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = species_colors) +
  scale_size_continuous(range = c(2, 7)) +
  scale_shape_manual(values = sig_shapes, drop = FALSE) +
  scale_fill_manual(values = sig_fills, drop = FALSE) +
  labs(
    title = "Observed versus expected heterozygosity across populations",
    subtitle = "Points below the 1:1 line indicate heterozygote deficiency; symbols show Table 4 heterozygosity-excess significance",
    x = "Expected heterozygosity (He)",
    y = "Observed heterozygosity (Ho)",
    color = "Species",
    size = "Sample size",
    shape = "Het. excess test",
    fill = "Het. excess test"
  ) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

save_plot(p_het_scatter, "Table_3_heterozygosity_scatter.png", width = 9, height = 7)

# 4. Fixation and null-allele diagnostics
diagnostic_long <- population_tbl %>%
  select(Population, Population_plot, species, F, FB, `% null`) %>%
  pivot_longer(cols = c(F, FB, `% null`), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("F", "FB", "% null"), labels = metric_labels[c("F", "FB", "% null")]),
    Population_plot = fct_reorder(Population_plot, value)
  )

p_diagnostics <- ggplot(diagnostic_long, aes(x = value, y = Population_plot, color = species)) +
  geom_vline(xintercept = 0, color = "gray80", linewidth = 0.4) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_x", ncol = 1) +
  scale_color_manual(values = species_colors) +
  labs(
    title = "Fixation and null-allele diagnostics",
    x = NULL,
    y = NULL,
    color = "Species"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

save_plot(p_diagnostics, "Table_3_fixation_null_diagnostics.png", width = 10, height = 11)

# 5. Bottleneck statistics from Table 4
bottleneck_long <- population_tbl %>%
  filter(species == "P. globosa", !is.na(het_excess_z)) %>%
  select(
    Population,
    Population_plot,
    `Geographic area`,
    het_excess_z,
    het_excess_w,
    mratio_z,
    mratio_w,
    het_excess_z_sig,
    het_excess_w_sig,
    mratio_z_sig,
    mratio_w_sig
  ) %>%
  pivot_longer(
    cols = c(het_excess_z, het_excess_w, mratio_z, mratio_w),
    names_to = "test",
    values_to = "statistic"
  ) %>%
  mutate(
    sig = case_when(
      test == "het_excess_z" ~ het_excess_z_sig,
      test == "het_excess_w" ~ het_excess_w_sig,
      test == "mratio_z" ~ mratio_z_sig,
      test == "mratio_w" ~ mratio_w_sig
    ),
    test = factor(
      test,
      levels = c("het_excess_z", "het_excess_w", "mratio_z", "mratio_w"),
      labels = c(
        "Het. excess Z-test",
        "Het. excess Wilcoxon",
        "M-ratio deficiency Z-test",
        "M-ratio deficiency Wilcoxon"
      )
    ),
    Population_plot = fct_reorder(Population_plot, statistic)
  )

p_bottleneck <- ggplot(bottleneck_long, aes(x = statistic, y = Population_plot)) +
  geom_vline(xintercept = 0, color = "gray82", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = statistic, y = Population_plot, yend = Population_plot), color = "gray75", linewidth = 0.35) +
  geom_point(aes(fill = sig, shape = sig), color = "black", size = 2.6, stroke = 0.3) +
  facet_wrap(~test, scales = "free_x", ncol = 1) +
  scale_shape_manual(values = sig_shapes[c("Not significant", "P < 0.01", "P < 0.001", "P < 0.0001")], drop = FALSE) +
  scale_fill_manual(values = sig_fills[c("Not significant", "P < 0.01", "P < 0.001", "P < 0.0001")], drop = FALSE) +
  labs(
    title = "Table 4 bottleneck test statistics for P. globosa",
    subtitle = "Symbols encode reported significance levels",
    x = "Test statistic",
    y = NULL,
    shape = "Significance",
    fill = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

save_plot(p_bottleneck, "Table_3_bottleneck_statistics.png", width = 10, height = 11)

# 6. Standardized metric heatmap
heatmap_long <- population_tbl %>%
  select(Population, Population_plot, species, A, AR, AP, Ho, He, F, `% null`, FB) %>%
  pivot_longer(cols = c(A, AR, AP, Ho, He, F, `% null`, FB), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(z_value = as.numeric(scale(value))) %>%
  ungroup() %>%
  mutate(
    metric = factor(metric, levels = c("A", "AR", "AP", "Ho", "He", "F", "% null", "FB"), labels = metric_labels[c("A", "AR", "AP", "Ho", "He", "F", "% null", "FB")])
  )

p_heatmap <- ggplot(heatmap_long, aes(x = metric, y = Population_plot, fill = z_value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Z-score") +
  labs(
    title = "Standardized microsatellite summary heatmap",
    subtitle = "Values are scaled within each metric to highlight relative highs and lows",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

save_plot(p_heatmap, "Table_3_metric_heatmap.png", width = 11, height = 8.5)

# 7. Species-level summary from manuscript averages
average_long <- average_tbl %>%
  select(species, A, AR, AP, Ho, He, F, `% null`, FB) %>%
  pivot_longer(cols = -species, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("A", "AR", "AP", "Ho", "He", "F", "% null", "FB"), labels = metric_labels[c("A", "AR", "AP", "Ho", "He", "F", "% null", "FB")])
  )

p_species_avg <- ggplot(average_long, aes(x = species, y = value, fill = species)) +
  geom_col(width = 0.72) +
  facet_wrap(~metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = species_colors, guide = "none") +
  labs(
    title = "Species-level averages reported in Table 3",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p_species_avg, "Table_3_species_averages.png", width = 12, height = 8)

# 8. P. globosa allelic richness versus private alleles
globosa_scatter_tbl <- population_tbl %>%
  filter(species == "P. globosa") %>%
  mutate(
    `Geographic area` = factor(`Geographic area`, levels = names(state_colors)),
    highlight_label = Population_plot %in% c("KY_02", "KY_05", "EX_02", "TNE_01", "TNW_02", "TNW_01")
  )

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
  scale_color_manual(values = state_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(
      floor(min(globosa_scatter_tbl$AP, na.rm = TRUE)),
      ceiling(max(globosa_scatter_tbl$AP, na.rm = TRUE)),
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

save_plot(p_globosa_ar_ap, "Table_3_P_globosa_AR_AP_scatter.png", width = 3.5, height = 3.5)

# 9. P. globosa expected versus observed heterozygosity colored by corrected fixation index
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
  ggrepel::geom_text_repel(
    data = globosa_scatter_tbl %>% filter(highlight_label),
    aes(label = Population_plot),
    size = 2.4,
    color = "#c00000",
    fontface = "bold",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "#fddbc7",
    high = "#b2182b",
    midpoint = 0
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

save_plot(p_globosa_he_ho_fb, "Table_3_P_globosa_He_Ho_FB_scatter.png", width = 3.5, height = 3.5)


p_globosa_scatter_combined <- cowplot::plot_grid(
  p_globosa_ar_ap,
  p_globosa_he_ho_fb,
  nrow = 1,
  align = "h",
  axis = "tb",
  labels = c("C", "D")
)

save_plot(p_globosa_scatter_combined, "Table_3_P_globosa_scatter_combined.png", width = 7, height = 3.5)

all_species_scatter_tbl <- population_tbl %>%
  mutate(
    Population_display = dplyr::coalesce(Wen_label, Population_plot)
  )

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

save_plot(p_all_species_ar_ap, "Table_3_all_species_AR_AP_scatter.png", width = 3.5, height = 3.5)

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
  scale_color_manual(values = species_colors, drop = FALSE) +
  labs(
    x = "Expected heterozygosity (He)",
    y = "Observed heterozygosity (Ho)",
    color = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none"
  )

save_plot(p_all_species_he_ho, "Table_3_all_species_He_Ho_scatter.png", width = 3.5, height = 3.5)

p_all_species_scatter_combined <- cowplot::plot_grid(
  p_all_species_ar_ap,
  p_all_species_he_ho,
  nrow = 1,
  align = "h",
  axis = "tb",
  labels = c("A", "B")
)

save_plot(p_all_species_scatter_combined, "Table_3_all_species_scatter_combined.png", width = 7, height = 3.5)

p_all_scatter_four_panel <- cowplot::plot_grid(
  p_all_species_ar_ap,
  p_all_species_he_ho,
  p_globosa_ar_ap,
  p_globosa_he_ho_fb,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  labels = c("A", "B", "C", "D")
)

save_plot(p_all_scatter_four_panel, "Table_3_all_scatter_four_panel.png", width = 7, height = 7)

readr::write_csv(population_tbl, file.path(output_dir, "Table_3_population_clean.csv"))
readr::write_csv(average_tbl, file.path(output_dir, "Table_3_species_averages_clean.csv"))
