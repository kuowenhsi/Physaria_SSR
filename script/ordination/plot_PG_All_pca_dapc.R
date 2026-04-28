library(tidyverse)
library(adegenet)
library(ggrepel)

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

input_file <- file.path(repo_root, "data", "PG_All.csv")
output_dir <- file.path(repo_root, "figure", "PG_All_ordination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
wen_lookup <- read_wen_label_lookup(repo_root)

read_pg_all <- function(path) {
  raw <- readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)

  locus_names <- raw[3, -(1:2)] %>%
    unlist(use.names = FALSE) %>%
    as.character()
  locus_names[is.na(locus_names)] <- ""

  allele_names <- make.unique(locus_names[locus_names != ""])
  allele_col_names <- paste0("allele_", seq_len(length(locus_names)))
  sample_tbl <- raw[-c(1:3), ] %>%
    setNames(c("sample", "site", allele_col_names)) %>%
    mutate(
      sample = as.character(sample),
      site = as.character(site),
      site = str_replace(site, "^p", ""),
      population = str_extract(sample, "^Pg\\d+"),
      species = "P. globosa",
      state = case_when(
        population == "Pg346" ~ "IN",
        population %in% c("Pg347", "Pg349", "Pg350", "Pg351", "Pg352", "Pg354") ~ "KY",
        population %in% c("Pg348", "Pg353", "Pg355") ~ "EX-S",
        population %in% c("Pg356", "Pg357") ~ "TNE",
        population %in% c("Pg358", "Pg359", "Pg360", "Pg361", "Pg362") ~ "TNW",
        TRUE ~ "Unknown"
      )
    )

  allele_cols <- allele_col_names

  allele_matrix <- sample_tbl %>%
    select(all_of(allele_cols)) %>%
    mutate(across(everything(), as.character))

  locus_tbl <- map_dfc(seq_along(allele_names), function(i) {
    allele_1 <- allele_matrix[[2 * i - 1]]
    allele_2 <- allele_matrix[[2 * i]]

    genotype <- if_else(
      is.na(allele_1) | is.na(allele_2) | allele_1 == "0" | allele_2 == "0",
      NA_character_,
      paste0(allele_1, "/", allele_2)
    )

    tibble(!!allele_names[i] := genotype)
  })

  metadata <- sample_tbl %>%
    mutate(population_label = apply_wen_label_from_pg(population, wen_lookup)) %>%
    select(sample, site, population, population_label, species, state)

  list(metadata = metadata, loci = locus_tbl)
}

build_genind <- function(parsed_tbl) {
  adegenet::df2genind(
    parsed_tbl$loci,
    sep = "/",
    ploidy = 2,
    ind.names = parsed_tbl$metadata$sample,
    pop = parsed_tbl$metadata$population,
    NA.char = "NA",
    type = "codom"
  )
}

run_pca <- function(genind_obj, metadata) {
  geno_tab <- adegenet::tab(genind_obj, NA.method = "mean")
  pca_obj <- prcomp(geno_tab, center = TRUE, scale. = FALSE)

  explained <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)

  ind_scores <- as_tibble(pca_obj$x[, 1:2, drop = FALSE]) %>%
    setNames(c("Axis1", "Axis2")) %>%
    bind_cols(metadata)

  pop_scores <- ind_scores %>%
    group_by(population, population_label, species, state) %>%
    summarise(
      Axis1 = mean(Axis1),
      Axis2 = mean(Axis2),
      n = n(),
      .groups = "drop"
    )

  list(
    ind_scores = ind_scores,
    pop_scores = pop_scores,
    explained = explained
  )
}

run_dapc <- function(genind_obj, metadata) {
  n_ind <- nInd(genind_obj)
  n_pop <- nPop(genind_obj)
  n_var <- ncol(adegenet::tab(genind_obj, NA.method = "mean"))

  n_pca <- max(2, min(40, n_ind - 1, n_var - 1))
  n_da <- max(1, min(10, n_pop - 1))

  dapc_obj <- adegenet::dapc(genind_obj, pop(genind_obj), n.pca = n_pca, n.da = n_da)

  ind_scores <- as_tibble(dapc_obj$ind.coord[, 1:2, drop = FALSE]) %>%
    setNames(c("Axis1", "Axis2")) %>%
    bind_cols(metadata)

  pop_scores <- ind_scores %>%
    group_by(population, population_label, species, state) %>%
    summarise(
      Axis1 = mean(Axis1),
      Axis2 = mean(Axis2),
      n = n(),
      .groups = "drop"
    )

  list(
    dapc_obj = dapc_obj,
    ind_scores = ind_scores,
    pop_scores = pop_scores,
    n_pca = n_pca,
    n_da = n_da
  )
}

state_colors <- c(
  "IN" = "#1b9e77",
  "KY" = "#d95f02",
  "EX-S" = "#7f7f7f",
  "TNE" = "#7570b3",
  "TNW" = "#e7298a",
  "Unknown" = "#666666"
)

plot_individual_ordination <- function(scores, method, axis_labels, subtitle) {
  ggplot(scores, aes(x = Axis1, y = Axis2, color = state)) +
    geom_point(size = 1.6, alpha = 0.25) +
    stat_ellipse(aes(group = population), linewidth = 0.35, alpha = 0.35, show.legend = FALSE) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = axis_labels[1],
      y = axis_labels[2],
      color = "State"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

plot_individual_ordination_with_population_overlay <- function(scores, pop_scores, method, axis_labels, subtitle) {
  ggplot(scores, aes(x = Axis1, y = Axis2)) +
    geom_point(
      data = scores %>% filter(state != "EX-S"),
      color = "gray90",
      size = 1.6,
      alpha = 0.25
    ) +
    geom_point(
      data = scores %>% filter(state == "EX-S"),
      color = "gray90",
      shape = 4,
      size = 1.9,
      stroke = 0.55,
      alpha = 0.25
    ) +
    geom_point(
      data = pop_scores,
      aes(x = Axis1, y = Axis2, fill = state),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.98,
      size = 3.2
    ) +
    ggrepel::geom_text_repel(
      data = pop_scores,
      aes(x = Axis1, y = Axis2, label = population_label, color = state),
      inherit.aes = FALSE,
      size = 3.1,
      fontface = "bold",
      max.overlaps = 50,
      seed = 42
    ) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    scale_fill_manual(values = state_colors, drop = FALSE, guide = "none") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = axis_labels[1],
      y = axis_labels[2],
      color = "State"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}

plot_individual_ordination_with_population_overlay_statecolor <- function(scores, pop_scores, method, axis_labels, subtitle) {
  ggplot(scores, aes(x = Axis1, y = Axis2)) +
    geom_point(
      data = scores %>% filter(state != "EX-S"),
      aes(color = state),
      size = 1.6,
      alpha = 0.25
    ) +
    geom_point(
      data = scores %>% filter(state == "EX-S"),
      aes(color = state),
      shape = 4,
      size = 1.9,
      stroke = 0.55,
      alpha = 0.25
    ) +
    geom_point(
      data = pop_scores,
      aes(x = Axis1, y = Axis2, fill = state),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.98,
      size = 3.2
    ) +
    ggrepel::geom_text_repel(
      data = pop_scores,
      aes(x = Axis1, y = Axis2, label = population_label, color = state),
      inherit.aes = FALSE,
      size = 3.1,
      fontface = "bold",
      max.overlaps = 50,
      seed = 42
    ) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    scale_fill_manual(values = state_colors, drop = FALSE, guide = "none") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = axis_labels[1],
      y = axis_labels[2],
      color = "State"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}

plot_population_ordination <- function(scores, method, axis_labels, subtitle) {
  ggplot(scores, aes(x = Axis1, y = Axis2, label = population_label, color = state)) +
    geom_hline(yintercept = 0, color = "gray85", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "gray85", linewidth = 0.35) +
    geom_point(aes(size = n), alpha = 0.95) +
    ggrepel::geom_text_repel(size = 3.1, max.overlaps = 50) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    scale_size_continuous(range = c(3, 8), name = "Individuals") +
    labs(
      title = paste(method, "at the population level"),
      subtitle = subtitle,
      x = axis_labels[1],
      y = axis_labels[2],
      color = "State"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

save_plot <- function(plot_obj, filename, width = 3.5, height = 3.5) {
  ggsave(
    file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

parsed <- read_pg_all(input_file)
genind_globosa <- build_genind(parsed)

pca_globosa <- run_pca(genind_globosa, parsed$metadata)
dapc_globosa <- run_dapc(genind_globosa, parsed$metadata)

pca_axis_labels <- c(
  sprintf("PC1 (%.1f%%)", 100 * pca_globosa$explained[1]),
  sprintf("PC2 (%.1f%%)", 100 * pca_globosa$explained[2])
)

dapc_axis_labels <- c("DA1", "DA2")

p_pca_ind <- plot_individual_ordination_with_population_overlay(
  pca_globosa$ind_scores,
  pca_globosa$pop_scores,
  method = "PCA",
  axis_labels = pca_axis_labels,
  subtitle = "PG_All.csv contains only P. globosa individuals; Wen-label population centroids overlaid"
)

p_pca_ind_statecolor <- plot_individual_ordination_with_population_overlay_statecolor(
  pca_globosa$ind_scores,
  pca_globosa$pop_scores,
  method = "PCA",
  axis_labels = pca_axis_labels,
  subtitle = "PG_All.csv contains only P. globosa individuals; individual points colored by state with Wen-label population centroids overlaid"
)

p_pca_pop <- plot_population_ordination(
  pca_globosa$pop_scores,
  method = "PCA",
  axis_labels = pca_axis_labels,
  subtitle = "Population centroids from individual PCA scores"
)

p_dapc_ind <- plot_individual_ordination(
  dapc_globosa$ind_scores,
  method = "DAPC",
  axis_labels = dapc_axis_labels,
  subtitle = sprintf("Retained %s PCs and %s discriminant axes", dapc_globosa$n_pca, dapc_globosa$n_da)
)

p_dapc_pop <- plot_population_ordination(
  dapc_globosa$pop_scores,
  method = "DAPC",
  axis_labels = dapc_axis_labels,
  subtitle = "Population centroids from individual DAPC scores"
)

save_plot(p_pca_ind, "PG_All_P_globosa_PCA_individual.png")
save_plot(p_pca_ind_statecolor, "PG_All_P_globosa_PCA_individual_statecolor.png")
save_plot(p_pca_pop, "PG_All_P_globosa_PCA_population.png")
save_plot(p_dapc_ind, "PG_All_P_globosa_DAPC_individual.png")
save_plot(p_dapc_pop, "PG_All_P_globosa_DAPC_population.png")

readr::write_csv(parsed$metadata, file.path(output_dir, "PG_All_metadata_clean.csv"))
