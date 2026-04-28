library(tidyverse)
library(clue)
library(ggh4x)
library(ggokabeito)

# setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/Manuscript/Physaria_SSR/script")

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

structure_dir <- file.path(repo_root, "results", "structureHarvester")
sample_file <- file.path(repo_root, "data", "PG_All_structure.txt")
output_dir <- file.path(repo_root, "figure", "STRUCTURE_individual_allK")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
wen_lookup <- read_wen_label_lookup(repo_root)

read_sample_metadata <- function(path) {
  lines <- readr::read_lines(path)
  sample_lines <- lines[-1]

  tibble(raw_line = sample_lines) %>%
    mutate(
      parts = str_split(str_squish(raw_line), "\\s+"),
      sample_id = map_chr(parts, 1),
      population = str_extract(sample_id, "Pg\\d+"),
      population_label = apply_wen_label_from_pg(population, wen_lookup),
      state = case_when(
        population == "Pg346" ~ "IN",
        population %in% c("Pg347", "Pg349", "Pg350", "Pg351", "Pg352", "Pg354") ~ "KY",
        population %in% c("Pg348", "Pg353", "Pg355") ~ "EX-S",
        population %in% c("Pg356", "Pg357") ~ "TNE",
        population %in% c("Pg358", "Pg359", "Pg360", "Pg361", "Pg362") ~ "TNW",
        TRUE ~ "Unknown"
      ),
      sample_order = row_number()
    ) %>%
    select(sample_id, population, population_label, state, sample_order)
}

read_indfile <- function(path, k) {
  lines <- readr::read_lines(path) %>%
    keep(~ str_detect(.x, ":"))

  parsed <- map_dfr(lines, function(line) {
    line <- str_squish(line)
    pieces <- str_split_fixed(line, " : ", 2)
    left <- str_split(pieces[1], "\\s+")[[1]]
    qvals <- as.numeric(str_split(pieces[2], "\\s+")[[1]])

    tibble(
      individual_index = as.integer(left[1]),
      pop_code = as.integer(left[length(left)]),
      q = list(qvals)
    )
  })

  qmat <- do.call(rbind, parsed$q)
  colnames(qmat) <- paste0("cluster_", seq_len(k))

  bind_cols(parsed %>% select(individual_index, pop_code), as_tibble(qmat))
}

align_and_average_runs <- function(ind_tbl, sample_meta, k, n_runs = 5) {
  n_samples <- nrow(sample_meta)
  stopifnot(nrow(ind_tbl) == n_samples * n_runs)

  run_list <- map(seq_len(n_runs), function(run_id) {
    start_idx <- (run_id - 1) * n_samples + 1
    end_idx <- run_id * n_samples
    qmat <- ind_tbl[start_idx:end_idx, paste0("cluster_", seq_len(k))] %>%
      as.matrix()
    list(run_id = run_id, qmat = qmat)
  })

  ref <- run_list[[1]]$qmat
  aligned <- list(ref)

  if (k > 1) {
    for (run_id in 2:n_runs) {
      current <- run_list[[run_id]]$qmat
      similarity <- cor(ref, current, use = "pairwise.complete.obs")
      similarity[is.na(similarity)] <- -1
      similarity <- similarity - min(similarity) + 1e-9
      perm <- solve_LSAP(similarity, maximum = TRUE)
      aligned[[run_id]] <- current[, perm, drop = FALSE]
    }
  } else {
    for (run_id in 2:n_runs) {
      aligned[[run_id]] <- run_list[[run_id]]$qmat
    }
  }

  avg_q <- Reduce(`+`, aligned) / n_runs
  row_totals <- rowSums(avg_q, na.rm = TRUE)
  needs_rescale <- !is.na(row_totals) & abs(row_totals - 1) > sqrt(.Machine$double.eps)

  if (any(needs_rescale)) {
    avg_q[needs_rescale, ] <- avg_q[needs_rescale, , drop = FALSE] / row_totals[needs_rescale]
  }

  colnames(avg_q) <- paste0("cluster_", seq_len(k))

  bind_cols(sample_meta, as_tibble(avg_q))
}

build_barplot <- function(avg_tbl, k) {
  cluster_cols <- paste0("cluster_", seq_len(k))

  population_levels <- avg_tbl %>%
    mutate(
      state = factor(state, levels = unique(state)),
      population_number = readr::parse_number(population)
    ) %>%
    distinct(state, population_number, population_label) %>%
    arrange(state, population_number) %>%
    pull(population_label)

  plot_tbl <- avg_tbl %>%
    mutate(
      state = factor(state, levels = unique(state)),
      population_number = readr::parse_number(population)
    ) %>%
    arrange(state, population_number, sample_order) %>%
    mutate(
      population_label = factor(population_label, levels = population_levels, labels = str_split_i(population_levels, "_", 2)),
      sample_id = factor(sample_id, levels = rev(sample_id))
    )

  long_tbl <- plot_tbl %>%
    pivot_longer(
      cols = all_of(cluster_cols),
      names_to = "cluster",
      values_to = "membership"
    ) %>%
    mutate(
      cluster = factor(cluster, levels = cluster_cols)
    )

  cluster_palette <- grDevices::hcl.colors(k, palette = "Set 3", rev = FALSE)
  names(cluster_palette) <- cluster_cols

  ggplot(long_tbl, aes(x = sample_id, y = membership, fill = cluster)) +
    geom_col(width = 1, color = NA, show.legend = FALSE) +
    ggh4x::facet_nested(
      cols = vars(state, population_label),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_manual(values = cluster_palette, name = "Cluster") +
    scale_y_continuous(
      limits = c(0, 1.001),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Mean ancestry composition"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y.left = element_text(size = 7),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      panel.spacing.x = grid::unit(0.05, "lines"),
      strip.placement = "outside",
      strip.text.x.top = element_text(angle = 0),
      panel.spacing.y = unit(0, "lines")
    )
}

sample_meta <- read_sample_metadata(sample_file)

for (k in 1:13) {

  ind_path <- file.path(structure_dir, paste0("K", k, ".indfile"))
  ind_tbl <- read_indfile(ind_path, k)
  avg_tbl <- align_and_average_runs(ind_tbl, sample_meta, k = k, n_runs = 5)
  p <- build_barplot(avg_tbl, k)
  p

  ggsave(
    file.path(output_dir, sprintf("STRUCTURE_K%02d_individual_mean.png", k)),
    plot = p,
    width = 7.25,
    height = 2,
    dpi = 600
  )

  readr::write_csv(
    avg_tbl,
    file.path(output_dir, sprintf("STRUCTURE_K%02d_individual_mean.csv", k))
  )
}

build_multi_k_barplot <- function(avg_tbl_list, k_values) {
  max_k <- max(k_values)
  cluster_cols <- paste0("cluster_", seq_len(max_k))

  population_levels <- avg_tbl_list[[1]] %>%
    mutate(
      state = factor(state, levels = unique(state)),
      population_number = readr::parse_number(population)
    ) %>%
    distinct(state, population_number, population_label) %>%
    arrange(state, population_number) %>%
    pull(population_label)

  long_tbl <- purrr::map2_dfr(avg_tbl_list, k_values, function(avg_tbl, k) {
    current_cluster_cols <- paste0("cluster_", seq_len(k))

    avg_tbl %>%
      mutate(
        state = factor(state, levels = unique(state)),
        population_number = readr::parse_number(population)
      ) %>%
      arrange(state, population_number, sample_order) %>%
      mutate(
        population_label = factor(
          population_label,
          levels = population_levels,
          labels = str_split_i(population_levels, "_", 2)
        ),
        sample_id = factor(sample_id, levels = rev(sample_id)),
        K = factor(paste0("K = ", k), levels = paste0("K = ", k_values))
      ) %>%
      pivot_longer(
        cols = all_of(current_cluster_cols),
        names_to = "cluster",
        values_to = "membership"
      )
  }) %>%
    mutate(
      cluster = factor(cluster, levels = cluster_cols)
    )

  ggplot(long_tbl, aes(x = sample_id, y = membership, fill = cluster)) +
    geom_col(width = 1, color = NA, show.legend = FALSE) +
    ggh4x::facet_nested(
      cols = vars(state, population_label),
      rows = vars(K),
      scales = "free_x",
      space = "free_x"
    ) +
    ggokabeito::scale_fill_okabe_ito(name = "Cluster", drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1.001),
      breaks = seq(0.2, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Mean ancestry composition"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y.left = element_text(size = 7),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      panel.spacing.x = grid::unit(0.05, "lines"),
      panel.spacing.y = grid::unit(0.15, "lines"),
      strip.placement = "outside",
      strip.text.x.top = element_text(angle = 0)
    )
}

k_values_combined <- 2:5
avg_tbl_list_combined <- purrr::map(k_values_combined, function(k) {
  ind_path <- file.path(structure_dir, paste0("K", k, ".indfile"))
  ind_tbl <- read_indfile(ind_path, k)
  align_and_average_runs(ind_tbl, sample_meta, k = k, n_runs = 5)
})

p_multi_k <- build_multi_k_barplot(avg_tbl_list_combined, k_values_combined)

ggsave(
  file.path(output_dir, "STRUCTURE_K02_K05_individual_mean_combined.png"),
  plot = p_multi_k,
  width = 7.25,
  height = 4,
  dpi = 600
)

########

k_values_combined <- 2:8
avg_tbl_list_combined <- purrr::map(k_values_combined, function(k) {
  ind_path <- file.path(structure_dir, paste0("K", k, ".indfile"))
  ind_tbl <- read_indfile(ind_path, k)
  align_and_average_runs(ind_tbl, sample_meta, k = k, n_runs = 5)
})

p_multi_k <- build_multi_k_barplot(avg_tbl_list_combined, k_values_combined)

ggsave(
  file.path(output_dir, "STRUCTURE_K02_K08_individual_mean_combined.png"),
  plot = p_multi_k,
  width = 7.25,
  height = 8.5,
  dpi = 600
)
