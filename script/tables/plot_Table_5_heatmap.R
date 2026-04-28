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

input_file <- file.path(repo_root, "data", "PG_All.csv")
output_png <- file.path(repo_root, "figure", "Table_5_pairwise_FST_GST_heatmap.png")
output_csv <- file.path(repo_root, "figure", "Table_5_pairwise_FST_GST_recomputed.csv")
wen_lookup <- read_wen_label_lookup(repo_root)

order_wen_labels <- function(labels) {
  state_order <- c("IN", "KY", "EX", "TNE", "TNW")

  tibble(label = labels) %>%
    mutate(
      state_group = stringr::str_extract(label, "^[A-Z]+"),
      state_group = factor(state_group, levels = state_order),
      numeric_id = readr::parse_number(label)
    ) %>%
    arrange(state_group, numeric_id, label) %>%
    pull(label)
}

read_pg_all <- function(path) {
  raw <- readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)

  locus_names <- raw[3, -(1:2)] %>%
    unlist(use.names = FALSE) %>%
    as.character()
  locus_names[is.na(locus_names)] <- ""

  allele_col_names <- paste0("allele_", seq_len(length(locus_names)))

  sample_tbl <- raw[-c(1:3), ] %>%
    setNames(c("sample", "site", allele_col_names)) %>%
    mutate(
      sample = as.character(sample),
      site = as.character(site),
      population = str_extract(sample, "^Pg\\d+")
    ) %>%
    left_join(
      wen_lookup %>%
        select(pg_population, canonical_population_code, wen_label),
      by = c("population" = "pg_population")
    ) %>%
    filter(!is.na(canonical_population_code))

  list(
    sample_tbl = sample_tbl,
    n_loci = length(locus_names) / 2
  )
}

calc_pair_metrics <- function(sample_tbl, pop1, pop2, n_loci) {
  hs_values <- c()
  ht_values <- c()

  for (i in seq_len(n_loci)) {
    allele_cols <- c(paste0("allele_", 2 * i - 1), paste0("allele_", 2 * i))

    pop_alleles <- map(c(pop1, pop2), function(pop) {
      sample_tbl %>%
        filter(canonical_population_code == pop) %>%
        select(all_of(allele_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.character() %>%
        discard(~ is.na(.x) || .x == "" || .x == "0")
    })
    names(pop_alleles) <- c(pop1, pop2)

    if (any(map_int(pop_alleles, length) == 0)) {
      next
    }

    pop_counts <- map(pop_alleles, table)
    all_alleles <- reduce(map(pop_counts, names), union)

    freq_vecs <- map(pop_counts, function(counts) {
      freqs <- setNames(rep(0, length(all_alleles)), all_alleles)
      freqs[names(counts)] <- as.numeric(counts) / sum(counts)
      freqs
    })

    h_values <- map_dbl(freq_vecs, ~ 1 - sum(.x^2))
    hs_locus <- mean(h_values)
    total_freq <- reduce(freq_vecs, `+`) / 2
    ht_locus <- 1 - sum(total_freq^2)

    hs_values <- c(hs_values, hs_locus)
    ht_values <- c(ht_values, ht_locus)
  }

  if (length(hs_values) == 0 || sum(ht_values, na.rm = TRUE) == 0) {
    return(c(fst = NA_real_, gst_prime = NA_real_))
  }

  hs_bar <- mean(hs_values, na.rm = TRUE)
  ht_bar <- mean(ht_values, na.rm = TRUE)
  fst <- (ht_bar - hs_bar) / ht_bar
  gst_prime <- if ((1 - hs_bar) > 0) {
    fst * (1 + hs_bar) / (1 - hs_bar)
  } else {
    NA_real_
  }

  c(
    fst = max(0, min(1, fst)),
    gst_prime = ifelse(is.na(gst_prime), NA_real_, max(0, min(1, gst_prime)))
  )
}

build_pairwise_table <- function(parsed_pg_all) {
  populations <- parsed_pg_all$sample_tbl %>%
    distinct(canonical_population_code, wen_label) %>%
    mutate(wen_label = factor(wen_label, levels = order_wen_labels(wen_label))) %>%
    arrange(wen_label)

  pop_codes <- populations$canonical_population_code
  pop_labels <- populations$wen_label %>% as.character()
  metric_mat <- matrix(NA_real_, nrow = length(pop_codes), ncol = length(pop_codes))
  rownames(metric_mat) <- pop_labels
  colnames(metric_mat) <- pop_labels

  for (i in seq_along(pop_codes)) {
    for (j in seq_along(pop_codes)) {
      if (i == j) {
        next
      }

      metrics <- calc_pair_metrics(
        sample_tbl = parsed_pg_all$sample_tbl,
        pop1 = pop_codes[i],
        pop2 = pop_codes[j],
        n_loci = parsed_pg_all$n_loci
      )

      metric_mat[i, j] <- if (i < j) metrics["fst"] else metrics["gst_prime"]
    }
  }

  as_tibble(metric_mat, rownames = "population")
}

tidy_pairwise_table <- function(tbl) {
  population_levels <- tbl$population

  tbl %>%
    pivot_longer(
      cols = -population,
      names_to = "population_2",
      values_to = "value"
    ) %>%
    mutate(
      population_1 = population,
      row_id = match(population_1, population_levels),
      col_id = match(population_2, population_levels),
      metric = case_when(
        row_id < col_id ~ "FST",
        row_id > col_id ~ "GST",
        TRUE ~ "Diagonal"
      )
    ) %>%
    filter(metric != "Diagonal", !is.na(value)) %>%
    transmute(
      metric,
      population_1 = factor(population_1, levels = rev(population_levels)),
      population_2 = factor(population_2, levels = population_levels),
      value
    )
}

build_heatmap <- function(tbl_long, fill_option = "inferno") {
  ggplot(tbl_long, aes(x = population_2, y = population_1, fill = value)) +
    geom_tile(color = "white", linewidth = 0) +
    scale_fill_viridis_c(
      option = fill_option,
      name = NULL
    ) +
    labs(
      title = "Upper diagonal: FST | Lower diagonal: GST",
      x = NULL,
      y = NULL
    ) +
    coord_fixed() +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

table_5 <- build_pairwise_table(read_pg_all(input_file))
table_5_long <- tidy_pairwise_table(table_5)
readr::write_csv(table_5, output_csv)

color_options <- c("inferno", "magma", "cividis", "plasma")

for (fill_option in color_options) {
  heatmap_plot <- build_heatmap(table_5_long, fill_option = fill_option)
  out_file <- if (fill_option == "inferno") {
    output_png
  } else {
    file.path(repo_root, "figure", paste0("Table_5_pairwise_FST_GST_heatmap_", fill_option, ".png"))
  }

  ggsave(out_file, plot = heatmap_plot, width = 3.7, height = 3.5, dpi = 600)
}
