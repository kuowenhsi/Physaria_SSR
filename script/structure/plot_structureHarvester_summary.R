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

summary_file <- file.path(repo_root, "results", "structureHarvester", "summary.txt")
evanno_file <- file.path(repo_root, "results", "structureHarvester", "evanno.txt")
output_dir <- file.path(repo_root, "figure", "structureHarvester")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

extract_table_block <- function(path, header_pattern) {
  lines <- readr::read_lines(path)
  header_idx <- grep(header_pattern, lines)
  if (length(header_idx) == 0) {
    stop("Could not find table header matching: ", header_pattern)
  }
  start_idx <- header_idx[1] + 1
  end_candidates <- which(seq_along(lines) > start_idx & (lines == "" | str_detect(lines, "^##########") | str_detect(lines, "^Software written")))
  end_idx <- if (length(end_candidates) > 0) min(end_candidates) - 1 else length(lines)
  lines[start_idx:end_idx]
}

read_summary_k_table <- function(path) {
  block <- extract_table_block(path, "^# K\\tReps\\tmean est\\. LnP\\(Data\\)")
  readr::read_tsv(
    I(block),
    col_names = c("K", "Reps", "mean_lnP", "sd_lnP"),
    show_col_types = FALSE
  )
}

read_summary_run_table <- function(path) {
  block <- extract_table_block(path, "^# File name\\tRun #\\tK\\tEst\\. Ln prob\\. of data")
  readr::read_tsv(
    I(block),
    col_names = c("file_name", "run_number", "K", "est_lnP", "mean_llh", "var_llh"),
    show_col_types = FALSE
  ) %>%
    mutate(run_number = as.integer(run_number))
}

read_evanno_table <- function(path) {
  block <- extract_table_block(path, "^# K\\tReps\\tMean LnP\\(K\\)")
  readr::read_tsv(
    I(block),
    col_names = c("K", "Reps", "Mean_LnP_K", "Stdev_LnP_K", "Ln_prime", "Abs_Ln_double_prime", "Delta_K"),
    na = c("NA"),
    show_col_types = FALSE
  )
}

save_plot <- function(plot_obj, filename, width = 9, height = 6) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

summary_k <- read_summary_k_table(summary_file)
summary_runs <- read_summary_run_table(summary_file)
evanno <- read_evanno_table(evanno_file)

best_delta_k <- evanno %>%
  filter(!is.na(Delta_K)) %>%
  slice_max(order_by = Delta_K, n = 1, with_ties = FALSE)

best_mean_lnp <- summary_k %>%
  slice_max(order_by = mean_lnP, n = 1, with_ties = FALSE)

p_mean_lnp <- ggplot(summary_k, aes(x = K, y = mean_lnP)) +
  geom_line(color = "#2c7fb8", linewidth = 0.7) +
  geom_point(color = "#2c7fb8", size = 2.2) +
  geom_errorbar(aes(ymin = mean_lnP - sd_lnP, ymax = mean_lnP + sd_lnP), width = 0.15, color = "#2c7fb8") +
  geom_point(
    data = best_mean_lnp,
    color = "#d95f02",
    size = 3
  ) +
  labs(
    title = "Mean estimated LnP(Data) across K",
    subtitle = "Error bars show standard deviation across 5 replicates",
    x = "K",
    y = "Mean estimated LnP(Data)"
  ) +
  theme_minimal(base_size = 11)

p_run_lnp <- ggplot(summary_runs, aes(x = K, y = est_lnP)) +
  geom_jitter(width = 0.14, height = 0, size = 2, alpha = 0.75, color = "#7570b3") +
  geom_line(
    data = summary_k,
    aes(x = K, y = mean_lnP),
    inherit.aes = FALSE,
    color = "#1b9e77",
    linewidth = 0.7
  ) +
  geom_point(
    data = summary_k,
    aes(x = K, y = mean_lnP),
    inherit.aes = FALSE,
    color = "#1b9e77",
    size = 2.2
  ) +
  labs(
    title = "Replicate-specific estimated LnP(Data)",
    subtitle = "Purple points are individual runs; green line is the mean by K",
    x = "K",
    y = "Estimated LnP(Data)"
  ) +
  theme_minimal(base_size = 11)

p_delta_k <- ggplot(evanno, aes(x = K, y = Delta_K)) +
  geom_col(fill = "#d95f02", alpha = 0.9) +
  geom_point(data = best_delta_k, color = "black", size = 3) +
  labs(
    title = "Evanno Delta K profile",
    subtitle = "Highest Delta K is highlighted",
    x = "K",
    y = "Delta K"
  ) +
  theme_minimal(base_size = 11)

p_sd_lnp <- ggplot(summary_k, aes(x = K, y = sd_lnP)) +
  geom_col(fill = "#636363", alpha = 0.9) +
  labs(
    title = "Run-to-run variability by K",
    subtitle = "Standard deviation of estimated LnP(Data) across replicates",
    x = "K",
    y = "SD of estimated LnP(Data)"
  ) +
  theme_minimal(base_size = 11)

save_plot(p_mean_lnp, "structureHarvester_mean_lnP.png")
save_plot(p_run_lnp, "structureHarvester_run_lnP.png")
save_plot(p_delta_k, "structureHarvester_deltaK.png")
save_plot(p_sd_lnp, "structureHarvester_sd_lnP.png")

readr::write_csv(summary_k, file.path(output_dir, "structureHarvester_summary_by_K.csv"))
readr::write_csv(summary_runs, file.path(output_dir, "structureHarvester_summary_by_run.csv"))
readr::write_csv(evanno, file.path(output_dir, "structureHarvester_evanno.csv"))
