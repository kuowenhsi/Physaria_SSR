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

metadata_dir <- file.path(repo_root, "data", "metadata")
old_counts_file <- file.path(metadata_dir, "Physaria_globosa_population_counts_from_individuals.csv")
old_individual_file <- file.path(metadata_dir, "Physaria_globosa_individual_metadata_clean.csv")
old_crosswalk_file <- file.path(metadata_dir, "Physaria_population_name_crosswalk.csv")
old_table_4_file <- file.path(repo_root, "data", "Table_4.csv")

new_counts_file <- file.path(metadata_dir, "Physaria_globosa_population_counts_from_individuals_relabeled.csv")
new_individual_file <- file.path(metadata_dir, "Physaria_globosa_individual_metadata_clean_relabeled.csv")
new_crosswalk_file <- file.path(metadata_dir, "Physaria_population_name_crosswalk_relabeled.csv")
new_table_4_file <- file.path(repo_root, "data", "Table_4_relabeled.csv")

relabel_wen_code <- function(x) {
  recode(
    x,
    "EX_01" = "TR_01",
    "EX_02" = "TR_02",
    "EX_03" = "EX_01",
    .default = x
  )
}

counts <- readr::read_csv(old_counts_file, show_col_types = FALSE) %>%
  mutate(
    Wen_Label_original = Wen_Label,
    Wen_Label = relabel_wen_code(Wen_Label)
  )

readr::write_csv(counts, new_counts_file)

individual_metadata <- readr::read_csv(old_individual_file, show_col_types = FALSE) %>%
  left_join(
    counts %>%
      select(canonical_population_code, Wen_Label, Wen_Label_original),
    by = "canonical_population_code"
  ) %>%
  relocate(Wen_Label, Wen_Label_original, .after = canonical_population_code)

readr::write_csv(individual_metadata, new_individual_file)

population_crosswalk <- readr::read_csv(old_crosswalk_file, show_col_types = FALSE) %>%
  left_join(
    counts %>%
      select(canonical_population_code, Wen_Label, Wen_Label_original),
    by = c("table1_population_code" = "canonical_population_code")
  ) %>%
  relocate(Wen_Label, Wen_Label_original, .after = table1_population_code)

readr::write_csv(population_crosswalk, new_crosswalk_file)

raw_lines <- readr::read_lines(old_table_4_file)
comment_lines <- raw_lines[str_detect(raw_lines, '^"?#')]
data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

table_4 <- readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
  mutate(
    Wen_label_original = Wen_label,
    Wen_label = relabel_wen_code(Wen_label)
  ) %>%
  relocate(Wen_label_original, .after = Wen_label)

if (length(comment_lines) > 0) {
  tmp_data_file <- tempfile(fileext = ".csv")
  readr::write_csv(table_4, tmp_data_file)
  readr::write_lines(c(comment_lines, "", readr::read_lines(tmp_data_file)), new_table_4_file)
} else {
  readr::write_csv(table_4, new_table_4_file)
}

message("Wrote relabeled metadata duplicates:")
message("  ", new_counts_file)
message("  ", new_individual_file)
message("  ", new_crosswalk_file)
message("  ", new_table_4_file)
