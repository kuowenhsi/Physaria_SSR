read_wen_label_lookup <- function(repo_root) {
  label_file <- file.path(
    repo_root,
    "data",
    "metadata",
    "Physaria_globosa_population_counts_from_individuals.csv"
  )

  readr::read_csv(label_file, show_col_types = FALSE) %>%
    transmute(
      canonical_population_code = canonical_population_code,
      wen_label = Wen_Label,
      pg_population = paste0(
        "Pg",
        stringr::str_pad(readr::parse_number(canonical_population_code), width = 3, pad = "0")
      )
    )
}

apply_wen_label_from_pg <- function(x, lookup) {
  dplyr::coalesce(lookup$wen_label[match(x, lookup$pg_population)], x)
}

apply_wen_label_from_canonical <- function(x, lookup) {
  dplyr::coalesce(lookup$wen_label[match(x, lookup$canonical_population_code)], x)
}
