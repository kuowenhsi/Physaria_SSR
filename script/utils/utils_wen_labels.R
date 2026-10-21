read_wen_label_lookup <- function(repo_root, label_file = NULL) {
  env_label_file <- Sys.getenv("PHYSARIA_WEN_LABEL_FILE", unset = "")
  if (is.null(label_file) && nzchar(env_label_file)) {
    label_file <- env_label_file
  }
  if (is.null(label_file)) {
    label_file <- file.path(
      repo_root,
      "data",
      "metadata",
      "Physaria_globosa_population_counts_from_individuals.csv"
    )
  }

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

use_wen_groups <- function() {
  env_value <- Sys.getenv("PHYSARIA_GROUP_FROM_WEN_LABEL", unset = "")
  if (nzchar(env_value)) {
    return(tolower(env_value) %in% c("1", "true", "yes"))
  }

  nzchar(physaria_output_suffix())
}

wen_group_from_label <- function(label) {
  stringr::str_extract(label, "^[A-Za-z]+")
}

apply_wen_group_from_pg <- function(x, fallback, lookup) {
  if (!use_wen_groups()) {
    return(fallback)
  }

  label <- apply_wen_label_from_pg(x, lookup)
  dplyr::coalesce(wen_group_from_label(label), fallback)
}

apply_wen_group_from_canonical <- function(x, fallback, lookup) {
  if (!use_wen_groups()) {
    return(fallback)
  }

  label <- apply_wen_label_from_canonical(x, lookup)
  dplyr::coalesce(wen_group_from_label(label), fallback)
}

physaria_output_suffix <- function() {
  Sys.getenv("PHYSARIA_OUTPUT_SUFFIX", unset = "")
}

suffix_filename <- function(filename, suffix = physaria_output_suffix()) {
  if (!nzchar(suffix)) {
    return(filename)
  }

  ext <- tools::file_ext(filename)
  stem <- if (nzchar(ext)) {
    sub(paste0("\\.", ext, "$"), "", filename)
  } else {
    filename
  }

  if (nzchar(ext)) {
    paste0(stem, suffix, ".", ext)
  } else {
    paste0(filename, suffix)
  }
}
