library(tidyverse)
library(readxl)

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

pg_all_file <- file.path(repo_root, "data", "PG_All.csv")
pg_all_structure_file <- file.path(repo_root, "data", "PG_All_structure.txt")
table_1_file <- file.path(repo_root, "data", "Table_1.csv")
table_3_file <- file.path(repo_root, "data", "Table_3.csv")
coords_file <- file.path(repo_root, "data", "coordsandadmixprpo.pglo.xlsx")
pglo_coords_file <- file.path(repo_root, "data", "pglo.coordinates.txt")
pca_metadata_file <- file.path(repo_root, "figure", "PG_All_ordination", "PG_All_metadata_clean.csv")
structure_k11_file <- file.path(repo_root, "figure", "STRUCTURE_individual_allK", "STRUCTURE_K11_individual_mean.csv")
output_dir <- file.path(repo_root, "data", "metadata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

state_from_pg_population <- function(population) {
  case_when(
    population == "Pg346" ~ "IN",
    population %in% c("Pg347", "Pg349", "Pg350", "Pg351", "Pg352", "Pg354") ~ "KY",
    population %in% c("Pg348", "Pg353", "Pg355") ~ "EX-S",
    population %in% c("Pg356", "Pg357") ~ "TNE",
    population %in% c("Pg358", "Pg359", "Pg360", "Pg361", "Pg362") ~ "TNW",
    TRUE ~ "Unknown"
  )
}

code_suffix_from_state <- function(state_region) {
  recode(
    state_region,
    "IN" = "IN",
    "KY" = "KY",
    "EX-S" = "EXS",
    "TNE" = "TNE",
    "TNW" = "TNW",
    .default = NA_character_
  )
}

clean_table_1 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
    rename(
      table1_population_code_raw = `Population code`,
      collector = Collector,
      species = `Species name`,
      locality = `State, County, Locality (Element Occurrence number)`,
      ecoregion = Ecoregion,
      lat_raw = Lat,
      long_raw = Long,
      n_samples_table1 = `Number of samples`
    ) %>%
    mutate(
      species = str_replace_all(species, "[“”]", ""),
      species = str_squish(species),
      table1_population_code = str_remove(table1_population_code_raw, "\\*$"),
      population_numeric_code = parse_number(table1_population_code),
      state_region = case_when(
        str_detect(table1_population_code, "^Pgl") ~ str_extract(table1_population_code, "(?<=_)[A-Z]+$"),
        table1_population_code == "PoAR01" ~ "AR",
        table1_population_code == "PoAR02" ~ "AR",
        table1_population_code == "PoAR03" ~ "AR",
        table1_population_code == "PfMO6" ~ "MO",
        table1_population_code == "PfMO10" ~ "MO",
        table1_population_code == "PfMO11" ~ "MO",
        table1_population_code == "PgrTX1" ~ "TX",
        TRUE ~ NA_character_
      ),
      voucher_code_clean = case_when(
        str_detect(table1_population_code, "^Pgl") ~ as.character(population_numeric_code),
        table1_population_code == "PoAR01" ~ "AR01",
        table1_population_code == "PoAR02" ~ "AR02",
        table1_population_code == "PoAR03" ~ "AR03",
        table1_population_code == "PfMO6" ~ "MO6",
        table1_population_code == "PfMO10" ~ "MO10",
        table1_population_code == "PfMO11" ~ "MO11",
        table1_population_code == "PgrTX1" ~ "TX1",
        TRUE ~ NA_character_
      ),
      lat_table1 = readr::parse_number(as.character(lat_raw)),
      long_table1 = readr::parse_number(as.character(long_raw)),
      n_samples_table1 = readr::parse_number(as.character(n_samples_table1))
    ) %>%
    select(
      table1_population_code_raw,
      table1_population_code,
      voucher_code_clean,
      population_numeric_code,
      species,
      state_region,
      collector,
      locality,
      ecoregion,
      lat_table1,
      long_table1,
      n_samples_table1
    )
}

clean_table_3 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), na = c("", "NA"), show_col_types = FALSE) %>%
    mutate(
      across(c(n, A, AR, AP, Ho, He, F, `% null`, FB), as.numeric),
      row_type = case_when(
        str_detect(Population, "^Average across") ~ "average",
        if_all(all_of(c("n", "A", "AR", "AP", "Ho", "He", "F", "% null", "FB")), is.na) ~ "section",
        TRUE ~ "population"
      ),
      species_header = if_else(row_type == "section", Population, NA_character_)
    ) %>%
    fill(species_header, .direction = "down") %>%
    mutate(
      species = case_when(
        str_detect(Population, "^Average across P\\. globosa") ~ "P. globosa",
        str_detect(Population, "^Average across P\\. ouachitensis") ~ "P. ouachitensis",
        str_detect(Population, "^Average across P\\. filiformis") ~ "P. filiformis",
        str_detect(Population, "^Average across P\\. gracilis") ~ "P. gracilis",
        TRUE ~ species_header
      )
    ) %>%
    filter(row_type == "population") %>%
    transmute(
      table3_population_raw = Population,
      species,
      geographic_area = `Geographic area`,
      n_samples_table3 = n,
      table3_population_code = case_when(
        str_detect(table3_population_raw, "^Pgl") ~ table3_population_raw,
        table3_population_raw == "PoAR1" ~ "PoAR01",
        table3_population_raw == "PoAR2" ~ "PoAR02",
        table3_population_raw == "PoAR3" ~ "PoAR03",
        TRUE ~ table3_population_raw
      ),
      voucher_code_clean = case_when(
        str_detect(table3_population_raw, "^Pgl") ~ as.character(parse_number(table3_population_raw)),
        table3_population_raw == "PoAR1" ~ "AR01",
        table3_population_raw == "PoAR2" ~ "AR02",
        table3_population_raw == "PoAR3" ~ "AR03",
        table3_population_raw == "PfMO6" ~ "MO6",
        table3_population_raw == "PfMO10" ~ "MO10",
        table3_population_raw == "PfMO11" ~ "MO11",
        table3_population_raw == "PgrTX1" ~ "TX1",
        TRUE ~ NA_character_
      ),
      table3_has_zero_padding_mismatch = table3_population_raw %in% c("PoAR1", "PoAR2", "PoAR3")
    )
}

read_pg_all <- function(path) {
  raw <- readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)

  header_population_labels <- raw[2, -(1:3)] %>%
    unlist(use.names = FALSE) %>%
    as.character()
  header_population_labels <- header_population_labels[header_population_labels != ""]

  sample_tbl <- raw[-c(1:3), ] %>%
    setNames(c("sample_pg_all_raw", "site_pg_all_raw", paste0("allele_", seq_len(ncol(raw) - 2)))) %>%
    transmute(
      sample_pg_all_raw = as.character(sample_pg_all_raw),
      site_pg_all_raw = as.character(site_pg_all_raw),
      sample_pg_all_clean = sample_pg_all_raw,
      population = str_extract(sample_pg_all_raw, "^Pg\\d+"),
      sample_index_within_population = parse_number(str_extract(sample_pg_all_raw, "\\d+$")),
      site_pg_all_clean = str_remove(site_pg_all_raw, "^p"),
      site_matches_population = site_pg_all_clean == str_remove(population, "^Pg"),
      pg_all_site_issue = case_when(
        str_detect(site_pg_all_raw, "^p") ~ "Site code used a p-prefix and was normalized",
        !site_matches_population ~ "Site code differs from population code",
        TRUE ~ NA_character_
      )
    )

  population_order <- sample_tbl %>%
    distinct(population) %>%
    pull(population)

  header_tbl <- tibble(
    population = population_order,
    pg_all_header_site_label_raw = header_population_labels[seq_along(population_order)]
  ) %>%
    mutate(
      pg_all_header_site_label_clean = str_remove(pg_all_header_site_label_raw, "^p"),
      pg_all_header_matches_population = pg_all_header_site_label_clean == str_remove(population, "^Pg"),
      pg_all_header_issue = case_when(
        !pg_all_header_matches_population ~ "Header label differs from population code",
        TRUE ~ NA_character_
      )
    )

  sample_tbl %>%
    left_join(header_tbl, by = "population")
}

read_structure_samples <- function(path) {
  lines <- readr::read_lines(path)
  sample_lines <- lines[-1]
  sample_lines <- sample_lines[str_squish(sample_lines) != ""]

  tibble(raw_line = sample_lines) %>%
    mutate(
      parts = str_split(str_squish(raw_line), "\\s+"),
      sample_id_structure_raw = map_chr(parts, 1),
      structure_population_index = map_chr(parts, ~ .x[2]) %>% as.integer(),
      sample_order = row_number(),
      sample_pg_clean = str_remove(sample_id_structure_raw, "^\\d+_"),
      sample_pg_clean = str_remove(sample_pg_clean, "\\.fsa$"),
      population = str_extract(sample_pg_clean, "^Pg\\d+"),
      sample_index_within_population = parse_number(str_extract(sample_pg_clean, "\\d+$")),
      voucher_code_clean = str_remove(population, "^Pg")
    ) %>%
    select(
      sample_order,
      sample_id_structure_raw,
      sample_pg_clean,
      population,
      voucher_code_clean,
      sample_index_within_population,
      structure_population_index
    )
}

read_coords_workbook <- function(path) {
  readxl::read_excel(path) %>%
    transmute(
      workbook_voucher_code_raw = `Voucher number/ code`,
      workbook_voucher_code = str_remove(workbook_voucher_code_raw, "\\*$"),
      lat_workbook = readr::parse_number(Lat),
      long_workbook = as.numeric(long),
      workbook_population_code = case_when(
        str_detect(workbook_voucher_code, "^\\d+$") ~ paste0(
          "Pgl",
          workbook_voucher_code,
          "_",
          code_suffix_from_state(state_from_pg_population(paste0("Pg", workbook_voucher_code)))
        ),
        workbook_voucher_code == "AR01" ~ "PoAR01",
        workbook_voucher_code == "AR02" ~ "PoAR02",
        workbook_voucher_code == "AR03" ~ "PoAR03",
        workbook_voucher_code == "MO6" ~ "PfMO6",
        workbook_voucher_code == "MO10" ~ "PfMO10",
        workbook_voucher_code == "MO11" ~ "PfMO11",
        workbook_voucher_code == "TX1" ~ "PgrTX1",
        TRUE ~ NA_character_
      )
    )
}

read_pglo_coordinates <- function(path) {
  readr::read_tsv(
    path,
    col_names = c("voucher_code_clean", "lat_pglo_coordinates", "long_pglo_coordinates"),
    show_col_types = FALSE
  ) %>%
    mutate(voucher_code_clean = as.character(voucher_code_clean))
}

table_1 <- clean_table_1(table_1_file)
table_3 <- clean_table_3(table_3_file)
pg_all <- read_pg_all(pg_all_file)
structure_samples <- read_structure_samples(pg_all_structure_file)
coords_workbook <- read_coords_workbook(coords_file)
pglo_coords <- read_pglo_coordinates(pglo_coords_file)
pca_metadata <- readr::read_csv(pca_metadata_file, show_col_types = FALSE) %>%
  rename(
    sample_pg_clean = sample,
    state_from_pca_dapc = state
  ) %>%
  select(sample_pg_clean, state_from_pca_dapc)
structure_k11 <- readr::read_csv(structure_k11_file, show_col_types = FALSE) %>%
  rename(
    sample_id_structure_raw = sample_id,
    state_from_structure_plot = state
  ) %>%
  select(sample_id_structure_raw, state_from_structure_plot, starts_with("cluster_"))

pg_globosa_population <- table_1 %>%
  filter(species == "P. globosa") %>%
  left_join(
    table_3 %>% select(table3_population_raw, table3_population_code, voucher_code_clean, geographic_area, n_samples_table3),
    by = c("table1_population_code" = "table3_population_code", "voucher_code_clean")
  ) %>%
  left_join(coords_workbook, by = c("voucher_code_clean" = "workbook_voucher_code")) %>%
  left_join(pglo_coords, by = "voucher_code_clean") %>%
  mutate(
    population = paste0("Pg", str_pad(voucher_code_clean, width = 3, pad = "0")),
    state_region = recode(state_region, "EXS" = "EX-S")
  )

individual_metadata <- structure_samples %>%
  left_join(
    pg_all,
    by = c(
      "sample_pg_clean" = "sample_pg_all_clean",
      "population",
      "sample_index_within_population"
    )
  ) %>%
  left_join(pg_globosa_population, by = c("population", "voucher_code_clean")) %>%
  left_join(pca_metadata, by = "sample_pg_clean") %>%
  left_join(structure_k11, by = "sample_id_structure_raw") %>%
  mutate(
    sample_id_consistent = sample_pg_clean,
    individual_name_is_consistent = sample_pg_clean == sample_pg_all_raw,
    canonical_population_code = paste0(
      "Pgl",
      voucher_code_clean,
      "_",
      code_suffix_from_state(state_region)
    ),
    figure_group_state = state_region,
    population_name_inconsistency = case_when(
      !pg_all_header_matches_population ~ "PG_All header uses 253 where sample/population labels use 353",
      TRUE ~ NA_character_
    ),
    individual_name_notes = case_when(
      !individual_name_is_consistent ~ "PG_All and STRUCTURE sample labels differ",
      !is.na(pg_all_site_issue) ~ pg_all_site_issue,
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    sample_order,
    sample_id_structure_raw,
    sample_pg_all_raw,
    sample_id_consistent,
    sample_pg_clean,
    individual_name_is_consistent,
    sample_index_within_population,
    population,
    voucher_code_clean,
    structure_population_index,
    site_pg_all_raw,
    site_pg_all_clean,
    site_matches_population,
    pg_all_site_issue,
    pg_all_header_site_label_raw,
    pg_all_header_site_label_clean,
    pg_all_header_matches_population,
    pg_all_header_issue,
    canonical_population_code,
    table1_population_code_raw,
    table1_population_code,
    table3_population_raw,
    species,
    state_region,
    geographic_area,
    figure_group_state,
    collector,
    locality,
    ecoregion,
    lat_table1,
    long_table1,
    lat_pglo_coordinates,
    long_pglo_coordinates,
    lat_workbook,
    long_workbook,
    n_samples_table1,
    n_samples_table3,
    state_from_pca_dapc,
    state_from_structure_plot,
    starts_with("cluster_"),
    population_name_inconsistency,
    individual_name_notes
  ) %>%
  arrange(sample_order)

stopifnot(nrow(individual_metadata) == 352)

population_crosswalk <- table_1 %>%
  left_join(table_3, by = c("table1_population_code" = "table3_population_code", "voucher_code_clean", "species")) %>%
  left_join(coords_workbook, by = c("voucher_code_clean" = "workbook_voucher_code")) %>%
  left_join(pglo_coords, by = "voucher_code_clean") %>%
  mutate(
    pg_population_label = case_when(
      species == "P. globosa" ~ paste0("Pg", str_pad(voucher_code_clean, width = 3, pad = "0")),
      TRUE ~ NA_character_
    ),
    state_region = recode(state_region, "EXS" = "EX-S"),
    alias_notes = case_when(
      str_detect(table1_population_code_raw, "\\*$") ~ "Table 1 includes an asterisk suffix for previously published outgroup samples",
      table3_has_zero_padding_mismatch ~ "Table 3 drops leading zero in P. ouachitensis population labels",
      species == "P. globosa" & voucher_code_clean == "353" ~ "PG_All.csv header row shows 253, while sample labels and other tables use 353",
      TRUE ~ NA_character_
    ),
    population_aliases = pmap_chr(
      list(
        table1_population_code_raw,
        table1_population_code,
        table3_population_raw,
        workbook_voucher_code_raw,
        pg_population_label
      ),
      function(...) {
        vals <- c(...)
        vals <- vals[!is.na(vals) & vals != ""]
        paste(unique(vals), collapse = " | ")
      }
    )
  ) %>%
  select(
    species,
    state_region,
    voucher_code_clean,
    population_numeric_code,
    table1_population_code_raw,
    table1_population_code,
    table3_population_raw,
    geographic_area,
    workbook_voucher_code_raw,
    workbook_population_code,
    pg_population_label,
    collector,
    locality,
    ecoregion,
    lat_table1,
    long_table1,
    lat_pglo_coordinates,
    long_pglo_coordinates,
    lat_workbook,
    long_workbook,
    n_samples_table1,
    n_samples_table3,
    population_aliases,
    alias_notes
  ) %>%
  arrange(species, state_region, voucher_code_clean)

readr::write_csv(
  individual_metadata,
  file.path(output_dir, "Physaria_globosa_individual_metadata_clean.csv")
)

readr::write_csv(
  population_crosswalk,
  file.path(output_dir, "Physaria_population_name_crosswalk.csv")
)

readr::write_csv(
  individual_metadata %>%
    count(
      state_region,
      canonical_population_code,
      name = "n_individuals"
    ),
  file.path(output_dir, "Physaria_globosa_population_counts_from_individuals.csv")
)
