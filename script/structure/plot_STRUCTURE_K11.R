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
source(file.path(repo_root, "script", "utils", "utils_wen_labels.R"))

input_file <- file.path(repo_root, "data", "coordsandadmixprpo.pglo.xlsx")
output_png <- file.path(repo_root, "figure", "STRUCTURE_K11_barplot.png")
wen_lookup <- read_wen_label_lookup(repo_root)

read_structure_k11 <- function(path) {
  readxl::read_excel(path, sheet = 1) %>%
    rename(
      voucher_code = `Voucher number/ code`,
      latitude = Lat,
      longitude = long
    ) %>%
    mutate(
      voucher_code = str_replace_all(voucher_code, fixed("*"), ""),
      numeric_code = readr::parse_number(voucher_code),
      population_label = case_when(
        str_detect(voucher_code, "^AR") ~ paste0("Po", voucher_code),
        str_detect(voucher_code, "^MO") ~ paste0("Pf", voucher_code),
        str_detect(voucher_code, "^TX") ~ paste0("Pgr", voucher_code),
        TRUE ~ paste0("Pgl", voucher_code)
      ),
      species = case_when(
        str_detect(voucher_code, "^AR") ~ "P. ouachitensis",
        str_detect(voucher_code, "^MO") ~ "P. filiformis",
        str_detect(voucher_code, "^TX") ~ "P. gracilis",
        TRUE ~ "P. globosa"
      ),
      group = case_when(
        voucher_code == "346" ~ "IN",
        voucher_code %in% c("347", "349", "350", "351", "352", "354") ~ "KY",
        voucher_code %in% c("348", "353", "355") ~ "EX-S",
        voucher_code %in% c("356", "357") ~ "TNE",
        voucher_code %in% c("358", "359", "360", "361", "362") ~ "TNW",
        str_detect(voucher_code, "^AR") ~ "AR",
        str_detect(voucher_code, "^MO") ~ "MO",
        str_detect(voucher_code, "^TX") ~ "TX",
        TRUE ~ "Other"
      )
    ) %>%
    mutate(
      canonical_population_code = case_when(
        species == "P. globosa" ~ paste0("Pgl", voucher_code, "_", recode(group, "EX-S" = "EXS", .default = group)),
        TRUE ~ population_label
      ),
      population_label = apply_wen_label_from_canonical(canonical_population_code, wen_lookup)
    ) %>%
    mutate(
      order_id = row_number(),
      across(`1`:`11`, as.numeric)
    )
}

build_structure_plot <- function(tbl) {
  cluster_palette <- c(
    "1" = "#1b9e77",
    "2" = "#d95f02",
    "3" = "#7570b3",
    "4" = "#e7298a",
    "5" = "#66a61e",
    "6" = "#e6ab02",
    "7" = "#a6761d",
    "8" = "#1f78b4",
    "9" = "#b2df8a",
    "10" = "#fb9a99",
    "11" = "#cab2d6"
  )

  long_tbl <- tbl %>%
    pivot_longer(
      cols = `1`:`11`,
      names_to = "cluster",
      values_to = "membership"
    ) %>%
    mutate(
      population_label = factor(population_label, levels = tbl$population_label),
      cluster = factor(cluster, levels = as.character(1:11))
    )

  species_breaks <- tbl %>%
    distinct(species, order_id) %>%
    group_by(species) %>%
    summarise(
      xmin = min(order_id) - 0.5,
      xmax = max(order_id) + 0.5,
      xmid = mean(c(xmin, xmax)),
      .groups = "drop"
    )

  separator_x <- head(species_breaks$xmax, -1)

  ggplot(long_tbl, aes(x = population_label, y = membership, fill = cluster)) +
    geom_col(width = 0.92, color = NA) +
    geom_vline(xintercept = separator_x, color = "gray35", linewidth = 0.45) +
    annotate(
      "text",
      x = species_breaks$xmid,
      y = 1.04,
      label = species_breaks$species,
      size = 3.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = cluster_palette, name = "Cluster") +
    scale_y_continuous(
      limits = c(0, 1.08),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    labs(
      title = "STRUCTURE ancestry coefficients at K = 11",
      subtitle = "Population-level membership from coordsandadmixprpo.pglo.xlsx",
      x = NULL,
      y = "Ancestry coefficient"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

structure_tbl <- read_structure_k11(input_file)
p <- build_structure_plot(structure_tbl)

ggsave(output_png, plot = p, width = 12, height = 6.5, dpi = 600)
