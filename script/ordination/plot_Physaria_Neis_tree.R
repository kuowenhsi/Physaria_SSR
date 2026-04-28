library(tidyverse)
library(ape)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(repo_root, "script", "utils", "utils_wen_labels.R"))

input_file <- file.path(repo_root, "data", "Physaria_Neis.txt")
output_dir <- file.path(repo_root, "figure", "Physaria_Neis_tree")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
wen_lookup <- read_wen_label_lookup(repo_root)
tip_key_file <- file.path(output_dir, "Physaria_Neis_tip_label_key.csv")

normalize_neis_label <- function(x) {
  case_when(
    x == "Pg348_KY_EX" ~ "Pgl348_EXS",
    x %in% c("Pg353_KY_EX", "Pg353_KY") ~ "Pgl353_EXS",
    x == "Pg355_KY_EX" ~ "Pgl355_EXS",
    x == "Pg356_TN_E" ~ "Pgl356_TNE",
    x == "Pg357_TN_E" ~ "Pgl357_TNE",
    x %in% c("Pg358_TN_C", "Pg358_TN_W") ~ "Pgl358_TNW",
    x %in% c("Pg359_TN_C", "Pg359_TN_W") ~ "Pgl359_TNW",
    x %in% c("Pg360_TN_C", "Pg360_TN_W") ~ "Pgl360_TNW",
    x == "Pg361_TN_W" ~ "Pgl361_TNW",
    x == "Pg362_TN_W" ~ "Pgl362_TNW",
    str_detect(x, "^Pg\\d+_") ~ str_replace(x, "^Pg", "Pgl"),
    x == "AR1" ~ "PoAR01",
    x == "AR2" ~ "PoAR02",
    x == "AR3" ~ "PoAR03",
    TRUE ~ x
  )
}

display_tree_label <- function(x) {
  canonical <- normalize_neis_label(x)
  wen_label <- apply_wen_label_from_canonical(canonical, wen_lookup)

  case_when(
    str_detect(canonical, "^Pgl") & canonical != wen_label ~ wen_label,
    str_detect(canonical, "^PoAR") ~ paste("P. ouachitensis", canonical),
    str_detect(canonical, "^MO") ~ paste("P. filiformis", canonical),
    str_detect(canonical, "^TX") ~ paste("P. gracilis", canonical),
    TRUE ~ canonical
  )
}

read_nei_matrix <- function(path) {
  lines <- readr::read_lines(path)
  matrix_start <- which(str_detect(lines, "^matrix\\s*$"))
  matrix_end <- which(
    seq_along(lines) > matrix_start[1] &
      str_detect(lines, "^\\s*;?\\s*end;\\s*$")
  )[1]

  if (length(matrix_start) == 0 || length(matrix_end) == 0) {
    stop("Could not find the distance matrix block in Physaria_Neis.txt.")
  }

  matrix_lines <- lines[(matrix_start + 1):(matrix_end - 1)]
  matrix_lines <- matrix_lines[str_squish(matrix_lines) != ""]

  parsed_rows <- purrr::map(matrix_lines, function(line) {
    pieces <- str_split(str_squish(line), "\\s+")[[1]]
    list(label = pieces[1], values = as.numeric(pieces[-1]))
  })

  taxa <- purrr::map_chr(parsed_rows, "label")
  n_taxa <- length(taxa)

  dist_mat <- matrix(NA_real_, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))

  for (i in seq_len(n_taxa)) {
    vals <- parsed_rows[[i]]$values
    dist_mat[i, seq_along(vals)] <- vals
  }

  dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
  diag(dist_mat) <- 0

  if (any(is.na(dist_mat))) {
    stop("The Nei matrix still contains missing values after parsing.")
  }

  dist_mat
}

dist_mat <- read_nei_matrix(input_file)
nei_dist <- as.dist(dist_mat)

nj_tree <- ape::nj(nei_dist)
nj_tree <- ape::root(nj_tree, outgroup = "TX1", resolve.root = TRUE)
nj_tree <- ape::ladderize(nj_tree, right = TRUE)

negative_edges <- which(nj_tree$edge.length < 0)
if (length(negative_edges) > 0) {
  nj_tree$edge.length[negative_edges] <- 0
}

label_key <- tibble(
  original_label = rownames(dist_mat),
  canonical_label = normalize_neis_label(original_label),
  species = case_when(
    str_detect(canonical_label, "^Pgl") ~ "P. globosa",
    str_detect(canonical_label, "^PoAR") ~ "P. ouachitensis",
    str_detect(canonical_label, "^MO") ~ "P. filiformis",
    str_detect(canonical_label, "^TX") ~ "P. gracilis",
    TRUE ~ NA_character_
  ),
  plot_label = display_tree_label(original_label)
) %>%
  distinct()

if (file.exists(tip_key_file)) {
  existing_tip_key <- readr::read_csv(tip_key_file, show_col_types = FALSE) %>%
    select(any_of(c("original_label", "Wen_plot_label")))

  label_key <- label_key %>%
    left_join(existing_tip_key, by = "original_label") %>%
    mutate(final_plot_label = coalesce(Wen_plot_label, plot_label))
} else {
  label_key <- label_key %>%
    mutate(
      Wen_plot_label = plot_label,
      final_plot_label = plot_label
    )
}

nj_tree$tip.label <- label_key$final_plot_label[match(nj_tree$tip.label, label_key$original_label)]

png(
  filename = file.path(output_dir, "Physaria_Neis_neighbor_joining_tree.png"),
  width = 1800,
  height = 1200,
  res = 300
)
par(mar = c(0.6, 0.6, 1.2, 6.2))
plot(
  nj_tree,
  type = "phylogram",
  cex = 0.56,
  font = 1,
  no.margin = FALSE,
  direction = "downwards",
  label.offset = 0.06
)
add.scale.bar(length = 0.2, lwd = 1)
dev.off()

write.tree(nj_tree, file = file.path(output_dir, "Physaria_Neis_neighbor_joining_tree.nwk"))
readr::write_csv(
  label_key %>% select(original_label, canonical_label, species, plot_label, Wen_plot_label, final_plot_label),
  tip_key_file
)
readr::write_csv(
  tibble(
    edge_index = seq_along(nj_tree$edge.length),
    parent_node = nj_tree$edge[, 1],
    child_node = nj_tree$edge[, 2],
    child_label = if_else(
      nj_tree$edge[, 2] <= length(nj_tree$tip.label),
      nj_tree$tip.label[nj_tree$edge[, 2]],
      NA_character_
    ),
    edge_length = nj_tree$edge.length
  ),
  file.path(output_dir, "Physaria_Neis_tree_edge_lengths.csv")
)
readr::write_csv(
  as_tibble(dist_mat, rownames = "taxon"),
  file.path(output_dir, "Physaria_Neis_matrix_parsed.csv")
)
