library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(scatterpie)

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

input_file <- file.path(repo_root, "data", "Table_1.csv")
structure_dir <- file.path(repo_root, "figure", "STRUCTURE_individual_allK")
largewater_file <- file.path(repo_root, "data", "P_globosa_map", "P_glob_largewater.shp")
smallwater_file <- file.path(repo_root, "data", "P_globosa_map", "P_glob_smallwater.shp")
wen_lookup <- read_wen_label_lookup(repo_root)
k_values <- 2:8

read_table_1 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
    rename(
      population_code = `Population code`,
      species_name = `Species name`,
      ecoregion = Ecoregion,
      lat_raw = Lat,
      long_raw = Long,
      n_samples = `Number of samples`
    ) %>%
    mutate(
      species_name = str_replace_all(species_name, "[“”]", ""),
      species_name = str_squish(species_name),
      lat = readr::parse_number(as.character(lat_raw)),
      long = readr::parse_number(as.character(long_raw)),
      n_samples = readr::parse_number(as.character(n_samples))
    ) %>%
    filter(!is.na(lat), !is.na(long))
}

read_structure_k <- function(path, k) {
  cluster_cols <- paste0("cluster_", seq_len(k))

  readr::read_csv(path, show_col_types = FALSE) %>%
    group_by(population_label) %>%
    summarise(across(all_of(cluster_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

compute_map_extent <- function(tbl) {
  tbl <- tbl %>% filter(species_name == "P. globosa")
  site_sf <- st_as_sf(tbl, coords = c("long", "lat"), crs = 4326, remove = FALSE)
  site_sf_3857 <- st_transform(site_sf, 3857)
  bbox <- st_bbox(site_sf_3857)
  x_span <- bbox$xmax - bbox$xmin
  y_span <- bbox$ymax - bbox$ymin
  x_pad <- x_span * 0.10 + x_span * 0.14
  y_pad <- y_span * 0.10 + y_span * 0.16
  st_bbox(c(
    xmin = as.numeric(bbox["xmin"]) - as.numeric(x_pad),
    ymin = as.numeric(bbox["ymin"]) - as.numeric(y_pad),
    xmax = as.numeric(bbox["xmax"]) + as.numeric(x_pad),
    ymax = as.numeric(bbox["ymax"]) + as.numeric(y_pad)
  ), crs = st_crs(site_sf_3857))
}

generate_random_offsets <- function(point_tbl, map_bbox, seed = 42, frac = 0.20) {
  set.seed(seed)
  x_limit <- (map_bbox["xmax"] - map_bbox["xmin"]) * frac
  y_limit <- (map_bbox["ymax"] - map_bbox["ymin"]) * frac

  point_tbl %>%
    st_drop_geometry() %>%
    transmute(
      point_id,
      dx = round(runif(n(), -x_limit, x_limit), 1),
      dy = round(runif(n(), -y_limit, y_limit), 1)
    )
}

tbl <- read_table_1(input_file)
map_bbox <- compute_map_extent(tbl)
map_bbox_sfc <- st_as_sfc(map_bbox)

globosa_tbl <- tbl %>%
  filter(species_name == "P. globosa") %>%
  mutate(point_id = apply_wen_label_from_canonical(population_code, wen_lookup))

site_sf_3857 <- st_as_sf(globosa_tbl, coords = c("long", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(3857)

# Use this helper if you want to regenerate a fresh offset table:
# generate_random_offsets(site_sf_3857, map_bbox)
offset_tbl <- tribble(
  ~point_id, ~population_code, ~dx, ~dy,
  "IN_01", "Pgl346_IN", -50000, -10000,
  "KY_02", "Pgl347_KY", -40000, 10000,
  "EX_01", "Pgl348_EXS", 55789.7, 35000,
  "KY_01", "Pgl349_KY", -115000, 50000,
  "KY_03", "Pgl350_KY", -80000, -45000,
  "KY_04", "Pgl351_KY", -25000, -60000,
  "KY_05", "Pgl352_KY", -40000, -90000.0,
  "EX_02", "Pgl353_EXS", 30000, -70000,
  "KY_06", "Pgl354_KY", 50000, -40000,
  "EX_03", "Pgl355_EXS", 50002.7, 20208.8,
  "TNE_02", "Pgl356_TNE", 130000, -15000,
  "TNE_01", "Pgl357_TNE", -8000, -61579.9,
  "TNW_04", "Pgl358_TNW", 107868.5, 50000,
  "TNW_05", "Pgl359_TNW", -50000, -40000,
  "TNW_03", "Pgl360_TNW", 81983.0, 90000,
  "TNW_02", "Pgl361_TNW", 50000, 60000,
  "TNW_01", "Pgl362_TNW", -8547.2, 60000
)

state_sf <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(postal %in% c("IN", "KY", "TN", "IL")) %>%
  mutate(
    state_fill = recode(
      postal,
      "IL" = "#f2f2f2",
      "IN" = "#e6e6e6",
      "KY" = "#d9d9d9",
      "TN" = "#cccccc"
    )
  ) %>%
  st_transform(3857) %>%
  st_intersection(map_bbox_sfc)

state_labels <- tibble(
  label = c("Illinois", "Indiana", "Kentucky", "Tennessee"),
  long = c(-88.5, -86.2, -86, -86.75),
  lat = c(38.35, 38.4, 37.25, 35.75),
  angle = c(70, 0, 0, 0)
) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3857) %>%
  cbind(st_coordinates(.)) %>%
  as_tibble()

largewater_sf <- st_read(largewater_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(3857) %>%
  st_intersection(map_bbox_sfc)

smallwater_sf <- st_read(smallwater_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(3857) %>%
  st_intersection(map_bbox_sfc)

point_xy_base <- site_sf_3857 %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  left_join(offset_tbl, by = c("point_id", "population_code")) %>%
  mutate(
    X_offset = X + dx,
    Y_offset = Y + dy
  )

if (any(is.na(point_xy_base$dx) | is.na(point_xy_base$dy))) {
  stop("offset_tbl is missing dx/dy values for one or more points.")
}

line_tbl <- point_xy_base %>%
  select(point_id, population_code, X_anchor = X, Y_anchor = Y, X_offset, Y_offset) %>%
  pivot_longer(
    cols = c(X_anchor, Y_anchor, X_offset, Y_offset),
    names_to = c(".value", "position"),
    names_pattern = "([XY])_(anchor|offset)"
  ) %>%
  mutate(position = factor(position, levels = c("anchor", "offset"))) %>%
  arrange(point_id, position)

base_pie_radius <- min(
  as.numeric(map_bbox["xmax"] - map_bbox["xmin"]),
  as.numeric(map_bbox["ymax"] - map_bbox["ymin"])
) * 0.028

for (k in k_values) {
  structure_file <- file.path(structure_dir, sprintf("STRUCTURE_K%02d_individual_mean.csv", k))
  structure_tbl <- read_structure_k(structure_file, k)
  cluster_cols <- paste0("cluster_", seq_len(k))

  point_xy <- point_xy_base %>%
    left_join(structure_tbl, by = c("point_id" = "population_label")) %>%
    mutate(
      pie_radius = scales::rescale(
        n_samples,
        to = c(base_pie_radius * 0.9, base_pie_radius * 1.7)
      )
    )

  if (any(!stats::complete.cases(point_xy[, cluster_cols, drop = FALSE]))) {
    stop(sprintf("K = %d STRUCTURE values are missing for one or more Wen labels.", k))
  }

  p <- ggplot() +
    geom_sf(
      data = state_sf,
      aes(fill = state_fill),
      color = NA
    ) +
    geom_sf(
      data = largewater_sf,
      fill = "#1f78b4",
      color = "#1f78b4",
      linewidth = 0.25,
      alpha = 1
    ) +
    geom_sf(
      data = smallwater_sf,
      color = "#1f78b4",
      linewidth = 0.45,
      alpha = 1,
      lineend = "round"
    ) +
    geom_text(
      data = state_labels,
      aes(x = X, y = Y, label = label, angle = angle),
      inherit.aes = FALSE,
      size = 5.2,
      color = "#4d4d4d",
      fontface = "bold",
      family = "sans"
    ) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    geom_line(
      data = line_tbl,
      aes(x = X, y = Y, group = point_id),
      color = "gray30",
      linewidth = 0.35,
      alpha = 0.9
    ) +
    geom_point(
      data = point_xy,
      aes(x = X, y = Y),
      color = "black",
      fill = "black",
      shape = 16,
      size = 1.7,
      stroke = 0
    ) +
    scatterpie::geom_scatterpie(
      data = point_xy,
      aes(x = X_offset, y = Y_offset, r = pie_radius),
      cols = cluster_cols,
      color = "black",
      linewidth = 0.2,
      alpha = 1,
      show.legend = FALSE
    ) +
    geom_text(
      data = point_xy,
      aes(x = X_offset, y = Y_offset - pie_radius * 1.3, label = point_id),
      fontface = "bold",
      show.legend = FALSE,
      size = 3.4
    ) +
    ggokabeito::scale_fill_okabe_ito(
      breaks = cluster_cols,
      labels = paste("Cluster", seq_len(k)),
      name = sprintf("STRUCTURE K = %d", k)
    ) +
    annotation_scale(location = "br", width_hint = 0.22) +
    coord_sf(
      xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
      ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
      expand = FALSE,
      crs = 3857
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#f2f2f2", color = NA),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.text = element_text(size = 9, color = "gray25"),
      axis.ticks = element_line(color = "gray35", linewidth = 0.3),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )

  ggsave(
    file.path(repo_root, "figure", sprintf("P_globosa_water_overlay_scatterpie_K%02d.png", k)),
    plot = p,
    width = 10,
    height = 8,
    dpi = 900
  )
}
