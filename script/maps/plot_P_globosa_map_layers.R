library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

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
shape_dir <- file.path(repo_root, "data", "P_globosa_map")
output_dir <- file.path(repo_root, "figure", "P_globosa_map_layers")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
wen_lookup <- read_wen_label_lookup(repo_root)

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
      n_samples = readr::parse_number(as.character(n_samples)),
      plot_population_label = apply_wen_label_from_canonical(population_code, wen_lookup)
    ) %>%
    filter(!is.na(lat), !is.na(long))
}

compute_map_extent <- function(tbl) {
  tbl <- tbl %>% filter(species_name == "P. globosa")
  site_sf <- st_as_sf(tbl, coords = c("long", "lat"), crs = 4326, remove = FALSE)
  site_sf_3857 <- st_transform(site_sf, 3857)
  bbox <- st_bbox(site_sf_3857)
  x_pad <- (bbox$xmax - bbox$xmin) * 0.14
  y_pad <- (bbox$ymax - bbox$ymin) * 0.16
  map_bbox <- st_bbox(c(
    xmin = as.numeric(bbox["xmin"]) - as.numeric(x_pad),
    ymin = as.numeric(bbox["ymin"]) - as.numeric(y_pad),
    xmax = as.numeric(bbox["xmax"]) + as.numeric(x_pad),
    ymax = as.numeric(bbox["ymax"]) + as.numeric(y_pad)
  ), crs = st_crs(site_sf_3857))

  list(
    site_sf_3857 = site_sf_3857,
    bbox = map_bbox,
    bbox_sfc = st_as_sfc(map_bbox)
  )
}

style_layer <- function(geom_type) {
  if (any(str_detect(geom_type, "POINT"))) {
    list(kind = "point", color = "#6a3d9a", fill = "#cab2d6")
  } else if (any(str_detect(geom_type, "LINE"))) {
    list(kind = "line", color = "#1f78b4")
  } else {
    list(kind = "polygon", color = "#4d4d4d", fill = "#a6cee3")
  }
}

plot_single_layer <- function(layer_sf, layer_name, state_sf, site_sf_3857, map_bbox) {
  geom_type <- unique(as.character(st_geometry_type(layer_sf)))
  layer_style <- style_layer(geom_type)

  base_plot <- ggplot() +
    geom_sf(data = state_sf, fill = "gray98", color = "gray60", linewidth = 0.4)

  if (layer_style$kind == "point") {
    base_plot <- base_plot +
      geom_sf(data = layer_sf, color = layer_style$color, fill = layer_style$fill, shape = 21, size = 2.2, stroke = 0.3)
  } else if (layer_style$kind == "line") {
    base_plot <- base_plot +
      geom_sf(data = layer_sf, color = layer_style$color, linewidth = 0.5, alpha = 0.95)
  } else {
    base_plot <- base_plot +
      geom_sf(data = layer_sf, fill = layer_style$fill, color = layer_style$color, linewidth = 0.25, alpha = 0.8)
  }

  base_plot +
    geom_sf(
      data = site_sf_3857,
      shape = 21,
      fill = "#f03b20",
      color = "black",
      stroke = 0.25,
      size = 2.1,
      alpha = 0.95
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
      y = NULL,
      title = layer_name,
      subtitle = paste("Geometry:", paste(geom_type, collapse = ", "))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

tbl <- read_table_1(input_file)
extent_info <- compute_map_extent(tbl)

state_sf <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(postal %in% c("IN", "KY", "TN")) %>%
  st_transform(3857) %>%
  st_intersection(extent_info$bbox_sfc)

shape_files <- list.files(shape_dir, pattern = "\\.shp$", full.names = TRUE) %>%
  sort()

for (shape_file in shape_files) {
  layer_name <- tools::file_path_sans_ext(basename(shape_file))
  layer_sf <- st_read(shape_file, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(3857) %>%
    st_intersection(extent_info$bbox_sfc)

  if (nrow(layer_sf) == 0) {
    next
  }

  p <- plot_single_layer(
    layer_sf = layer_sf,
    layer_name = layer_name,
    state_sf = state_sf,
    site_sf_3857 = extent_info$site_sf_3857,
    map_bbox = extent_info$bbox
  )

  ggsave(
    filename = file.path(output_dir, paste0(layer_name, ".png")),
    plot = p,
    width = 10,
    height = 8,
    dpi = 900
  )
}
