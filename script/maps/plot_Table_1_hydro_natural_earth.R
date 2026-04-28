library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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

input_file <- file.path(repo_root, "data", "Table_1.csv")
output_png <- file.path(repo_root, "figure", "Table_1_hydro_natural_earth.png")
rivers_base_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_rivers_lake_centerlines",
  "ne_10m_rivers_lake_centerlines.shp"
)
rivers_na_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_rivers_north_america",
  "ne_10m_rivers_north_america.shp"
)
lakes_base_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_lakes",
  "ne_10m_lakes.shp"
)
lakes_na_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_lakes_north_america",
  "ne_10m_lakes_north_america.shp"
)
old_rivers_file <- file.path(
  "/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden",
  "General - IMLS National Leadership Grant 2023",
  "GIS_material",
  "Rivers",
  "USA_Rivers_and_Streams-shp",
  "9ae73184-d43c-4ab8-940a-c8687f61952f2020328-1-r9gw71.0odx9.shp"
)

combine_sf_common <- function(...) {
  layers <- list(...)
  common_cols <- Reduce(intersect, lapply(layers, names))
  bind_rows(lapply(layers, function(x) x[, common_cols]))
}

read_table_1 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
    rename(
      species_name = `Species name`,
      lat_raw = Lat,
      long_raw = Long
    ) %>%
    mutate(
      species_name = str_replace_all(species_name, "[“”]", ""),
      species_name = str_squish(species_name),
      lat = readr::parse_number(as.character(lat_raw)),
      long = readr::parse_number(as.character(long_raw))
    ) %>%
    filter(species_name == "P. globosa", !is.na(lat), !is.na(long))
}

build_hydro_map <- function(tbl) {
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
  map_bbox_sfc <- st_as_sfc(map_bbox)

  state_sf <- ne_states(country = "United States of America", returnclass = "sf") %>%
    filter(postal %in% c("IN", "KY", "TN")) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc)

  rivers_sf <- combine_sf_common(
    st_read(rivers_base_file, quiet = TRUE),
    st_read(rivers_na_file, quiet = TRUE)
  ) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc) %>%
    filter(featurecla %in% c("River", "Lake Centerline"), scalerank <= 12)

  lakes_sf <- combine_sf_common(
    st_read(lakes_base_file, quiet = TRUE),
    st_read(lakes_na_file, quiet = TRUE)
  ) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc) %>%
    filter(featurecla %in% c("Lake", "Reservoir"))

  site_buffer <- site_sf_3857 %>%
    st_union() %>%
    st_buffer(10000)

  old_rivers_sf <- st_read(old_rivers_file, quiet = TRUE) %>%
    filter(
      State %in% c("IN", "KY", "TN", "IL"),
      Feature %in% c("Stream", "Artificial Path", "Canal/Ditch", "Connector")
    ) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc)

  old_rivers_sf <- old_rivers_sf[
    st_intersects(old_rivers_sf, site_buffer, sparse = FALSE)[, 1],
  ]

  ggplot() +
    geom_sf(data = state_sf, fill = "gray97", color = "gray70", linewidth = 0.45) +
    geom_sf(data = old_rivers_sf, color = "#9fd6ef", linewidth = 0.3, alpha = 0.95, lineend = "round") +
    geom_sf(data = rivers_sf, color = "#2c7fb8", linewidth = 0.5, alpha = 0.95, lineend = "round") +
    geom_sf(data = lakes_sf, aes(fill = featurecla), color = "#66a9d1", linewidth = 0.25, alpha = 0.95) +
    scale_fill_manual(values = c("Lake" = "#9fd3f2", "Reservoir" = "#7fb8dd"), name = NULL) +
    coord_sf(
      xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
      ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
      expand = FALSE,
      crs = 3857
    ) +
    labs(
      title = "Hydrography in the Table 1 map extent",
      subtitle = "Lakes, 10m rivers, and local fine river paths within 10 km of sampling sites",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

table_1 <- read_table_1(input_file)
hydro_map <- build_hydro_map(table_1)

ggsave(output_png, plot = hydro_map, width = 10, height = 8, dpi = 900)
