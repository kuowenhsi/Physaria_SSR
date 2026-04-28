library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)

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
output_png <- file.path(repo_root, "figure", "Table_1_map_natural_earth.png")
rivers_base_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_rivers_lake_centerlines",
  "ne_10m_rivers_lake_centerlines.shp"
)
wen_lookup <- read_wen_label_lookup(repo_root)
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
region_polys_file <- file.path(
  repo_root,
  "data",
  "natural_earth_10m",
  "ne_10m_geography_regions_polys",
  "ne_10m_geography_regions_polys.shp"
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
      population_code = `Population code`,
      collector = Collector,
      species_name = `Species name`,
      locality = `State, County, Locality (Element Occurrence number)`,
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
      state_group = case_when(
        str_detect(population_code, "_IN$") ~ "IN",
        str_detect(population_code, "_KY$|_EXS$") ~ "KY",
        str_detect(population_code, "_TNE$|_TNW$") ~ "TN",
        TRUE ~ NA_character_
      ),
      plot_population_label = apply_wen_label_from_canonical(population_code, wen_lookup)
    ) %>%
    filter(!is.na(lat), !is.na(long))
}

build_map <- function(tbl) {
  tbl <- tbl %>%
    filter(species_name == "P. globosa")

  site_sf <- st_as_sf(tbl, coords = c("long", "lat"), crs = 4326, remove = FALSE)
  site_sf_3857 <- st_transform(site_sf, 3857)

  label_df <- site_sf %>%
    st_drop_geometry() %>%
    group_by(population_code, plot_population_label, ecoregion) %>%
    summarise(
      long = mean(long),
      lat = mean(lat),
      n_samples = sum(n_samples),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(3857) %>%
    mutate(label = paste0(plot_population_label, " (n=", n_samples, ")"))

  label_xy <- cbind(st_drop_geometry(label_df), st_coordinates(label_df))

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
    filter(
      featurecla %in% c("River", "Lake Centerline"),
      !is.na(name),
      scalerank <= 12
    )

  lakes_sf <- combine_sf_common(
    st_read(lakes_base_file, quiet = TRUE),
    st_read(lakes_na_file, quiet = TRUE)
  ) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc) %>%
    filter(featurecla %in% c("Lake", "Reservoir"), !is.na(name))

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

  river_label_xy <- rivers_sf %>%
    filter(featurecla == "River", name == "Ohio") %>%
    st_centroid() %>%
    mutate(label = "Ohio River") %>%
    cbind(st_coordinates(.))

  region_sf <- st_read(region_polys_file, quiet = TRUE) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc) %>%
    filter(FEATURECLA %in% c("Range/mtn", "Plain", "Plateau")) %>%
    filter(NAME %in% c("APPALACHIAN MTS.", "CENTRAL LOWLAND", "CUMBERLAND PLAT."))

  region_fill <- c(
    "APPALACHIAN MTS." = "#d8c6a5",
    "CENTRAL LOWLAND" = "#e6efd7",
    "CUMBERLAND PLAT." = "#ddd2ec"
  )

  region_labels <- region_sf %>%
    st_point_on_surface() %>%
    mutate(label = str_to_title(NAME)) %>%
    cbind(st_coordinates(.))

  fill_values <- c(
    "Bluegrass" = "#d95f02",
    "Bluegrass/ partially cultivated" = "#fdb863",
    "Ex situ" = "#bdbdbd",
    "Highland Rim" = "#1b9e77",
    "Nashville Basin" = "#7570b3",
    "Shawnee Hills" = "#e7298a"
  )

  ggplot() +
    geom_sf(data = region_sf, aes(fill = NAME), color = NA, alpha = 0.18) +
    geom_sf(data = state_sf, fill = "gray98", color = "gray55", linewidth = 0.45) +
    geom_sf(
      data = old_rivers_sf,
      color = "#9fd6ef",
      linewidth = 0.3,
      alpha = 0.95,
      lineend = "round"
    ) +
    geom_sf(
      data = rivers_sf,
      color = "#2c7fb8",
      linewidth = 0.55,
      alpha = 0.9,
      lineend = "round"
    ) +
    geom_sf(
      data = lakes_sf,
      fill = "#a9d8f3",
      color = "#5f9fc9",
      linewidth = 0.2,
      alpha = 0.9
    ) +
    geom_text(
      data = river_label_xy,
      aes(x = X, y = Y, label = label),
      color = "#2c7fb8",
      size = 3.1,
      family = "sans",
      fontface = "italic"
    ) +
    geom_text(
      data = region_labels,
      aes(x = X, y = Y, label = label),
      color = "gray35",
      size = 2.8,
      family = "sans",
      fontface = "bold",
      alpha = 0.8
    ) +
    geom_sf(
      data = site_sf_3857,
      aes(fill = ecoregion, size = n_samples),
      shape = 21,
      color = "black",
      stroke = 0.3,
      alpha = 0.95
    ) +
    geom_text_repel(
      data = label_xy,
      aes(x = X, y = Y, label = label),
      size = 3,
      min.segment.length = 0.1,
      seed = 42,
      max.overlaps = 50,
      box.padding = 0.25,
      point.padding = 0.2,
      segment.color = "gray40"
    ) +
    scale_fill_manual(
      values = c(region_fill, fill_values),
      breaks = names(fill_values),
      drop = FALSE,
      name = "Ecoregion"
    ) +
    scale_size_continuous(range = c(2.2, 5.5), name = "Sample size") +
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
      title = "Physaria globosa sampling map",
      subtitle = "Lakes, 10m rivers, and local fine river paths near sampling sites"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

tbl <- read_table_1(input_file)
p <- build_map(tbl)

ggsave(output_png, plot = p, width = 10, height = 8, dpi = 900)
