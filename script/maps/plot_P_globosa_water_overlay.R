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

input_file <- file.path(repo_root, "data", "Table_1.csv")
largewater_file <- file.path(repo_root, "data", "P_globosa_map", "P_glob_largewater.shp")
smallwater_file <- file.path(repo_root, "data", "P_globosa_map", "P_glob_smallwater.shp")
output_png <- file.path(repo_root, "figure", "P_globosa_water_overlay.png")
output_png_points <- file.path(repo_root, "figure", "P_globosa_water_overlay_points.png")

read_table_1 <- function(path) {
  raw_lines <- readr::read_lines(path)
  data_lines <- raw_lines[!str_detect(raw_lines, '^"?#') & raw_lines != ""]

  readr::read_csv(I(data_lines), show_col_types = FALSE) %>%
    rename(
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

tbl <- read_table_1(input_file)
map_bbox <- compute_map_extent(tbl)
map_bbox_sfc <- st_as_sfc(map_bbox)
globosa_tbl <- tbl %>% filter(species_name == "P. globosa")
site_sf_3857 <- st_as_sf(globosa_tbl, coords = c("long", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(3857)

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

# Hard-coded state label coordinates in lon/lat so they can be edited easily.
state_labels <- tibble(
  label = c("Illinois", "Indiana", "Kentucky", "Tennessee"),
  long = c(-88.5, -86.2, -86, -86.5),
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

fill_values <- c(
  "Bluegrass" = "#d95f02",
  "Bluegrass/ partially cultivated" = "#fdb863",
  "Ex situ" = "#bdbdbd",
  "Highland Rim" = "#1b9e77",
  "Nashville Basin" = "#7570b3",
  "Shawnee Hills" = "#e7298a"
)

base_plot <- ggplot() +
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
    title = "P. globosa water overlay",
    subtitle = "P_glob_largewater polygons and P_glob_smallwater lines"
  ) +
  scale_fill_identity() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#f2f2f2", color = NA),
    axis.text = element_text(size = 9, color = "gray25"),
    axis.ticks = element_line(color = "gray35", linewidth = 0.3),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave(output_png, plot = base_plot, width = 10, height = 8, dpi = 900)

p_points <- base_plot +
  ggnewscale::new_scale_fill() +
  geom_sf(
    data = site_sf_3857,
    aes(fill = ecoregion, size = n_samples),
    shape = 21,
    color = "black",
    stroke = 0.3,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = names(fill_values),
    drop = FALSE,
    name = "Ecoregion"
  ) +
  scale_size_continuous(range = c(2.2, 5.5), name = "Sample size") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )+
  coord_sf(
    xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
    ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
    expand = FALSE,
    crs = 3857
  )

ggsave(output_png_points, plot = p_points, width = 10, height = 8, dpi = 900)
