library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(tigris)
# ensure robust geometry operations and normalize common shapefile attribute names
sf::sf_use_s2(FALSE)

# wrap st_read to normalize attribute names so later filtering (State, Miles) works
.original_st_read <- sf::st_read
st_read <- function(dsn, ..., quiet = TRUE) {
  sf_obj <- .original_st_read(dsn, quiet = quiet, ...)
  if (nrow(sf_obj) == 0) return(sf_obj)

  nm <- names(sf_obj)
  up <- toupper(nm)

  # normalize state column -> "State"
  state_idx <- which(up %in% c("STATE", "ST", "ST_ABBR", "STATE_ABBR", "STUSPS"))
  if (length(state_idx) && !"State" %in% nm) {
    sf_obj$State <- sf_obj[[state_idx[1]]]
  }

  # normalize miles/length column -> "Miles"
  miles_idx <- which(up %in% c("MILES", "MILE", "LENGTH", "LEN_MI", "MI"))
  if (length(miles_idx) && !"Miles" %in% nm) {
    # coerce to numeric (remove non-numeric chars)
    sf_obj$Miles <- as.numeric(gsub("[^0-9.\\-]", "", as.character(sf_obj[[miles_idx[1]]])))
  }

  sf_obj
}
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

input_file <- file.path(repo_root, "data", "Table_1.csv")
output_png <- file.path(repo_root, "figure", "Table_1_map_globosa_boundaries_rivers.png")
river_file <- file.path(
  "/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/GIS_material/Rivers/USA_Rivers_and_Streams-shp/9ae73184-d43c-4ab8-940a-c8687f61952f2020328-1-r9gw71.0odx9.shp"
)
wen_lookup <- read_wen_label_lookup(repo_root)

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
    filter(postal %in% c("IN", "KY", "TN", "IL")) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc)

  options(tigris_class = "sf")
  county_sf <- bind_rows(
    counties(state = "IN", cb = TRUE, year = 2023),
    counties(state = "KY", cb = TRUE, year = 2023),
    counties(state = "TN", cb = TRUE, year = 2023),
    counties(state = "IL", cb = TRUE, year = 2023)
  ) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc)

  rivers_all <- st_read(river_file, quiet = TRUE) %>%
    filter(State %in% c("IN", "KY", "TN", "IL")) %>%
    st_transform(3857) %>%
    st_intersection(map_bbox_sfc)

  # ~20 km buffer around sampled sites
  site_buf <- st_union(st_buffer(st_union(site_sf_3857), 20000))

  # rivers within ~20 km of any sampled site
  near_rivers <- rivers_all[st_intersects(rivers_all, site_buf, sparse = FALSE), ]

  # define "major" rivers by length (>= 50 km) and keep those connected to the near_rivers
  rivers_all$len_m <- as.numeric(st_length(rivers_all))
  major_candidates <- filter(rivers_all, len_m >= 50000)
  connected_major <- major_candidates[st_intersects(major_candidates, st_union(near_rivers), sparse = FALSE), ]
  # group major candidates into connected components (allow small gaps), then keep components that touch near_rivers
  if (nrow(major_candidates) > 0 && nrow(near_rivers) > 0) {
    # adjacency: direct intersection or small gap (1 km) to account for fragmented polylines
    adj_intersect <- st_intersects(major_candidates, major_candidates, sparse = FALSE)
    adj_close <- st_is_within_distance(major_candidates, major_candidates, dist = 1000, sparse = FALSE)
    adj <- adj_intersect | adj_close
    diag(adj) <- FALSE

    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    comps <- igraph::components(g)$membership

    # which major segments touch (or are within 1 km of) the near_rivers
    touch_near <- apply(st_is_within_distance(major_candidates, near_rivers, dist = 1000, sparse = FALSE), 1, any) |
                  apply(st_intersects(major_candidates, near_rivers, sparse = FALSE), 1, any)

    if (any(touch_near)) {
      comp_ids <- unique(comps[touch_near])
      connected_major <- major_candidates[comps %in% comp_ids, ]
    } else {
      # fallback: include any major segment within 5 km of near_rivers
      close_idx <- apply(st_is_within_distance(major_candidates, near_rivers, dist = 5000, sparse = FALSE), 1, any)
      connected_major <- major_candidates[close_idx, ]
    }
  } else {
    connected_major <- major_candidates[integer(0), ]
  }

  rivers_sf <- bind_rows(near_rivers, connected_major) %>% distinct()
  if (nrow(rivers_sf) == 0) rivers_sf <- rivers_all

  # fill by state
  fill_values <- c(
    "IN" = "#1b9e77",
    "KY" = "#d95f02",
    "TN" = "#7570b3"
  )
  # use state_group as the fill variable (keeps existing aes(fill = ecoregion) unchanged)
  site_sf_3857$ecoregion <- factor(site_sf$state_group, levels = c("IN", "KY", "TN"))

  ggplot() +
    geom_sf(data = county_sf, fill = "gray99", color = "gray82", linewidth = 0.18) +
    geom_sf(data = state_sf, fill = NA, color = "gray35", linewidth = 0.25) +
    geom_sf(
      data = site_sf_3857,
      aes(fill = ecoregion, size = n_samples),
      shape = 21,
      color = "black",
      stroke = 0.3,
      alpha = 0.95
    ) +
    geom_sf(data = rivers_sf, color = "#4aa3df", linewidth = 0.55, alpha = 0.9) +
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
    scale_fill_manual(values = fill_values, drop = FALSE, name = "Ecoregion") +
    scale_size_continuous(range = c(2.2, 5.5), name = "Sample size") +
    annotation_scale(location = "br", width_hint = 0.22) +
    coord_sf(
      xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
      ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
      expand = FALSE,
      crs = 3857
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
    )
}

tbl <- read_table_1(input_file)
p <- build_map(tbl)

ggsave(output_png, plot = p, width = 10, height = 8, dpi = 900)
