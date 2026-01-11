# ============================================================
# Demo: Global distribution of Bathynomus observations
# ============================================================
# Purpose:
#   - Visualize global distributions of Bathynomus species
#   - Integrate morphological and molecular observations
#   - Demonstrate conservative filtering of questionable records
#
# Inputs:
#   - ../data/obs_db.csv
#   - ../data/barcodes_db.csv
#
# Output:
#   - ../figures/global_distribution_bathynomus.png
# ============================================================

# ------------------------
# 0. Libraries
# ------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(ggtext)
library(cowplot)
library(scales)

# ------------------------
# 1. Load data
# ------------------------
obs <- read.csv("../data/obs_db.csv", stringsAsFactors = FALSE)
barcodes <- read.csv("../data/barcodes_db.csv", stringsAsFactors = FALSE)

# ------------------------
# 2. Basic sanity checks
# ------------------------
required_obs_cols <- c(
  "latitude", "longitude", "species",
  "date", "aquarium_inp", "id_correct"
)

missing_cols <- setdiff(required_obs_cols, colnames(obs))
if (length(missing_cols) > 0) {
  stop("obs_db.csv is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

required_bc_cols <- c(
  "latitude", "longitude", "species_name", "date"
)

missing_bc <- setdiff(required_bc_cols, colnames(barcodes))
if (length(missing_bc) > 0) {
  stop("barcodes_db.csv is missing required columns: ",
       paste(missing_bc, collapse = ", "))
}

# ------------------------
# 3. Remove inland records
# ------------------------
obs_clean <- obs %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  sharkPulseR::removeInlands3(buffer = -10) %>%
  filter(inland == FALSE)

# ------------------------
# 4. Filter observations by confidence
# ------------------------
obs_id <- obs_clean %>%
  filter(
    (is.na(aquarium_inp) | aquarium_inp == FALSE),
    (id_correct == '' | id_correct == "yes"),
    species != "Bathynomus paracelensis"   # likely a junior synonym of B. vaderi (in review)
  )

# ------------------------
# 5. Morphological observations
# ------------------------
obs_morph <- obs_id %>%
  mutate(
    parsed_date = parse_date_time(
      date,
      orders = c(
        "m/d/y",   # <- critical for your data
        "m/d/Y",
        "Y-m-d", "Y/m/d",
        "d-m-Y", "d/m/Y",
        "b-Y", "Y"
      ),
      truncated = 1
    )
  ) %>%
  transmute(
    species,
    latitude,
    longitude,
    date = as.Date(parsed_date),
    depth = suppressWarnings(as.numeric(depth)),
    obs_type = "morphological"
  )

# ------------------------
# 6. Molecular (barcode) observations
# ------------------------
obs_mol <- barcodes %>%
  mutate(
    parsed_date = parse_date_time(
      date,
      orders = c(
        "d-b-y",   # 24-Jan-25
        "d-b-Y",   # 24-Jan-2025 (just in case)
        "b-y",     # Feb-19
        "b-Y",     # Feb-2019
        "Y-m-d"    # fallback ISO
      ),
      truncated = 1
    )
  ) %>%
  transmute(
    species   = species_name,
    latitude,
    longitude,
    date      = as.Date(parsed_date),
    depth     = NA_real_,
    obs_type  = "molecular"
  )

# ------------------------
# 7. Combine datasets
# ------------------------
plot_df <- bind_rows(obs_morph, obs_mol) %>%
  filter(
    !is.na(species),
    !is.na(latitude),
    !is.na(longitude)
  )

# ------------------------
# 8. World basemap
# ------------------------
world <- ne_countries(scale = "large", returnclass = "sf")

pos_jitter <- position_jitter(width = 0.5, height = 0.5)

# ------------------------
# 9. Helper functions
# ------------------------
abbrev_species <- function(sp) {
  if (grepl("^Bathynomus\\s+sp", sp, ignore.case = TRUE)) {
    "*Bathynomus* spp."
  } else if (grepl("^Bathynomus\\s+", sp)) {
    paste0("*B. ", sub("^Bathynomus\\s+", "", sp), "*")
  } else {
    sp
  }
}

depth_plot <- function(df) {
  if (all(is.na(df$depth))) return(NULL)
  ggplot(df %>% filter(!is.na(depth)), aes(x = 1, y = depth)) +
    geom_boxplot(width = 0.3, outlier.size = 1.2) +
    scale_x_continuous(limits = c(0.5, 1.5)) +
    labs(y = "Depth (m)", x = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

time_plot <- function(df) {
  if (all(is.na(df$date))) return(NULL)
  ggplot(df, aes(x = date)) +
    geom_histogram(bins = 20, fill = "grey60", color = "grey30") +
    labs(x = "Year", y = "n") +
    theme_minimal(base_size = 10)
}

plot_species <- function(sp, df) {
  
  df_sp <- df %>%
    filter(species == sp)
  
  if (nrow(df_sp) == 0) return(NULL)
  
  xlim <- range(df_sp$longitude)
  ylim <- range(df_sp$latitude)
  
  p_map <- ggplot() +
    geom_sf(data = world, fill = "grey93", color = "grey60") +
    geom_point(
      data = df_sp,
      aes(
        x = longitude,
        y = latitude,
        shape = obs_type,
        color = obs_type
      ),
      size = 2.3,
      position = pos_jitter
    ) +
    scale_shape_manual(values = c(morphological = 21, molecular = 24)) +
    scale_color_manual(values = c(morphological = "grey40",
                                  molecular = "black")) +
    coord_sf(
      xlim = xlim + c(-12, 12),
      ylim = ylim + c(-12, 12),
      expand = FALSE
    ) +
    labs(title = abbrev_species(sp)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = ggtext::element_markdown()
    )
  
  p_map
}

# ------------------------
# 10. Generate plots
# ------------------------
plot_df_gt50 <- plot_df %>%
  add_count(species) %>%
  filter(n > 50)

species_list <- sort(unique(plot_df_gt50$species))

plots <- lapply(species_list, plot_species, df = plot_df)
names(plots) <- species_list

# ------------------------
# 11. Assemble panel figure
# ------------------------
bathy_spp <- species_list[
  grepl("^Bathynomus\\s+sp", species_list, ignore.case = TRUE)
]

top_species <- setdiff(species_list, bathy_spp)[1:4]

wrap_panel <- function(p, label) {
  cowplot::ggdraw() +
    cowplot::draw_plot(p) +
    cowplot::draw_label(
      label, x = 0.01, y = 0.99,
      hjust = 0, vjust = 1,
      fontface = "bold", size = 14
    )
}

pA <- wrap_panel(plots[[bathy_spp]], "A")
pB <- wrap_panel(plots[[top_species[1]]], "B")
pC <- wrap_panel(plots[[top_species[2]]], "C")
pD <- wrap_panel(plots[[top_species[3]]], "D")
pE <- wrap_panel(plots[[top_species[4]]], "E")

final_plot <- plot_grid(
  pA,
  plot_grid(pB, pC, pD, pE, ncol = 2),
  ncol = 1,
  rel_heights = c(1, 2)
)
final_plot
# ------------------------
# 12. Save output
# ------------------------
# dir.create("../figures", showWarnings = FALSE)

ggsave(
  "../figures/global_distribution_bathynomus.jpg",
  final_plot,
  width = 11.5,
  height = 14,
  dpi = 300
)

message("Global distribution figure saved.")
