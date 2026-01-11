# ============================================================
# Demo: Global distribution of rare Bathynomus species
# (< 50 observations)
# ============================================================
# Purpose:
#   - Visualize sparse global distributions for rare species
#   - Integrate morphological and molecular observations
#   - Separate Atlantic vs Indo-Pacific hemispheres
#
# Inputs:
#   - ../data/obs_db.csv
#   - ../data/barcodes_db.csv
#
# Output:
#   - ../figures/global_distribution_sparse_bathynomus.png
# ============================================================

# ------------------------
# 0. Libraries
# ------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(ggtext)
library(Polychrome)

# ------------------------
# 1. Load data
# ------------------------
obs <- read.csv("../data/obs_db.csv", stringsAsFactors = FALSE)
barcodes <- read.csv("../data/barcodes_db.csv", stringsAsFactors = FALSE)

# ------------------------
# 2. QC check
# ------------------------
stopifnot(
  all(c("species","latitude","longitude","aquarium_inp","id_correct") %in% colnames(obs)),
  all(c("species_name","latitude","longitude") %in% colnames(barcodes))
)

# ------------------------
# 3. Conservative filtering (same logic as main map)
# ------------------------
obs_clean <- obs %>%
  filter(
    !is.na(latitude), !is.na(longitude),
    latitude >= -90, latitude <= 90,
    longitude >= -180, longitude <= 180,
    (is.na(aquarium_inp) | aquarium_inp == FALSE),
    (id_correct == '' | id_correct == "yes"),
    species != "Bathynomus paracelensis" # likely junior synonym of B. vaderi (in review)
  )

# ------------------------
# 4. Morphological observations
# ------------------------
morph <- obs_clean %>%
  transmute(
    species   = as.character(species),
    latitude  = as.numeric(latitude),
    longitude = as.numeric(longitude),
    obs_type  = "morphological"
  ) %>%
  filter(!is.na(species), species != "")

# ------------------------
# 5. Molecular (barcode) observations
# ------------------------
mol <- barcodes %>%
  transmute(
    species   = as.character(species_name),
    latitude  = as.numeric(latitude),
    longitude = as.numeric(longitude),
    obs_type  = "molecular"
  ) %>%
  filter(!is.na(species), species != "")

# ------------------------
# 6. Combine and annotate
# ------------------------
plot_df <- bind_rows(morph, mol) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(
    hemisphere = ifelse(longitude < 0, "Western", "Eastern"),
    obs_type   = factor(obs_type, levels = c("morphological", "molecular"))
  )

# ------------------------
# 7. Count observations per species & hemisphere
# ------------------------
species_counts <- morph %>%
  mutate(hemisphere = ifelse(longitude < 0, "Western", "Eastern")) %>%
  count(species, hemisphere, name = "n_obs")

plot_df <- plot_df %>%
  left_join(species_counts, by = c("species", "hemisphere")) %>%
  mutate(n_obs = replace_na(n_obs, 0)) %>%
  filter(species != "Bathynomus sp.") %>%
  filter(n_obs <= 50)

# ------------------------
# 8. Color palettes
# ------------------------
species_levels <- sort(unique(plot_df$species))

pal <- Polychrome::createPalette(
  length(species_levels),
  seedcolors = c("#1b9e77", "#d95f02", "#7570b3")
)

species_cols <- setNames(pal, species_levels)

# ------------------------
# 9. Basemap & helpers
# ------------------------
world <- ne_countries(scale = "large", returnclass = "sf")
pos_jitter <- position_jitter(width = 0.75, height = 0.75)

format_species_legend <- function(x) {
  ifelse(
    grepl("^Bathynomus\\s+", x),
    paste0("*B. ", sub("^Bathynomus\\s+", "", x), "*"),
    x
  )
}

# ------------------------
# 10. Hemisphere map function
# ------------------------
make_hemi_map <- function(hemi = c("Western","Eastern")) {
  
  hemi <- match.arg(hemi)
  
  if (hemi == "Western") {
    xlim <- c(-120, 0)
    ylim <- c(-50, 50)
    title <- "Atlantic"
  } else {
    xlim <- c(25, 170)
    ylim <- c(-50, 50)
    title <- "Indo-Pacific"
  }
  
  df_hemi <- plot_df %>% filter(hemisphere == hemi)
  
  ggplot() +
    ggtitle(title) +
    geom_sf(data = world, fill = "grey95", color = "grey70") +
    
    # Morphological (filled)
    geom_point(
      data = df_hemi %>% filter(obs_type == "morphological"),
      aes(x = longitude, y = latitude, fill = species),
      shape = 21, size = 2.5,
      stroke = 0.1, color = "grey20",
      position = pos_jitter
    ) +
    
    # Molecular (outlined)
    geom_point(
      data = df_hemi %>% filter(obs_type == "molecular"),
      aes(x = longitude, y = latitude, fill = species),
      shape = 21, size = 2.7,
      stroke = 0.9, color = "black",
      position = pos_jitter
    ) +
    
    scale_fill_manual(
      values = species_cols,
      labels = format_species_legend,
      drop = FALSE
    ) +
    
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = NULL, y = NULL, fill = "Species") +
    
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = ggtext::element_markdown(size = 10),
      legend.key.size = unit(0.35, "cm"),
      plot.margin = margin(4,4,4,4)
    ) +
    guides(
      fill = guide_legend(
        nrow = 3, byrow = TRUE,
        override.aes = list(
          shape = 21, color = NA, stroke = 0, size = 3
        )
      )
    )
}

# ------------------------
# 11. Build panels
# ------------------------
wrap_panel <- function(p, label) {
  cowplot::ggdraw() +
    cowplot::draw_plot(p) +
    cowplot::draw_label(
      label, x = 0.01, y = 0.99,
      hjust = 0, vjust = 1,
      fontface = "bold", size = 14
    )
}

p_west <- wrap_panel(make_hemi_map("Western"), "F")
p_east <- wrap_panel(make_hemi_map("Eastern"), "G")

final_plot <- cowplot::plot_grid(
  p_west, p_east, ncol = 2, align = "hv"
)
final_plot
# ------------------------
# 12. Save output
# ------------------------
dir.create("../figures", showWarnings = FALSE)

ggsave(
  "../figures/global_distribution_sparse_bathynomus.jpg",
  final_plot,
  width = 15,
  height = 8,
  dpi = 300
)

message("Sparse-species global distribution figure saved.")
