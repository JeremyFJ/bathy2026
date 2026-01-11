library(ape)
library(ggtree)
library(tidyverse)
library(phytools)

acc_map <- read.csv("../data/COX1_accession_species.csv", stringsAsFactors = FALSE)
tree    <- read.tree("../data/COX1_aln_trimmed.fasta.contree")
tips    <- tree$tip.label

clean_tip <- function(x) {
  x <- gsub("^_R_", "", x)          # drop any leading '_R_'
  x <- sub("\\.\\d+$", "", x)       # drop version suffix .1, .2, etc
  x
}

tips_clean         <- clean_tip(tips)
acc_map$acc_clean  <- clean_tip(acc_map$accession_number)

# Match each tip to a species_name
species_matched <- acc_map$species_name[match(tips_clean, acc_map$acc_clean)]

# If only your own sequence is NA, assign that; otherwise you might want to be stricter.
# Here I'll assume there is one NA corresponding to your Bsp_COI tip:
na_idx <- which(is.na(species_matched))
species_matched[na_idx] <- "Bathynomus cf. giganteus"

# Create a pure species vector (what each tip "is")
species_vec <- species_matched
# For any remaining NAs (things not in the DB), fall back to original label
species_vec[is.na(species_vec)] <- tips[is.na(species_vec)]

# Midpoint root the original tree
tree_mid <- phytools::midpoint.root(tree)

# Keep the first occurrence of each species
keep <- !duplicated(species_vec)

# Drop all other accessions
tree_species <- drop.tip(tree_mid, tree_mid$tip.label[!keep])

# Give each remaining tip the species name only
tree_species$tip.label <- species_vec[keep]

# Make labels unique just in case (e.g., "Bathynomus giganteus", "Bathynomus giganteus.1")
tree_species$tip.label <- make.unique(tree_species$tip.label)

abbrev_species <- function(x) {
  x <- gsub("^Bathynomus\\s+", "B. ", x)
  x
}

format_tip_label <- function(x) {
  # Abbreviate genus
  x <- gsub("^Bathynomus\\s+", "B. ", x)
  
  # Escape single quotes for safety
  x <- gsub("'", "\\\\'", x)
  
  # Case 1: B. sp. or B. spp.
  if (grepl("^B\\.\\s+sp\\.$", x)) {
    return("paste(italic('B.'), plain(' sp.'))")
  }
  if (grepl("^B\\.\\s+spp\\.$", x)) {
    return("paste(italic('B.'), plain(' spp.'))")
  }
  
  # Case 2: everything else → fully italic
  paste0("italic('", x, "')")
}


tree_species$tip.label <- vapply(tree_species$tip.label, format_tip_label, character(1))
highlight_label <- format_tip_label("Bathynomus cf. giganteus")

library(ggplot2)

p_species_cox <- ggtree(tree_species) +
  geom_tree() +
  ggtitle("COX1") +
  geom_tiplab(
    aes(color = label == highlight_label),
    size = 3.5,
    parse = TRUE
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    guide = "none"
  ) +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    plot.margin = margin(5.5, 80, 5.5, 5.5, "pt")  # big right margin
  ) +
  coord_cartesian(clip = "off")                    # <-- key line

p_species_cox

ggsave("../figures/COX1_phylogeny.pdf", p_species, width = 13, height = 10)