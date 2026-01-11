# ============================================================
# Demo: Sourcing and standardizing GBIF Bathynomus records
# ============================================================
# Purpose:
#   - Demonstrate retrieval of GBIF records
#   - Illustrate conservative species-name cleaning
#   - Output a standardized table suitable for database import
#
# Requirements:
#   - No authentication required
#   - Uses only public GBIF API
#
# Output:
#   - data/gbif_bathynomus_demo.csv
# ============================================================

library(rgbif)
library(dplyr)
library(purrr)
library(tibble)
library(lubridate)

#-------------------------
# 1. Resolve GBIF taxonKey
#-------------------------
bathynomus_key <- name_backbone(
  name  = "Bathynomus",
  genus = "Bathynomus"
)$usageKey

message("Bathynomus GBIF taxonKey: ", bathynomus_key)

#-------------------------
# 2. Download all GBIF records (paged)
#-------------------------
get_all_gbif_records <- function(taxon_key, chunk_size = 300) {
  
  meta <- occ_search(
    taxonKey = taxon_key,
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE,
    limit = 0
  )
  
  total_n <- meta$meta$count
  starts  <- seq(0, total_n - 1, by = chunk_size)
  
  map_dfr(starts, function(st) {
    occ_search(
      taxonKey = taxon_key,
      hasCoordinate = TRUE,
      hasGeospatialIssue = FALSE,
      limit = chunk_size,
      start = st
    )$data %>%
      as_tibble()
  })
}

gbif_raw <- get_all_gbif_records(bathynomus_key)

#-------------------------
# 3. Conservative species cleaning
#-------------------------
clean_species <- function(x) {
  vapply(x, function(nm) {
    if (is.na(nm) || !startsWith(nm, "Bathynomus")) return("")
    parts <- strsplit(nm, "\\s+")[[1]]
    if (length(parts) < 2) return("")
    if (!grepl("^[a-z]+(-[a-z]+)?$", parts[2])) return("")
    paste(parts[1], parts[2])
  }, character(1))
}

#-------------------------
# 4. Standardize for database ingestion
#-------------------------
bathynomus_gbif <- gbif_raw %>%
  mutate(
    latitude  = as.numeric(decimalLatitude),
    longitude = as.numeric(decimalLongitude),
    date      = as.Date(eventDate),
    species   = clean_species(scientificName),
    n         = coalesce(individualCount, 1L),
    location  = coalesce(locality, waterBody, country),
    reference = coalesce(
      bibliographicCitation,
      paste0("Dataset: ", datasetName)
    ),
    url       = paste0("https://www.gbif.org/occurrence/", key),
    source_id = coalesce(occurrenceID, as.character(gbifID)),
    notes     = paste(
      "orig_name:", scientificName,
      "| basis:", basisOfRecord,
      "| institution:", institutionCode,
      "| catalog:", catalogNumber
    )
  ) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  transmute(
    latitude,
    longitude,
    date,
    species,
    url,
    reference,
    n = as.integer(n),
    location,
    depth = as.numeric(depth),
    notes,
    source = "GBIF"
  )

#-------------------------
# 5. Write reproducible output
#-------------------------
# dir.create("data", showWarnings = FALSE)

write.csv(
  bathynomus_gbif,
  "../data/gbif_bathynomus_demo.csv",
  row.names = FALSE
)

message("Data sourced in: data/gbif_bathynomus_demo.csv")
