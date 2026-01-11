# ============================================================
# Demo: Sourcing mitochondrial barcode data for Bathynomus
# ============================================================
# Purpose:
#   - Demonstrate reproducible retrieval of mitochondrial
#     barcode sequences from NCBI GenBank
#   - Extract basic metadata (locus, coordinates, date)
#   - Write standardized output to CSV
#
# Requirements:
#   - NCBI Entrez API key (free)
#     https://www.ncbi.nlm.nih.gov/account/settings/
#
# Output:
#   - data/ncbi_bathynomus_mito_demo.csv
# ============================================================

library(rentrez)
library(stringr)
library(dplyr)
library(tibble)

#-------------------------
# 0. Set NCBI API key
#-------------------------
# IMPORTANT:
# Replace "YOUR_API_KEY_HERE" with your own NCBI API key.
# This increases rate limits and ensures stable access.
#
# Obtain a key at:
# https://www.ncbi.nlm.nih.gov/account/settings/

Sys.setenv(ENTREZ_KEY = "YOUR_API_KEY_HERE")

#-------------------------
# 1. Helper: parse lat / lon
#-------------------------
parse_lat_lon <- function(coord) {
  if (is.na(coord) || coord == "") return(NA_real_)
  value <- as.numeric(gsub("[^0-9.-]", "", coord))
  dir   <- gsub("[^NSEW]", "", coord)
  if (dir %in% c("S", "W")) value <- -value
  value
}

#-------------------------
# 2. Extract locus from GenBank record
#-------------------------
extract_locus <- function(record_text) {
  
  # Extract the INSDSeq_definition field
  def <- sub(
    ".*<INSDSeq_definition>(.*?)</INSDSeq_definition>.*",
    "\\1",
    record_text
  )
  
  if (identical(def, record_text) || is.na(def)) {
    return("unknown")
  }
  
  def <- tolower(def)
  
  # Canonical mitochondrial loci
  if (grepl("cytochrome c oxidase subunit i|cox1|coi", def)) return("cox1")
  if (grepl("cytochrome b|cytb", def)) return("cytb")
  if (grepl("16s", def)) return("16s_rrna")
  if (grepl("12s", def)) return("12s_rrna")
  if (grepl("nd1", def)) return("nd1")
  if (grepl("nd2", def)) return("nd2")
  if (grepl("nd4", def)) return("nd4")
  if (grepl("nd5", def)) return("nd5")
  if (grepl("subunit ribosomal", def)) return("ribosomal subunit")
  
  # Whole mitogenomes
  if (grepl("mitochondrial genome|complete genome", def))
    return("mitogenome")
  
  return("unknown")
}



#-------------------------
# 3. Fetch metadata from GenBank
#-------------------------
fetch_metadata <- function(accession) {
  
  gb <- entrez_fetch(
    db = "nuccore",
    id = accession,
    rettype = "gb",
    retmode = "text"
  )
  
  latlon <- str_extract(gb, '/lat_lon="[^"]+"')
  date   <- str_extract(gb, '/collection_date="[^"]+"')
  auth   <- str_extract(gb, "AUTHORS.+")
  
  if (!is.na(latlon)) {
    latlon <- gsub('/lat_lon="|"', "", latlon)
    parts  <- strsplit(latlon, " ")[[1]]
    lat <- parse_lat_lon(paste(parts[1:2], collapse=" "))
    lon <- parse_lat_lon(paste(parts[3:4], collapse=" "))
  } else {
    lat <- lon <- NA_real_
  }
  
  tibble(
    latitude  = lat,
    longitude = lon,
    date      = ifelse(is.na(date), NA, gsub('/collection_date="|"', "", date)),
    authors   = ifelse(is.na(auth), NA, gsub("AUTHORS\\s+", "", auth))
  )
}

#-------------------------
# 4. Search NCBI (mitochondrial only)
#-------------------------
species <- "Bathynomus"
n_records <- 20

search <- entrez_search(
  db   = "nucleotide",
  term = paste0(species, "[Organism] AND mitochondrion[Filter]"),
  retmax = n_records
)

ids <- search$ids
message("Retrieved ", length(ids), " mitochondrial records")

#-------------------------
# 5. Fetch records
#-------------------------
records <- lapply(ids, function(id) {
  
  Sys.sleep(0.3)  # polite rate limiting
  
  fasta <- entrez_fetch(
    db = "nucleotide",
    id = id,
    rettype = "fasta"
  )
  
  header <- strsplit(fasta, "\n")[[1]][1]
  header <- gsub("^>", "", header)
  
  acc <- strsplit(header, " ")[[1]][1]
  sp  <- paste(strsplit(header, " ")[[1]][2:3], collapse=" ")
  
  gb <- entrez_fetch(
    db = "nucleotide",
    id = id,
    rettype = "gbc",
    retmode = "text"
  )
  
  locus <- extract_locus(gb)
  meta  <- fetch_metadata(acc)
  
  tibble(
    accession_number = acc,
    species_name     = sp,
    locus            = locus,
    sequence_fasta   = fasta,
    latitude         = meta$latitude,
    longitude        = meta$longitude,
    date             = meta$date,
    study_ref        = meta$authors,
    seq_ref          = "NCBI GenBank",
    url              = paste0("https://www.ncbi.nlm.nih.gov/nuccore/", acc)
  )
})

df <- bind_rows(records)

#-------------------------
# 6. Write output
#-------------------------
# dir.create("data", showWarnings = FALSE)

write.csv(
  df,
  "../data/ncbi_bathynomus_mito_demo.csv",
  row.names = FALSE
)

message("Demo sourced in: data/ncbi_bathynomus_mito_demo.csv")
