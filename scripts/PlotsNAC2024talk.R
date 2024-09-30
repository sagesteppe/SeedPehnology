library(tidyverse)
library(sf)

modelled <- list.files('../data/processed/timestamps/')

t <- read.csv('../data/SOS/Scouting.csv') %>% 
  select(ABBREV_NAME, taxa, scoutDate, percentDormant:percentPost, LATITUDE_DECIMAL, LONGITUDE_DECIMAL) |>
  st_as_sf(coords = c(x = 'LONGITUDE_DECIMAL', y = 'LATITUDE_DECIMAL'), crs = 4326) |>
  mutate(taxa = str_replace_all(taxa, ' ', '_')) |>
  mutate(
    taxa = case_when(
    taxa == 'Machaeranthera_canescens' ~ 'Dieteria_canescens', 
    taxa != 'Machaeranthera_canescens' ~ taxa
  ),
  scoutDate = lubridate::yday( as.Date(scoutDate, tryFormats = c('%m/%d/%Y')))) |>
  filter(taxa %in% modelled) %>%
  # we will remove Phacelia crenulata and Eriogonum fusiforme which did not flower this year #
  filter(!taxa %in% c('Phacelia_crenulata', 'Eriogonum_fusiforme')) 

# load and extract raster values to ground verified point.

extractPredictions <- function(x){
  
  taxon <- sf::st_drop_geometry(x$taxa[1])
  fp <- paste0('../data/processed/timestamps/', taxon, '/doy_preds')
  doys <- data.frame(
    doy = gsub('.//tif', '', list.files(fp)), 
    path = fp
  )

  # each species had its own biweekly intervals generated, these do not align. 
  # we will add the nearest DOY prediction to each file. If the nearest date is #
  # before the start, or after the end. We will convert these into predictions of 
  # 0 for the phenophase. 
  
  file <- 
    dplyr::inner_join(x, join_by(scoutDate == doy, closest(scoutDate == doy))) |>
    dplyr::mutate(fp = dplyr::if_else(abs(scoutDate - doy) >= 14, 0, fp))
   
  
}

?data.table::roll

?closest
