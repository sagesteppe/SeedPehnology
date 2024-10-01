library(tidyverse)
library(sf)

modelled <- list.files('../data/processed/timestamps/')

t <- read.csv('../data/SOS/Scouting.csv') %>% 
  select(taxa, scoutDate, percentDormant:percentPost, LATITUDE_DECIMAL, LONGITUDE_DECIMAL) |>
  st_as_sf(coords = c(x = 'LONGITUDE_DECIMAL', y = 'LATITUDE_DECIMAL'), crs = 4326) |>
  st_transform(5070) |>
  mutate(taxa = str_replace_all(taxa, ' ', '_')) |>
  mutate(
    taxa = case_when(
    taxa == 'Machaeranthera_canescens' ~ 'Dieteria_canescens', 
    taxa != 'Machaeranthera_canescens' ~ taxa
  ),
  scoutDate = lubridate::yday( as.Date(scoutDate, tryFormats = c('%m/%d/%Y'))),
  sDate_join = as.numeric(scoutDate), .before = geometry) |>
  filter(taxa %in% modelled) %>%
  # we will remove Phacelia crenulata and Eriogonum fusiforme which did not flower this year #
  filter(!taxa %in% c('Phacelia_crenulata', 'Eriogonum_fusiforme', 'Eriogonum_shockleyi')) 

# load and extract raster values to ground verified point.

extractPredictions <- function(x){
  
  taxon <- sf::st_drop_geometry(x$taxa[1])
  fp <- paste0('../data/processed/timestamps/', taxon, '/doy_preds')
  doys <- data.frame(
    sDate_join = as.numeric(gsub('.tif', '', list.files(fp))), 
    path = paste0(fp, '/', list.files(fp))
  )

  # each species had its own biweekly intervals generated, these do not align. 
  # we will add the nearest DOY prediction to each file. If the nearest date is #
  # before the start, or after the end. We will convert these into predictions of 
  # 0 for the phenophase. 
  
  data.table::setDT(doys); data.table::setkey(doys, sDate_join)
  data.table::setDT(x); data.table::setkey(x, sDate_join)
  set_merged <- doys[x, roll = "nearest"]
  
  set_merged <- set_merged |> 
    dplyr::mutate(
      sDate_join = as.numeric(gsub('.tif', '', basename(path))),
      path = dplyr::if_else(abs(scoutDate - sDate_join) >= 14, NA, path)) |>
    sf::st_as_sf()
  
  # now split each of these objects by the raster for the date we will be comparing them to. 
  # we can load the raster once, and compare it to several scouting points simultaneously. 
  
  NAindexes <- which(is.na(set_merged$path))

  if(length(NAindexes)>0){
    noRasters <- set_merged[NAindexes,]
    set_merged <- set_merged[-NAindexes,]
  }
  
  # and now run the process for the populated cells. 
  obbies <- split(set_merged, f = set_merged$sDate_join)

  extractByDate <- function(x){
    
    Path <- sf::st_drop_geometry(x) |> 
      dplyr::pull(path) |> 
      unique() 

    surface <- terra::rast(Path)
    extracted_values <- terra::extract(surface, terra::vect(x)) |>
      dplyr::select(-ID) 
    
    extracted_values <- data.frame(
      doy = colnames(extracted_values),
      prob = extracted_values[,1]
    )
    
    
    # combbine the extracted values with the features of the input data set.
    x <- dplyr::bind_cols(x, extracted_values)
    
    # if any of the rows are NA, use a feature engineered probability of flowering for them
    
    if(any(is.na(x$prob))){

      # now proceed row wise 
      need_imputed <- x[which(is.na(x$prob)),]
      
      for(i in 1:nrow(need_imputed)){
      # mask any raster sites > 10k from population #
      mask <- sf::st_buffer(need_imputed[i,], 10000)
      surface_m <- terra::mask(surface, mask, inverse = FALSE)
      m <- terra::global(surface_m, mean, na.rm = TRUE)
      need_imputed[i,'prob'] <- m 
      }
      
      x <- dplyr::bind_rows(
        x[which(!is.na(x$prob)),],
        need_imputed)
      
    }
      return(x)
  }
  
  
  ebd <- lapply(obbies, extractByDate)
  ebd <- dplyr::bind_rows(ebd)
  if(exists('noRasters')){
    
    if(nrow(ebd)> 0){
      noRasters <- dplyr::bind_cols(
        noRasters, 
        data.frame(
          doy = as.character(noRasters$sDate_join), 
          prob = 999))
    ebd <- dplyr::bind_rows(ebd, noRasters) |>
      dplyr::arrange(scoutDate)} else {
      ebd <- noRasters
      }}
  
  # now make values outside of the range 0-1 - 0 these are our not-flowering dates. 
  return(ebd)
  
}


dat <- split(t, f = t$taxa)
ob <- lapply(dat, extractPredictions) |>
  bind_rows() 

