library(sf)
library(tidyverse)
library(terra)
library(mgcv)

setwd('/media/steppe/hdd/SeedPhenology/scripts')
#setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

f <- paste0('../results/PresAbs/', list.files('../results/PresAbs', pattern = 'shp$'))
spp <- lapply(f, st_read, quiet = TRUE) |>
  bind_rows() |>
  select(-accessr) |>  
  pivot_longer(cols = doy:cessAbs, values_to = 'doy', names_to = 'flowering') |> 
  mutate(flowering = if_else(flowering == 'doy', 1, 0)) # don't convert to factor

# now we can import variables which we believe correlate with flowering. 
p <- '../data/spatial/processed/'
preds <- rast(lapply(paste0(p, list.files(p, pattern = '.tif')), rast))
names(preds)[19] <- 'latitude'

vif_r <- usdm::vifstep(spatSample(preds, 5000))
preds <- subset(preds, vif_r@results[['Variables']])
rm(vif_r)

spp <- terra::extract(preds, spp, bind = TRUE) |> 
  st_as_sf() |> 
  mutate( 
    across( 
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) 

splicies <- split(spp, f = spp$scntfcnm)
# lapply(X = splicies, modeller) 

f1 <- file.path('../results/models', list.files('../results/models'))

# play with prediction
splicies <- bind_rows(splicies)

sdm_surfs <- list.files('../data/spatial/PhenPredSurfaces') # need to align rasters
focal_surf <- terra::rast(
  file.path( '../data/spatial/PhenPredSurfaces', sdm_surfs[1])
)
preds <- project(preds, crs(focal_surf))

remove <- setdiff(gsub('[.]rds', '', basename(f1)), gsub('[.]tif', '', sdm_surfs))
f1 <- f1[ ! grepl(paste(remove, collapse = '|'), f1) ]
splicies <- bind_rows(splicies)
lapply(f1, spat_predict, spp = splicies)

# lomatium grayi, peak is less than initiation
# penstemon eatonii, not all cells covered

# erigeron bloomeri, not all cells recovered for initiation and peak
# euthamnia occidentalis, not all cells recovered for initiation and peak

f <- list.files('../data/processed/timestamps', recursive = T)
f <- paste0('../data/processed/timestamps/',  f[grep('doy_preds', f)])
f <- paste0(unique(dirname(f)), '/')

lapply(f, spat_summarize)

plot(flr_events, col = c('red', 'yellow', 'green', 'purple', 'brown', 'cyan', 'black'))


p <- '../data/processed/timestamps'
acsp <- rast(file.path(p, 'Lasthenia_gracilis/summary_doys/Lasthenia_gracilis.tif'))
plot(acsp, col = 'red')
plot(acsp, col = c('red', 'yellow', 'green', 'purple', 'brown', 'cyan', 'black'))
