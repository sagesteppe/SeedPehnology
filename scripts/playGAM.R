library(sf)
library(tidyverse)
library(terra)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

acle <- st_read('../results/PresAbs/Achnatherum_lemmonii.shp', quiet = TRUE) |>
  select(-accessr) %>% 
  pivot_longer(cols = doy:cessAbs, values_to = 'doy', names_to = 'flowering') |>
  mutate(flowering = if_else(flowering == 'doy', 1, 0)) # don't convert to factor

# now we can import variables which we believe correlate with flowering. 
p <- '../data/spatial/processed/'
preds <- rast(lapply(paste0(p, list.files(p, pattern = '.tif')), rast))
names(preds)[19] <- 'latitude'

acle <- terra::extract(preds, acle, bind = TRUE) |>
  st_as_sf() |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))

modeller(acle)

acle_model <- readRDS('../results/models/Achnatherum_lemmonii.rds')
# play with prediction

pred_df <- expand.grid(doy = seq(min(acle$doy), max(acle$doy)), 
                       vpd_mean = seq(min(acle$vpd_mean), max(acle$vpd_mean), by = 1), 
                       gdgfgd5 = seq(min(acle$gdgfgd5), max(acle$gdgfgd5), by = 1))

pred_df <- expand.grid(doy = seq(50, 200), 
                       bulk_density = mean(achy$bulk_density)) # toy around to see results

pred <- predict(acle_model, pred_df, type = 'response', se = TRUE)

plot(pred)
# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
