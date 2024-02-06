library(sf)
library(tidyverse)
library(terra)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

acle <- st_read('../results/PresAbs/Achnatherum_lemmonii.shp', quiet = TRUE) |>
  select(-accessr) %>% 
  pivot_longer(cols = doy:cessAbs, values_to = 'doy', names_to = 'flowering') |>
  mutate(flowering = as.factor(if_else(flowering == 'doy', 1, 0)))


# now we can import variables which we believe correlate with flowering. 
p <- '../data/spatial/processed/'
preds <- rast(lapply(paste0(p, list.files(p, pattern = '.tif')), rast))
names(preds)[19] <- 'latitude'

acle <- terra::extract(preds, acle, bind = TRUE) |>
  st_as_sf() |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))

acle1 <- acle[, 1:5]
ob <- modeller(acle)



























# play with prediction

pred_df <- expand.grid(doy = seq(min(achy$doy), max(achy$doy)), 
                       bulk_density = seq(min(achy$bulk_density), max(achy$bulk_density), by = 1))

pred_df <- expand.grid(doy = seq(50, 200), 
                       bulk_density = mean(achy$bulk_density)) # toy around to see results

pred <- predict(mod.aspatial, pred_df, type = 'response', se = TRUE)

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
