library(sf)
library(tidyverse)
library(terra)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')




# now add spatial attributes to data. These will be used to look up the 
# relevant independent variables
p <- '../data/processed/high_priority_sheets'
tab_data <- st_read(file.path(p, (paste0( basename(p), '.shp'))), quiet = T) %>% 
  filter(scntfcn == 'achnatherum hymenoides') %>% 
  dplyr::select(siteID, scientificname = scntfcn, doy)

achy <- left_join(achy, tab_data) |>
  st_as_sf() |>
  st_transform(4326)

rm(tab_data)
# we now have Day of year, and latitude. DOY will be a fixed predictor, and I anticipate
# will hold the most explanatory power, where other terms will be covariates for it



# now we can import variables which we believe correlate with flowering. 

p <- '../results/spatial/gddPCA.tif'
preds <- rast(p)

achy <- terra::extract(preds, text_result, bind = TRUE) |>
  st_as_sf()


## we don't have much scored data, so we will impute NA bulk density data ##
achy <- achy |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
         across(Pct_Bud:Pct_Dropped, ~ ifelse(.x > 0, 1, .x)))


set.seed(28)
ob <- modeller(achy)



# play with prediction

pred_df <- expand.grid(doy = seq(min(achy$doy), max(achy$doy)), 
                       bulk_density = seq(min(achy$bulk_density), max(achy$bulk_density), by = 1))

pred_df <- expand.grid(doy = seq(50, 200), 
                       bulk_density = mean(achy$bulk_density)) # toy around to see results

pred <- predict(mod.aspatial, pred_df, type = 'response', se = TRUE)

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
