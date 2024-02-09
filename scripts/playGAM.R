library(sf)
library(tidyverse)
library(terra)

setwd('~/Documents/SeedPhenology/scripts')
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
splicies <- Filter(function(y) nrow(y) >= 50, splicies)
lapply(X = splicies[2], modeller) # 2:40 ~ 500 occurrences

































f <- paste0('../results/models/', list.files('../results/models/'))
model <- readRDS(f[2])

gam.check(model)
summary(model)
plot(model)
mat <- splicies[[3]] |> st_transform(5070)








plot(model)
gam.check(model)
summary(model)

data.frame(mgcv:::k.check(model))
rownames(check_table[check_table[,4] < 0.05,]) # this variable requires a higher k'
names(splicies)

# play with prediction


# rip out the terms here to subset raster stack... 
acle_model$gam$formula  # do this

# subset the prediction surface to areas with suitable habitat

# extract the relevant variables for each occupied pixel

# create prediction data.frame using observed values and full DOY Range for each cell

# identify first day with > 5% probability of flowering, peak, and last day with >10% flowering

# write 9 rasters, 1) FIRST , 2) PEAK, 3) LAST, and each with it's upper and lower SE. 








pred_df <- expand.grid(doy = seq(min(acle$doy), max(acle$doy), length.out = 25), 
                       vpd_mean = seq(min(acle$vpd_mean), max(acle$vpd_mean), length.out = 20), 
                       gdgfgd5 = seq(min(acle$gdgfgd5), max(acle$gdgfgd5), length.out = 20))

pred <- predict(acle_model, newdata = pred_df, type = 'response', se = TRUE)
pred_df$fit <- pred$fit
pred_df$SE_low <- pred$fit - pred$se.fit
pred_df$SE_high <- pred$fit + pred$se.fit


pred_core <- pred_df[pred_df$fit > 0.2,]
start <- pred_core[ which.min(pred_core$doy), 'doy']
end <- pred_core[ which.max(pred_core$doy), 'doy']
peak <- pred_df[ which.max(pred_df$fit), 'doy'] # identify peak flower, day 131. 

ggplot(pred_df, aes(x = doy, y = fit)) + 
  geom_point() + 
  geom_smooth(se = FALSE) + 
  geom_smooth(aes(y = SE_low), se = FALSE, color = 'red') + 
  geom_smooth(aes(y = SE_high), se = FALSE, color = 'red') + 
  theme_classic()

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!


