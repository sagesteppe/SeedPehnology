library(sf)
library(tidyverse)
library(terra)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

aggl <- st_read('../results/PresAbs/Achnatherum_lemmonii.shp', quiet = TRUE) |>
  select(-accessr) %>% 
  pivot_longer(cols = doy:cessAbs, values_to = 'doy', names_to = 'flowering') |>
  mutate(flowering = if_else(flowering == 'doy', 1, 0)) # don't convert to factor

# now we can import variables which we believe correlate with flowering. 
p <- '../data/spatial/processed/'
preds <- rast(lapply(paste0(p, list.files(p, pattern = '.tif')), rast))
names(preds)[19] <- 'latitude'

vif_r <- usdm::vifstep(spatSample(preds, 5000))
preds <- subset(preds, vif_r@results[['Variables']])
rm(vif_r)

aggl <- terra::extract(preds, aggl, bind = TRUE) |>
  st_as_sf() |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))

rm(p)
modeller(aggl) # 10:13

acle_model <- readRDS('../results/models/Achnatherum_lemmonii.rds')
aggl_model <- readRDS('../results/models/Agoseris_glauca.rds')
acam_model <- readRDS('../results/models/Acmispon_americanus.rds')

concurvity(acle_model)
summary(acle_model)
out <- gam.check(acle_model)

check_table <- data.frame(mgcv:::k.check(acam_model))
rownames(check_table[check_table[,4] < 0.05,]) # this variable requires a higher k'

# by default 10 k are used. We will increase the value by x1.5

update(acam_model, )

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


