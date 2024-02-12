library(sf)
library(tidyverse)
library(terra)
library(mgcv)

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
lapply(X = splicies[40:50], modeller) 

f <- paste0('../results/models/', list.files('../results/models/'))
model <- readRDS(f[])

gam.check(model)
summary(model)
plot(model)
mat <- splicies[[3]] |> st_transform(5070)


f <- file.path('../results/models', list.files('../results/models'))
model <- readRDS(f[16])

plot(model)
abline(h=0)
gam.check(model)
summary(model)

# play with prediction

# 1) make an aspatial prediction, determine start dates and end dates. 
pred_df <- expand.grid(doy = seq(min(spp$doy), max(spp$doy), length.out = 25), 
                       latitude = seq(min(spp$latitude), max(spp$latitude), length.out = 20), 
                       bulk_density = seq(min(spp$bulk_density), max(spp$bulk_density), length.out = 20))
pred_df$fit <- pred_df$fit <- predict(model, newdata = pred_df, type = 'response', se = F)
pred_df$fit <- 1- pred_df$fit

 # predict between these days onwards
if({lowerDOY <- min(pred_df[ pred_df$fit > 0.1, ]$doy)} < 0){lowerDOY <- 0} else {
  lowerDOY <- floor(lowerDOY)}
if({upperDOY <- max(pred_df[ pred_df$fit > 0.1, ]$doy)} > 365){upperDOY <- 365} else {
  upperDOY <- ceiling(upperDOY)}

# predict only on these days. 
timeStamps <- round(seq(lowerDOY, upperDOY, by = 14)) # biweekly time stamps for prediction

# create raster layers for each time point. - this is memory hungry, so each raster
# will be written to disk, and then these will be reloaded. 

r1 <- rast(preds)[[1]]
# s <- rast(c(replicate(length(timeStamps), r1)))

for (i in seq_along(1:length(timeStamps))){
  r_fill <- setValues(r1, timeStamps[i])
  r_fill <- terra::mask(r_fill, preds[[1]])
  writeRaster(r_fill, )
  
}

plot(s)
names(r1) <- 'doy'

ts1 <- c(preds, r1)
doy149 <- terra::predict(ts1, model, type="response") 
doy149 <- app(doy149, function(x){1-x})


plot(doy149)
points(spp)
# rip out the terms here to subset raster stack... 
acle_model$gam$formula  # do this

# subset the prediction surface to areas with suitable habitat

# extract the relevant variables for each occupied pixel

# create prediction data.frame using observed values and full DOY Range for each cell

# identify first day with > 5% probability of flowering, peak, and last day with >10% flowering

# write 9 rasters, 1) FIRST , 2) PEAK, 3) LAST, and each with it's upper and lower SE. 
model

spp <- splicies[[which( gsub('.shp', '', basename(f[25])) == 
                          gsub(' ', '_', names(splicies))  )]]

f[25] # bidens frondosa. 

predictWRAP <- function(x, y){
  
  
}

ob145 <- terra::predict(preds, model, const = data.frame(doy = 145), type="response") 
ob152 <- terra::predict(preds, model, const = data.frame(doy = 152), type="response") 
ob159 <- terra::predict(preds, model, const = data.frame(doy = 159), type="response") 
ob166 <- terra::predict(preds, model, const = data.frame(doy = 166), type="response") 
ob225 <- terra::predict(preds, model, const = data.frame(doy = 225), type="response") 
ob280 <- terra::predict(preds, model, const = data.frame(doy = 280), type="response") 


ob_lyrs <-  terra::predict(preds, model, 
                           const = data.frame(doy = c(145, 159, 172)), type="response") 

plot(ob100) # these should all be 0
plot(ob225) # these should all be 1
plot(ob310) # these should all be 0

ob145_inv <- app(ob145, function(x){1-x})
ob152_inv <- app(ob152, function(x){1-x})
ob159_inv <- app(ob159, function(x){1-x})
ob166_inv <- app(ob166, function(x){1-x})
ob225_inv <- app(ob225, function(x){1-x})
ob280_inv <- app(ob280, function(x){1-x})

plot(ob145_inv)
plot(ob152_inv)
plot(ob159_inv)
plot(ob166_inv)
plot(ob225_inv)

# out <- 1- out 

plot(out)
points(vect(spp))


ggplot(pred_df, aes(x = doy, y = fit)) + 
  geom_point() + 
  geom_smooth(se = FALSE) + 
  geom_smooth(aes(y = SE_low), se = FALSE, color = 'red') + 
  geom_smooth(aes(y = SE_high), se = FALSE, color = 'red') + 
  theme_classic()

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!


