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
# lapply(X = splicies[40:50], modeller) 



f <- file.path('../results/models', list.files('../results/models'))
model <- readRDS(f[41])

plot(model)
abline(h=0)
gam.check(model)
summary(model)

spp <- splicies[[which( gsub('.rds', '', basename(f[41])) == 
                          gsub(' ', '_', names(splicies))  )]]

# play with prediction

ob <- model$terms[[3]]

spat_predict <- function(x){
  
  spp_range <- spp[spp$flowering==1,]
  hist(spp_range$doy)
  
  pred_df <- expand.grid(
    doy = seq(min(spp_range$doy)-14, max(spp_range$doy)+14, length.out = 25), # fixed, always here. 
    bio10 = seq(min(spp_range$bio10), max(spp_range$bio10)), 
    gddlgd10 = seq(min(spp_range$gddlgd10), max(spp_range$gddlgd10), length.out = 20)
    )
  
  pred_df$fit <- pred_df$fit <- predict(model, newdata = pred_df, type = 'response', se = F)
  plot(pred_df$doy, pred_df$fit)
  
}

# 1) make an aspatial prediction, determine start dates and end dates. 
spp_range <- spp[spp$flowering==1,]
hist(spp_range$doy)
pred_df <- expand.grid(doy = seq(min(spp_range$doy)-14, max(spp_range$doy)+14, length.out = 25), 
                       bio10 = seq(min(spp_range$bio10), max(spp_range$bio10)), 
                       gddlgd10 = seq(min(spp_range$gddlgd10), max(spp_range$gddlgd10), length.out = 20))

pred_df$fit <- pred_df$fit <- predict(model, newdata = pred_df, type = 'response', se = F)
plot(pred_df$doy, pred_df$fit)

# identify first day with > 5% probability of flowering, peak, and last day with >10% flowering
# predict between these days onwards
if({lowerDOY <- min(pred_df[pred_df$fit > 0.1, ]$doy) } < 0){lowerDOY <- 0} else {
  lowerDOY <- floor(lowerDOY)}
if({upperDOY <- max(pred_df[ pred_df$fit > 0.1, ]$doy) } > 365){upperDOY <- 365} else { 
  upperDOY <- ceiling(upperDOY)}

# predict only on these days. 
timeStamps <- round(seq(lowerDOY, upperDOY, by = 14)) # biweekly time stamps for prediction

# create raster layers for each time point. - this is memory hungry, so each raster
# will be written to disk, and then these will be reloaded. 

r1 <- rast(preds)[[1]]
names(r1) <- 'doy'
taxon <- gsub('.rds', '', basename(f[41]))
p1 <- file.path('../data/processed/timestamps', taxon)

dir.create(file.path(p1, 'doy_constants'), showWarnings = FALSE, recursive = T)
for (i in seq_along(1:length(timeStamps))){
  
  r_fill <- terra::setValues(r1, timeStamps[i])
  r_fill <- terra::mask(r_fill, preds[[1]])
  
  terra::writeRaster(r_fill, paste0(p1, '/doy_constants/', timeStamps[[i]],'.tif'), overwrite = T)
  rm(r_fill) 
}

# reload the rasters into a stack - ensure they are in ascending DOY Order. 

orderLoad <- function(path){
  doy_f <- list.files(path)
  doy_stack <- paste0(path, 
                      sort(as.numeric(gsub('.tif', '', 
                                             basename(doy_f)))), '.tif')
  doy_stack <- terra::rast(lapply(doy_stack, terra::rast))
  return(doy_stack)
}

doy_stack <- orderLoad(paste0(p1, '/doy_constants/'))

# predict the probability of the species flowering at each time point
dir.create(file.path(p1, 'doy_preds'), showWarnings = FALSE)
for (i in seq_along(1:length(timeStamps))){
  
  space_time <- c(preds, doy_stack[[i]])
  time_pred <- terra::predict(space_time, model, type="response") 
 # time_pred <- terra::app(time_pred, function(x){1-x})
  
  names(time_pred) <- timeStamps[[i]]
  msk <- terra::ifel(time_pred < 0.01, NA, time_pred)
  time_pred <- terra::mask(time_pred, msk)
  
  NACount <- freq(is.na(time_pred))
  NACount <- NACount[ NACount$value==1, 'count'] / (dim(time_pred)[1]*dim(time_pred)[2])
  
  if(NACount < 0.99){
  terra::writeRaster(time_pred, 
                     paste0(p1, '/doy_preds/', timeStamps[[i]],'.tif'), overwrite = T)
  } else {message('> 99% of Cells are NA, not writing layer ',  i, ' to disk.')}
  
  rm(time_pred)
  terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
}

pred_stack <- orderLoad(paste0(p1, '/doy_preds/'))
plot(pred_stack)

# isolate the initiation of flowering (10%), peak (max value), and cessation (90%) for each cell. 

# subset the prediction surface to areas with suitable habitat
















# write 9 rasters, 1) FIRST , 2) PEAK, 3) LAST, and each with it's upper and lower SE. 
model


f[25] # bidens frondosa. 



ggplot(pred_df, aes(x = doy, y = fit)) + 
  geom_point() + 
  geom_smooth(se = FALSE) + 
  geom_smooth(aes(y = SE_low), se = FALSE, color = 'red') + 
  geom_smooth(aes(y = SE_high), se = FALSE, color = 'red') + 
  theme_classic()

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!


