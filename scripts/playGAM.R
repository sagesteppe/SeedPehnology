library(sf)
library(tidyverse)
library(terra)
library(mgcv)

setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

f <- paste0('../results/PresAbs/', list.files('../results/PresAbs', pattern = 'shp$'))
spp <- lapply(f, st_read, quiet = TRUE) |>
  bind_rows() |>
  select(-accessr) #|>  
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

spat_predict <- function(x){
  
  # identify the taxon we are working with. 
  taxon <- gsub('.rds', '', basename(f[41]))
  
  #identify independent variables
  terms <- unlist(strsplit(split = ' [+] ', as.character(model$terms[[3]][[2]]))) 
  terms <- terms[ grep('[+]', terms, invert = TRUE)]
  terms <- c(terms, 'doy')
  spp_range <- sf::st_drop_geometry(spp[spp$flowering==1,terms])

  #create prediction grid of variables. 
  dfParameterValues <- data.frame(
    ParameterName = colnames(spp_range),
    seqFrom = apply(spp_range, MARGIN = 2, FUN = min),
    seqTo = apply(spp_range, MARGIN = 2, FUN = max), 
    lOut = c(15, 15, 15))
  
  pred_df <- setNames(
    expand.grid(
      data.frame( 
        apply(dfParameterValues[,c("seqFrom", "seqTo", 'lOut')], 1,
              function(x) seq(from = x['seqFrom'], to = x['seqTo'], length.out = x['lOut']))
      )
    ), dfParameterValues$ParameterName
  )
  
  # fit model
   pred_df$fit <- predict(model, newdata = pred_df, type = 'response', se = F)
  
  # identify first day with > 5% probability of flowering, peak, and last day with >10% flowering
  # predict between these days onwards
  if({lowerDOY <- min(pred_df[pred_df$fit > 0.55, ]$doy) } < 0){lowerDOY <- 0} else {
    lowerDOY <- floor(lowerDOY)}
  if({upperDOY <- max(pred_df[ pred_df$fit > 0.55, ]$doy) } > 365){upperDOY <- 365} else { 
    upperDOY <- ceiling(upperDOY)}
   
   # predict only on these days. 
   timeStamps <- round(seq(lowerDOY, upperDOY, by = 14)) # biweekly time stamps for prediction
   
   # predict only in areas with suitable habitat for the species. 
   sdm_surfs <- list.files('../data/raw/PhenPredSurfaces')
   focal_surf <- terra::rast( sdm_surfs[ grep(taxon, sdm_surfs) ] )
   preds_sub <- terra::mask(preds, focal_surf)
   
   # create raster layers for each time point. - this is memory hungry, so each raster
   # will be written to disk, and then these will be reloaded. 
   r1 <- rast(preds_sub)[[1]]
   names(r1) <- 'doy'
   p1 <- file.path('../data/processed/timestamps', taxon)
   
   dir.create(file.path(p1, 'doy_constants'), showWarnings = FALSE, recursive = T)
   for (i in seq_along(1:length(timeStamps))){
     
     r_fill <- terra::setValues(r1, timeStamps[i])
     r_fill <- terra::mask(r_fill, preds_sub[[1]])
     
     terra::writeRaster(r_fill, paste0(p1, '/doy_constants/', timeStamps[[i]],'.tif'), overwrite = T)
     rm(r_fill) 
   }
   
   doy_stack <- orderLoad(paste0(p1, '/doy_constants/'))
   
   # predict the probability of the species flowering at each time point
   dir.create(file.path(p1, 'doy_preds'), showWarnings = FALSE)
   for (i in seq_along(1:length(timeStamps))){
     
     space_time <- c(preds_sub, doy_stack[[i]])
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
   
}

# isolate the initiation of flowering (10%), peak (max value), and cessation (90%) for each cell. 

p1 <- file.path('../data/processed/timestamps', 'Dichelostemma_capitatum')
pred_stack <- orderLoad(paste0(p1, '/doy_preds/'))
plot(pred_stack)
out <- values(pred_stack)

# this identifies a 'peak' date. 
pred_stack <- ifel(is.na(pred_stack), -999, pred_stack) # we need to remove ocean NA's
p1 <- terra::app(pred_stack, which.max) # peak flower. 
p2 <- terra::subst(p1, from = 1:dim(pred_stack)[3], to = names(pred_stack)) 

plot(pred_stack)
plot(p2, col = c('red', 'blue', 'green', 'purple', 'brown'))

# this isolates an effective start date.  
pred_stack <- ifel(pred_stack < 0, 999, pred_stack)
p3 <- terra::app(pred_stack, which.min) # start date.  
p4 <- terra::subst(p3, from = 1:dim(pred_stack)[3], to = names(pred_stack)) 
plot(p4)




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


