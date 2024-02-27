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
lapply(X = splicies[8], modeller) 

f1 <- file.path('../results/models', list.files('../results/models'))

plot(model)
abline(h=0)
gam.check(model)
summary(model)

# play with prediction
splicies <- bind_rows(splicies)
head(splicies)

sdm_surfs <- list.files('../data/spatial/PhenPredSurfaces') # need to align rasters
focal_surf <- terra::rast(
  file.path( '../data/spatial/PhenPredSurfaces', sdm_surfs[1])
)
preds <- project(preds, crs(focal_surf))

#' predict a gam into space and time
#' @param x a vector of paths to models
#' @param spp a dataframe of all species and variables which are relevant to prediction.
spat_predict <- function(x, spp){
  
  # identify the taxon we are working with. 
  taxon <- gsub('.rds', '', basename(x))
  
  # read in the fitted model. 
  model <- readRDS(x)
  
  #identify independent variables
  terms <- unlist(strsplit(split = ' [+] ', as.character(model$terms[[3]][[2]]))) 
  terms <- terms[ grep('[+]', terms, invert = TRUE)]
  if(all(terms==1)){terms <- 'doy'} else {terms <- c(terms, 'doy')}
  
  spp_f <- spp[spp$scntfcnm== gsub('_', ' ', taxon),]
  spp_range <- sf::st_drop_geometry(spp_f[spp_f$flowering==1, terms])

  #create prediction grid of variables. 
  dfParameterValues <- data.frame(
    ParameterName = colnames(spp_range),
    seqFrom = apply(spp_range, MARGIN = 2, FUN = min),
    seqTo = apply(spp_range, MARGIN = 2, FUN = max),
    lOut = rep(15, times = ncol(spp_range)))
  
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
  
  # identify first day with > 0.6% probability of flowering, and last day with >.6% flowering
  # predict between these days onwards
  if({lowerDOY <- min(pred_df[pred_df$fit > 0.55, ]$doy) } < 0){lowerDOY <- 0} else {
    lowerDOY <- floor(lowerDOY)}
  if({upperDOY <- max(pred_df[ pred_df$fit > 0.6, ]$doy) } > 365){upperDOY <- 365} else { 
    upperDOY <- ceiling(upperDOY)}
   
   # if the upperdoy is greater than 1 month after the last observed flowering record, cull it to +31 days.
   
   # predict only on these days. 
   timeStamps <- round(seq(lowerDOY, upperDOY, by = 14)) # biweekly time stamps for prediction
   
   # predict only in areas with suitable habitat for the species. 
   sdm_surfs <- list.files('../data/spatial/PhenPredSurfaces')
   focal_surf <- terra::rast(
     file.path( '../data/spatial/PhenPredSurfaces',
                sdm_surfs[ grep(taxon, sdm_surfs)])
     )

   focal_surf <- resample(focal_surf, preds)
   preds_sub <- terra::crop(preds, focal_surf, mask = TRUE)
   
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
   
   # determine how many cells are suitable habitat. # use this to determine
   # whether it's worth writing out that time stamp
   n_cells <- (dim(focal_surf)[1] * dim(focal_surf)[2])
   prop_pop <- freq(is.na(focal_surf))$count[1]  / n_cells 
   
   for (i in seq_along(1:length(timeStamps))){
     
     space_time <- c(preds_sub, doy_stack[[i]])
     time_pred <- terra::predict(space_time, model, type="response") 
     
     names(time_pred) <- timeStamps[[i]]
     msk <- terra::ifel(time_pred < 0.5, NA, time_pred)
     time_pred <- terra::mask(time_pred, msk)
     
     populated <- freq(is.na(time_pred))
     populated <- populated[ populated$value==0, 'count']
     if(length(populated)==0){populated<-1}
     
     if(populated >= ((n_cells * prop_pop) * 0.05) ){
       terra::writeRaster(time_pred, 
                          paste0(p1, '/doy_preds/', timeStamps[[i]],'.tif'), overwrite = T)
     } else {message('> 95% of Cells are NA, not writing layer ',  i, ' to disk.')}
     
     rm(time_pred)
     terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
   }
   
   unlink(paste0(p1, '/doy_constants/')) # remove the constants. 
   pred_stack <- orderLoad(paste0(p1, '/doy_preds/'))
   
   # this identifies a 'peak' date. 
   peak_date <- terra::app(pred_stack, which.max) # peak flower. 
   peak_flr <- terra::subst(peak_date, 
                            from = 1:dim(pred_stack)[3], to = names(pred_stack)) 
   
   # this isolates an effective start date.  
   start_date <- terra::app(pred_stack, which.min) # start date.  
   start_flr <- terra::subst(start_date, 
                             from = 1:dim(pred_stack)[3], to = names(pred_stack)) 
   
   # this isolates an effective end of flowering date. 
   v_mat <- terra::values(pred_stack)
   v <- colnames(v_mat)[ max.col(!is.na(v_mat), ties.method = 'last') ]
   end_flr <- terra::setValues(focal_surf, values = v)
   rm(v_mat, v)
   
   flr_events <- c(start_flr, peak_flr, end_flr)
   names(flr_events) <- c('Initiation', 'Peak', 'Cessation')
   flr_events <- terra::mask(flr_events, focal_surf)
   flr_events <- terra::trim(flr_events)
   
   dir.create(file.path(p1, 'summary_doys'), showWarnings = FALSE)
   writeRaster(flr_events, file.path(p1, 'summary_doys', paste0(taxon, '.tif')))
   terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
   
}

summary( readRDS(f1[8]) )
spat_predict(f1[127], spp = splicies)
