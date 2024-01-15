library(sf)
library(tidyverse)
library(terra)
library(rgeoda)

setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

acle <- ridigbio::idig_search(rq = list(scientificname = 'Achnatherum lemmonii')) |> 
  dplyr::select(uuid, scientificname, datecollected, collector, geopoint.lon, geopoint.lat) |> 
  tidyr::drop_na(geopoint.lon, geopoint.lat, datecollected) |> 
  dplyr::mutate(datecollected = as.Date(datecollected), 
                doy = as.numeric(lubridate::yday(datecollected))) |> 
  dplyr::filter(datecollected > '1979-01-01') |> # sensing data comes on line here
  dplyr::arrange(datecollected) |>
  sf::st_as_sf(coords = c(x = 'geopoint.lon', y = 'geopoint.lat'), crs = 4326)

hist(acle$doy)

# extract the PCA axis to the values for clustering
pca1 <- rast('../results/spatial/gddPCA.tif')

acle <- extract(pca1, acle, bind = TRUE) %>% 
  st_as_sf()

# select points to cluster which are actually in meaningful geographic proximity /
# the points futher out will not require flowering absence records. 
nf <- acle[ st_nearest_feature(acle), ]
acle2clust <- acle[ as.numeric( st_distance(acle, nf, by_element = T) ) < 80000 , ] 

clusts <- geoClust_wrap(acle2clust, 'PC1')
rm(nf, acle2clust)

no_cores <- parallel::detectCores()
out <- initiation_cessation(clusts, n_cores = no_cores)


grps <- split(clusts, f = clusts$ClusterID)
estimates <- split(out, f = out$ClusterID)

lapply(estimates, doy_gen, event = 'cessation')

#' create pseudoabsences for phenophases
#' @param x the output of geoClust
#' @param estimates the output of initiation_cessation
pheno_abs <- function(x, estimates){
  
  grps <- split(x, f = x$ClusterID)
  estimates <- split(estimates, f = estimates$ClusterID)
  
  ed <- lapply(estimates, doy_gen, 'initiation')
  ld <- lapply(estimates, doy_gen, 'cessation')
  
  # PCA1 1 = hottest, 0 = coldest. 
  # Give the earliest initiation DOY absence to the warmest place, 
  # and the latest absences to the coldest places. 
   sort(PC1, decreasing = TRUE)
  
  #  Give the earliest cessation DOY absence to the warmest place, 
   # and the latest absences to the coldest place. 
   sort(PC1, decreasing = TRUE)
   
   # identify records in clusters with >= 1/3 of their neighbors more near another
   # cluster than their own. These records will have their original values wiped, 
   # and their values will be reassigned via interpolation of the PCA axis. 
   
   st_distance()[1:3]
   
   early_mod <- fields::Tps(x$PC1, x$doy_ed)
   late_mod <- fields::Tps(x$PC1, x$doy_ld)
   
   # subset the raster surface to only cover the points which could benefit
   # from interpolated values. 
   hull <- sf::st_convex_hull(x)
   surface <- terra::mask(pca1, hull)
   early_surf <- terra::interpolate(surface, early_mod)
   late_surf <- terra::interpolate(surface, late_mod)
   
}
?interpolate

??tps
fit <- Krig(ChicagoO3$x, ChicagoO3$y, aRange=20)  


# generate absences for 50% of the sites per cluster. 
# sample the 50% of absences and order them along the PCA variable
# to assign to points. 
# Interpolate the values for the remaining 50% of sites. 


# try another species. 
abvi <- ridigbio::idig_search(rq = list(scientificname = 'Abronia villosa')) |> 
  dplyr::select(uuid, scientificname, datecollected, collector, geopoint.lon, geopoint.lat) |> 
  tidyr::drop_na(geopoint.lon, geopoint.lat, datecollected) |> 
  dplyr::mutate(datecollected = as.Date(datecollected), 
                doy = as.numeric(lubridate::yday(datecollected))) |> 
  dplyr::filter(datecollected > '1979-01-01') |> # sensing data comes on line here
  dplyr::arrange(datecollected) |>
  dplyr::filter(geopoint.lat < 40) |>
  sf::st_as_sf(coords = c(x = 'geopoint.lon', y = 'geopoint.lat'), crs = 4326)

abvi <- extract(pca1, abvi, bind = TRUE) %>% 
  st_as_sf()
nf <- abvi[ st_nearest_feature(abvi), ]
abvi2clust <- abvi[ as.numeric( st_distance(abvi, nf, by_element = T) ) < 80000 , ] 



# note some species may begin flowering in upper DOY, and peak afterwards....
# see for example abronia villosa. 
hist(abvi$doy)
rm(abvi)s

 


## Perhaps some species are bimodal in flowering dates!!! 
# we can test bimodal distributions using a fn from LaplacesDemon 
LaplacesDemon::is.bimodal(data)
