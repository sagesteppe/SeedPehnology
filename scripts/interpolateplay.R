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

ggplot() + 
  geom_sf(data = clusts, aes(color = ClusterID)) 

no_cores <- parallel::detectCores()
out <- initiation_cessation(clusts, n_cores = no_cores)

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
