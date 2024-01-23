library(sf)
library(tidyverse)
library(terra)
library(rgeoda)

setwd('~/Documents/SeedPhenology/scripts')
source('functions.R')

fields2getrecord <- c(
  "data.dwc:scientificName",  "data.dwc:decimalLatitude", 
  "data.dwc:decimalLongitude", "collector", "uuid",
  "data.dwc:year", "data.dwc:month", "data.dwc:day")
names(fields2getrecord) <- c(
  'scientificname', 'geopoint.lat', 'geopoint.lon', 'collector', 'uuid', 'year', 'month', 'day')


napu <- specimen_sampler(x = 'Nassella pulchra', date = '1981-01-01') # grab records. 
napu <- rec_selec(napu) # drop dupes and far off points

# extract pca values

preds <- rast('../results/spatial/gddPCA.tif')

napu <- terra::extract(preds, napu, bind = TRUE) |>
  st_as_sf() %>% 
  drop_na(PC1)

# flag records for manual review. 
napu <- visual_review_flagger(napu)

ggplot() + 
  geom_point(data = napu, aes(x = PC1, y = doy, color = phen_flag))


# write out records here !!! and do it. 

# select points to cluster which are actually in meaningful geographic proximity /
# the points futher out will not require flowering absence records. 
nf <- acle[ st_nearest_feature(acle), ]
acle2clust <- acle[ as.numeric( st_distance(acle, nf, by_element = T) ) < 80000 , ] 
 
clusts <- geoClust_wrap(acle2clust, 'PC1') # identify clusters here. 
rm(nf, acle2clust)

no_cores <- parallel::detectCores()
ince <- initiation_cessation(clusts, n_cores = no_cores)

ll <- pheno_abs(clusts, ince)

ggplot() + 
  geom_sf(data = ll, aes(color = cessAbs))


ggplot(data = ll, aes(x = doy)) + 
  geom_point(aes(y = doy, shape = factor(ClusterID))) + 
  geom_point(aes(y = initAbs), color = 'red')  + 
  geom_point(aes(y = cessAbs), color = 'purple') 
  

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

