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


#' create pseudo-absences for phenophases 
#' @param x the output of geoClust 
#' @param estimates the output of initiation_cessation 
pheno_abs <- function(x, estimates){
  
  grps <- split(x, f = x$ClusterID)
  estimates <- split(estimates, f = estimates$ClusterID)
  
  ed <- unlist(lapply(estimates, doy_gen, 'initiation'))
  ld <- unlist(lapply(estimates, doy_gen, 'cessation'))
  
  # PCA1 1 = hottest, 0 = coldest. 
  # Give the earliest initiation DOY absence to the warmest place, 
  # and the latest absences to the coldest places. 
  x <- arrange(x, ClusterID, PC1)
  out <- data.frame(initAbs = ed, cessAbs = ld)
  out <- dplyr::bind_cols(x, data.frame(initAbs = ed, cessAbs = ld))
  
  # ensure that the cessation is > 1 weeks after observation and
  # initiation greater > 1 weeks before observation.
  out <- out |>
    dplyr::rowwise() |>
    dplyr::mutate(
      initAbs = if_else((doy - 10) < initAbs, (doy - 10), initAbs), 
      cessAbs = if_else((doy + 10) > cessAbs, (doy + 10), cessAbs)
    )
  
   # identify records in clusters which have a nearest neighbor belonging to a 
   # different cluster
   
   nearestClusterMismatch <- which( sf::st_drop_geometry(out[,'ClusterID']) 
          != sf::st_drop_geometry(out[dmat, 'ClusterID']) )
   
   tdat <- out[-nearestClusterMismatch,]
   # subset the raster surface to only cover the points which could benefit
   # from interpolated values. 
   
   tdat_hull <- vect(
     sf::st_convex_hull(
       sf::st_transform(
       sf::st_buffer(
         sf::st_transform(
           sf::st_union(
             out[nearestClusterMismatch,]), 
           5070), 
         750),
       crs(pca1))))
   
   surface <- terra::crop(pca1, tdat_hull, mask = TRUE)
   early_mod <- fields::Tps(st_coordinates(tdat), tdat$initAbs, Z = tdat$PC1)
   late_mod <- fields::Tps(st_coordinates(tdat), tdat$cessAbs, Z = tdat$PC1)
   
   early_surf <- terra::interpolate(surface, early_mod, fun=pfun)
   late_surf <- terra::interpolate(surface, late_mod, fun=pfun)
   names(early_surf) <- 'initAbs'
   names(late_surf) <- 'cessAbs'
   
   mismatch <- out[nearestClusterMismatch, !names(out) %in% c("initAbs", "cessAbs")]
   ef <- terra::extract(early_surf, mismatch, bind = TRUE)
   lf <- terra::extract(late_surf, ef, bind = TRUE) |>
     sf::st_as_sf() |>
     dplyr::mutate(across(initAbs:cessAbs, round))
   
   out <- dplyr::bind_rows(out, lf)
   return(out)
  # return(c(early_surf, late_surf))
   
}

ll <- pheno_abs(clusts, out)


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



