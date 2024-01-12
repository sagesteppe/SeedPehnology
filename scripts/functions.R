
#' collect all scanned images of a species from idigbio after a certain date
#' @param x Species name
#' @param date date of collection (e.g. '1979-01-01', the year many climate data sets start)
#' @param basepts number of base points
#' @param overpts number of oversample points
#' @param mindist minimum distance between points
#' @example aslo <- specimen_sampler('Astragalus lonchocarpus', '1979-01-01', 
#'      basepts = 70, overpts = 30, mindis = 2000)
#' @param ... further arguments passed onto spsurvey::grts 
specimen_sampler <- function(x, date, basepts, mindist, overpts, ...){
  
  text_result <- ridigbio::idig_search(rq = list(scientificname = x)) |> 
    dplyr::select(uuid, scientificname, datecollected, collector, geopoint.lon, geopoint.lat) |> 
    tidyr::drop_na(geopoint.lon, geopoint.lat, datecollected) |> 
    dplyr::mutate(datecollected = as.Date(datecollected)) |> 
    dplyr::filter(datecollected > date) |> # sensing data comes on line here
    dplyr::arrange(datecollected)
  
  if(nrow(text_result) > 2){
  
  media_result <- ridigbio::idig_search_media(
    rq = list(scientificname = x), mq = TRUE) |>
    dplyr::select(accessuri, records)
  
  media_result <- data.frame(apply(FUN = as.character, X = media_result, MARGIN = 2))
  
  records <- dplyr::inner_join(text_result, media_result, by = c('uuid' = 'records')) |>
    sf::st_as_sf(coords = c(x = 'geopoint.lon', y = 'geopoint.lat'), crs = 4326) |>
    dplyr::mutate(
      accessuri = gsub('[?]size=fullsize', '', accessuri), 
      collector = gsub("collector[(s):]", '', collector), 
      collector = gsub(';.*$|,.*$', '', collector),
      doy = as.numeric(lubridate::yday(datecollected))
    ) |>
    dplyr::select(scientificname, accessuri, datecollected, collector, doy, uuid)
  
    if(nrow(records) < 100){basepts <- (nrow(records) -1); overpts <- 1}
    
    samples <- spsurvey::grts(
      sf::st_transform(records, 5070), n_base = basepts, n_over = overpts, mindis = mindist, ...,
      seltype = "proportional", aux_var = 'doy')
    
    samples <- dplyr::bind_rows(
      samples[['sites_base']], samples[['sites_over']]
    ) |>
      dplyr::mutate(
        siteID = paste0(siteuse, gsub('[A-z]', '', siteID)), 
        collector = gsub(';.*$|,.*$', '', collector)
      ) |>
      dplyr::select(siteID, wgt, scientificname:geometry)
    
    return(samples)
  } else{
    'Taxon Not Found'
  }
  
}

