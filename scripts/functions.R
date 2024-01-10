
#' collect all scanned images of a species from idigbio after a certain date
#' @param x Species name
#' @param date date of collection (e.g. '1979-01-01', the year many climate data sets start)
#' @example aslo <- specimen_sampler('Astragalus lonchocarpus', '1979-01-01', 
#'      n_base = 70, n_over = 30, mindis = 2000)
#' @param ... further arguments passed onto spsurvey::grts 
specimen_sampler <- function(x, date, ...){
  
  text_result <- ridigbio::idig_search(rq = list(scientificname = x)) |> 
    dplyr::select(uuid, scientificname, datecollected, collector, geopoint.lon, geopoint.lat) |> 
    tidyr::drop_na(geopoint.lon, geopoint.lat) |> 
    dplyr::mutate(datecollected = as.Date(datecollected)) |> 
    dplyr::filter(datecollected > date) |> # sensing data comes on line here
    dplyr::arrange(datecollected)
  
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
  
  samples <- spsurvey::grts(
      sf::st_transform(records, 5070), ...,
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
  
}
