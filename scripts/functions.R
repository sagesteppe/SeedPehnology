
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
      sf::st_transform(records, 5070), n_base = basepts, n_over = overpts, 
      mindis = mindist, ...,
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


#' this function wraps around Clustgeo using sensible automatic selections for good results
#' @param x an sf/datamframe/tibble
#' @param clust_vars a vector of column names to be used for clustering
geoClust_wrap <- function(x, clust_vars){
  
  dat <- sf::st_drop_geometry(x) |>
    dplyr::select(any_of(clust_vars)) 
  D0 <- dist(dat)
  
  # first perform a clustering based on attributes of the variables 
  tree <- ClustGeo::hclustgeo(D0)
  op_K <- maptree::kgs(tree, D0, maxclus = 10) # an at least decent K. 
  K <- as.numeric(names(op_K[which(op_K == min(op_K))]))
  
  # now incorporate space to the existing clusters to realign them
  D1 <- as.dist(sf::st_distance(x))
  
  # the alpha parameter blends the contributions of the data attributes and
  # geographic distance
  cr <- ClustGeo::choicealpha(D0, D1, seq(0.0, 1, 0.1), 
                              K, graph = FALSE)
  
  tab <- data.frame(cr$Q)
  tab <- tab[tab$Q0 > tab$Q1,]
  opt_alpha <- as.numeric(
    gsub(
      'alpha=', '', rownames(tab[ which.max(tab$Q1), ])
    )
  )
  
  # cluster the trees
  tree <- ClustGeo::hclustgeo(D0, D1, alpha = opt_alpha)
  clusters <- cutree(tree, K)
  
  # return the modified input object
  ob <- dplyr::bind_cols(x, ClusterID = clusters) 
  return(ob)
}


#' calculate the initiation and cessation of a phenophase using phenesse 
#' this function is intended to calculate these variables for clusters of populations
#' of a species. For example the output of geoClust_wrap
#' @param x output of geoClust_wrap
#' @param iter the number of iterations to perform see phenesse::weib_percentile_ci. defaults to 250 here, 500 is their suggestion.
#' @param bs the number of bootstraps to perform to generate a confidence interval, note these BS take considerably longer than a typical bs sampling event, Sys.time(initiation_cessation, bs = 10) to get a feel. Defaults to the number of supplied cores. 
#' @param n_cores the number of cores to use while processing the records, defaults to all cores.
initiation_cessation <- function(x, iter, bs, n_cores){
  
  if(missing(iter)){iter <- 250; 
  message('# of iterations not supplied to `iter` defaults to 250')}
  if(missing(n_cores)){n_cores <- parallel::detectCores(); 
  message('# of cores not supplied to `n_cores` defaults to all, e.g. parallel::detectCores()')}
  if(missing(bs)){bs <- n_cores;
  message('# of bootstrap reps not supplied to `bs` defaults to n_cores')}
  
  grps <- split(x, f = x$ClusterID)
  grps <- lapply(grps, '[[', 'doy')
  
  initiation <- parallel::mclapply(
    X = grps, FUN = phenesse::weib_percentile_ci, mc.cores = n_cores,
    iterations = iter, percentile = 0.01, bootstraps = bs)
  cessation <- parallel::mclapply(
    X = grps, FUN = phenesse::weib_percentile_ci, mc.cores = n_cores,
    iterations = iter, percentile = 0.99, bootstraps = bs)
  
  results <- dplyr::bind_rows(
    data.frame(event = 'initiation', dplyr::bind_rows(initiation, .id = "ClusterID")), 
    data.frame(event = 'cessation', dplyr::bind_rows(cessation, .id = "ClusterID"))
  )
  no_obs <- data.frame(
    ClusterID = names(grps),
    observations = sapply(X = grps, FUN = length) 
  )
  results <- merge(results, no_obs, by = 'ClusterID')
  
  return(results)
}

