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
  
  fields2getrecord <- c(
    "data.dwc:scientificName",  "data.dwc:decimalLatitude", 
    "data.dwc:decimalLongitude", "collector", "uuid",
    "data.dwc:year", "data.dwc:month", "data.dwc:day")
  names(fields2getrecord) <- c(
    'scientificname', 'geopoint.lat', 'geopoint.lon', 'collector', 'uuid', 'year', 'month', 'day')
  
  text_result <- ridigbio::idig_search(rq = list(scientificname = x),
                                       fields = fields2getrecord) |> 
    dplyr::rename(all_of(fields2getrecord)) |>
    tidyr::drop_na(geopoint.lon, geopoint.lat, year, month, day) |> 
    tidyr::unite('datecollected', year, month, day, sep = '-') |>
    dplyr::mutate(datecollected = as.Date(datecollected), 
                  doy = lubridate::yday(datecollected)) |> 
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
    iterations = iter, percentile = 0.01, bootstraps = bs, conf = 0.9)
  cessation <- parallel::mclapply(
    X = grps, FUN = phenesse::weib_percentile_ci, mc.cores = n_cores,
    iterations = iter, percentile = 0.99, bootstraps = bs, conf = 0.9)
  
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


#' grab day of years before and after phenophase
#' @param x the output of initiation_cessation
#' @param event either 'initiation' or 'cessation'
doy_gen <- function(x, event){
  
  ci_event <- round(x[x$event == event, 'estimate'])
  
  if(event == 'cessation'){
    doys <- seq.int(ci_event + 28, ci_event)
  } else {
    doys <- seq.int(ci_event - 28, ci_event)
  }
  
  doys <- sample(doys, 
                 size = x[x$event == event, 'observations'],  replace = TRUE)
  
  return(sort(doys))
  
}


#' create an interpolated surface with a linear covariate
#' ?interpolate in terra for examples
pfun <- function(model, x, ...) {
  predict(model, x[,1:2], Z=x[,3], ...)
}

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
  
}


# assess the output from a function. 
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# only grab models which converged without any warnings. 
conv_ob <- function(x){
  if(is.null(x$warning) & is.null(x$warning) == TRUE){return(x[['value']])}
}

# match and run the top model. 
f_modeller <- function(model, type = c("mod.aspatial", "mod.corExp", "mod.corGaus",
                                   'mod.corSpher', 'mod.corRatio', 'mod.corLin'), data){
  type <- match.arg(type)
  switch(
    type,
    mod.aspatial = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial'),
    mod.corExp = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial',
                            correlation = nlme::corExp(form = cor_form, nugget=T)),
    mod.corGaus = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial',
                             correlation = nlme::corGaus(form = cor_form, nugget=)), 
    mod.corSpher = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial',
                              correlation = nlme::corSpher(form = cor_form, nugget=T)), 
    mod.corRatio = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial',
                              correlation = nlme::corRatio(form = cor_form, nugget=T)),
    mod.corLin = mgcv::gamm(formula, data = data, method = 'ML', family = 'binomial',
                            correlation = nlme::corLin(form = cor_form, nugget=T))
  )
}

#' fit some gams. 
modeller <- function(x){
  
  x <- x |>
    st_transform(5070) %>% 
    mutate(
      Latitude = sf::st_coordinates(.)[,2], 
      Longitude = sf::st_coordinates(.)[,1],
      .before = 'geometry') |>
    st_transform(4326)
  
  taxon <- sf::st_drop_geometry(x$scientificname)[1]
  
  # we will determine which features have utiltity in predicting whether a 
  # population in in a phenophase. 
  ctrl <- caret::rfeControl(functions = caret::gamFuncs,
                            method = "repeatedcv", number = 10,  repeats = 5,
                            verbose = FALSE, rerank = TRUE, 
                            allowParallel = TRUE, seeds = NA)
  
  inde <- as.matrix(sf::st_drop_geometry(x[,8:26]))
  de <- sf::st_drop_geometry(x$Pct_Anthesis)
  
  cl <- parallel::makeCluster(parallel::detectCores(), type='PSOCK')
  doParallel::registerDoParallel(cl)
  lmProfile <- caret::rfe(
    as.matrix(sf::st_drop_geometry(x[,8:26])), # independent variables
    sf::st_drop_geometry(x$Pct_Anthesis), # dependent variables
    sizes  = c(3, seq.int(from = 4, to = 10, by = 2)),
    rank = TRUE, 
    rfeControl = ctrl)
  ParallelLogger::stopCluster(cl)
  rm(cl)
  
  terms <- unique(c('doy', lmProfile[['optVariables']]))
  formula <- as.formula(
    paste(
      "Pct_Anthesis ~ ",  paste0("s(", terms, ", bs = 'tp')", collapse = " + ")
    )
  )
  cor_form <- as.formula("~ Latitude + Longitude")
  
  urGamm <- MuMIn::uGamm ; formals(urGamm)$na.action <- 'na.omit' 
  formals(urGamm)$family <- 'binomial'; formals(urGamm)$method <- 'REML'
  
  mod.aspatial <- conv_ob(myTryCatch(urGamm(formula, data = x)))
  mod.corExp <- conv_ob(myTryCatch(urGamm(
    formula, data = x, correlation = nlme::corExp(form = cor_form, nugget=T))))
  mod.corGaus <- conv_ob(myTryCatch(urGamm(
    formula, data = x, correlation = nlme::corGaus(form = cor_form, nugget=T))))
  mod.corSpher <- conv_ob(myTryCatch(urGamm(
    formula, data = x, correlation = nlme::corSpher(form = cor_form, nugget=T))))
  mod.corRatio <- conv_ob(myTryCatch(urGamm(
    formula, data = x, correlation = nlme::corRatio(form = cor_form, nugget=T))))
  mod.corLin <- conv_ob(myTryCatch(urGamm(
    formula, data = x, correlation = nlme::corLin(form = cor_form, nugget=T))))
  
  tf <- unlist(lapply(X = mget(ls(pattern = 'mod[.]')), FUN = is.null))
  rm(list = c(names(tf)[which(tf == TRUE)]))
  
  # select and refit the top model with method 'ML'
  msel_tab <- MuMIn::model.sel(mget(ls(pattern = 'mod[.]')))
  top_mod <- row.names(msel_tab[1,])
  mod.final <- f_modeller(model = top_mod, data = x)
  
  # write out results to local #
  msel_tab <- data.frame(msel_tab) |>
    tibble::rownames_to_column('model')
  
  saveRDS(mod.final, 
          file.path('../results/models', paste0(gsub(' ', '_', taxon), '.rds')))
  
  colnames(msel_tab) <- gsub('^s[.]', '', colnames(msel_tab))
  colnames(msel_tab) <- gsub('\\..*$', '', colnames(msel_tab))
  msel_tab <- data.frame(
    lapply(msel_tab, function(y) if(is.numeric(y)) round(y, 5) else y)) 
  write.csv(msel_tab, row.names = FALSE,
            file.path('../results/selection_tables', paste0(gsub(' ', '_', taxon), '.csv')))
  
  return(mod.final)
  
}


