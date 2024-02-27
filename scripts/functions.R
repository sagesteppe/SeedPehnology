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
    dplyr::filter(year > 1900) |>
    tidyr::unite('datecollected', year, month, day, sep = '-') |>
    dplyr::mutate(datecollected = as.Date(datecollected), 
                  doy = lubridate::yday(datecollected)) |> 
    dplyr::filter(datecollected > date, 
                  geopoint.lat > 1) |> # sensing data comes on line here
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
  
  
 #   if(nrow(records) < 100){basepts <- (nrow(records) -1); overpts <- 1}
    
 #   samples <- spsurvey::grts(
#      sf::st_transform(records, 5070), n_base = basepts, n_over = overpts, 
#      mindis = mindist, ...,
#      seltype = "proportional", aux_var = 'doy')
    
#    samples <- dplyr::bind_rows(
#      samples[['sites_base']], samples[['sites_over']]
#    ) |>
#      dplyr::mutate(
#        siteID = paste0(siteuse, gsub('[A-z]', '', siteID)), 
#        collector = gsub(';.*$|,.*$', '', collector)
#      ) |>
#      dplyr::select(siteID, wgt, scientificname:geometry)
    
#    return(samples)
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

  # if any cluster < 10 members, add them to the most geographically proximal group. 
  ob <- dplyr::bind_cols(x, ClusterID = clusters) 
  
  small_clust <- ob |>
    dplyr::group_by(ClusterID) |>
    dplyr::filter(n() <= 10) |>
    dplyr::pull(ClusterID)
  
  if(length(small_clust) > 0){
    obs_small <- dplyr::filter(ob, ClusterID %in% small_clust)
    ref <- dplyr::filter(ob, !ClusterID %in% small_clust)
    
    if(nrow(ref) == 0){
      ob$ClusterID <- 1
    } else {
    newClusts <- sf::st_drop_geometry(
      ref[ sf::st_nearest_feature(obs_small, ref), 'ClusterID'])
    ob <- dplyr::bind_rows(
      dplyr::bind_cols(select(obs_small, -ClusterID), newClusts) , ref) |>
      dplyr::mutate(ClusterID = as.numeric(ClusterID))
    rownames(ob) <- 1:nrow(ob)
    }
  }
  return(ob)   # return the modified input object
}


#' calculate the initiation and cessation of a phenophase using phenesse 
#' this function is intended to calculate these variables for clusters of populations
#' of a species. For example the output of geoClust_wrap
#' @param x output of geoClust_wrap 
#' @param iter the number of iterations to perform see phenesse::weib_percentile_ci. defaults to 250 here, 500 is their suggestion. 
#' @param bs the number of bootstraps to perform to generate a confidence interval, note these BS take  considerably longer than a typical bs sampling event, Sys.time(initiation_cessation, bs = 10) to get a feel. Defaults to the number of supplied cores.  
#' @param n_cores the number of cores to use while processing the records, defaults to all cores. 
initiation_cessation <- function(x, iter, bs, n_cores){
  
  if(missing(iter)){iter <- 250; 
  message('Number of iterations not supplied to `iter` defaults to 250')}
  if(missing(n_cores)){n_cores <- parallel::detectCores(); 
  message('Number of cores not supplied to `n_cores` defaults to all, e.g. `parallel::detectCores()`')}
  if(missing(bs)){bs <- n_cores;
  message('Number of bootstrap reps not supplied to `bs` defaults to n_cores')}
  
  if(length(unique(x$ClusterID)) > 1){
    
    grps <- split(x, f = x$ClusterID)
    grps <- lapply(grps, '[[', 'doy')
    
    initiation <- parallel::mclapply(
      X = grps, FUN = phenesse::weib_percentile_ci, mc.cores = n_cores,
      iterations = iter, percentile = 0.01, bootstraps = bs, conf = 0.9)
  
    cessation <- parallel::mclapply(
      X = grps, FUN = phenesse::weib_percentile_ci, mc.cores = n_cores,
      iterations = iter, percentile = 0.99, bootstraps = bs, conf = 0.9)
  } else {
    vals <- x$doy
    initiation <- 
      phenesse::weib_percentile_ci(vals, iterations = iter, percentile = 0.01, 
                                   bootstraps = bs, conf = 0.9)
    cessation <- 
      phenesse::weib_percentile_ci(vals, iterations = iter, percentile = 0.99, 
                                   bootstraps = bs, conf = 0.9)
  }
  
  results <- dplyr::bind_rows(
    data.frame(event = 'initiation', dplyr::bind_rows(initiation, .id = "ClusterID")), 
    data.frame(event = 'cessation', dplyr::bind_rows(cessation, .id = "ClusterID"))
  )
  
  if(length(unique(x$ClusterID)) > 1){
    no_obs <- data.frame(
      ClusterID = names(grps),
      observations = sapply(X = grps, FUN = length) 
    )} else {
    no_obs <- data.frame(
      ClusterID = 1, 
      observations = length(vals)
    )
  }
  
  results <- merge(results, no_obs, by = 'ClusterID')
  
  return(results)
}

#' grab day of years before and after phenophase
#' @param x the output of initiation_cessation
#' @param event either 'initiation' or 'cessation'
doy_gen <- function(x, eventA){

  if(eventA == 'cessation'){
    ci_event <- as.numeric(round(x[x$event == eventA, 'high_ci']))
    doys <- seq.int(ci_event + 28, ci_event+7)
  } else {
    eventA <- 'initiation'
    ci_event <- as.numeric(round(x[x$event == eventA, 'low_ci']))
    doys <- seq.int(ci_event - 28, ci_event-7)
  } 
  
  doys <- sample(doys, 
                 size = as.numeric(x[x$event == eventA, 'observations']),  replace = TRUE)
  
  if(eventA == 'cessation'){
    return(sort(doys, decreasing = TRUE))} else {
      return(sort(doys, decreasing = TRUE))}
  
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
  out1 <- data.frame(initAbs = ed, cessAbs = ld)
  if(nrow(out1) > nrow(x)){
    out1 <- out1[sample(nrow(out1), size = nrow(x)),]
  } else if (nrow(out1) < nrow(x)){
    message( (nrow(x) - nrow(out1)))
    
    out1 <- rbind(out1,
                  out1[ sample(nrow(out1), size = (nrow(x) - nrow(out1))),])
  }
  out <- dplyr::bind_cols(x, out1)
  
  # ensure that the cessation is > 3 days after observation and
  # initiation greater > 3 days before observation.
  out <- out |>
    dplyr::rowwise() |>
    dplyr::mutate(
      initAbs = if_else((doy - 7) < initAbs, (doy - 7), initAbs), 
      cessAbs = if_else((doy + 7) > cessAbs, (doy + 7), cessAbs)
    )
  
  # identify records in clusters which have a nearest neighbor belonging to a 
  # different cluster
  dmat <- sf::st_nearest_feature(out)
  nearestClusterMismatch <- which( sf::st_drop_geometry(out[,'ClusterID']) 
                                   != sf::st_drop_geometry(out[dmat, 'ClusterID']) )
  
  tdat <- out[-nearestClusterMismatch,]
  # subset the raster surface to only cover the points which could benefit
  # from interpolated values. 
  
  if(length(nearestClusterMismatch) > 0){
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
    names(early_surf) <- 'initAbs'; names(late_surf) <- 'cessAbs'
  
    mismatch <- out[nearestClusterMismatch, !names(out) %in% c("initAbs", "cessAbs")]
    ef <- terra::extract(early_surf, mismatch, bind = TRUE)
    lf <- terra::extract(late_surf, ef, bind = TRUE) |>
      sf::st_as_sf() |>
      dplyr::mutate(across(initAbs:cessAbs, round))
  
    out <- dplyr::bind_rows(out, lf)
  } 
  return(out)
}

#' generate phenology absences

pheno_abs_writer <- function(x){
  
  nf <- x[ st_nearest_feature(x), ]
  recs2clust <- x[ as.numeric( st_distance(x, nf, by_element = T) ) < 80000 , ] 
  clusts <- geoClust_wrap(recs2clust, 'PC1') # identify clusters here.
  
  r_taxon <- gsub(' ', '_', sf::st_drop_geometry(x$scntfcnm[1]))
  est <- dplyr::filter(tables, scntfcnm == r_taxon)
  out <- pheno_abs(clusts, est) |>
    dplyr::select(scntfcnm, accessr, doy, initAbs, cessAbs)
  
  sf::st_write(out, paste0('../results/PresAbs/', r_taxon, '.shp'), quiet = TRUE, append = FALSE)
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
  if(is.null(x$error) & is.null(x$warning) == TRUE){return(x[['value']])}
}

f_modeller <- function(model, m_form, type = c("mod.aspatial", "mod.corExp", "mod.corGaus",
                                   'mod.corSpher', 'mod.corRatio', 'mod.corLin'), data){
  
  type <- match.arg(type)
  switch(
    type,
    mod.aspatial = mgcv::gam(m_form, data, method = 'ML',
                             family = 'binomial', select = TRUE),
    mod.corExp = mgcv::gamm(m_form, data, method = 'ML', family = 'binomial',
                            correlation = nlme::corExp(form = cor_form, nugget=T)),
    mod.corGaus = mgcv::gamm(m_form, data, method = 'ML', family = 'binomial',
                             correlation = nlme::corGaus(form = cor_form, nugget=T)), 
    mod.corSpher = mgcv::gamm(m_form, data, method = 'ML', family = 'binomial',
                              correlation = nlme::corSpher(form = cor_form, nugget=T)), 
    mod.corRatio = mgcv::gamm(m_form, data, method = 'ML', family = 'binomial',
                              correlation = nlme::corRatio(form = cor_form, nugget=T)),
    mod.corLin = mgcv::gamm(m_form, data, method = 'ML', family = 'binomial',
                            correlation = nlme::corLin(form = cor_form, nugget=T))
  )
}



#' fit some gams. 
modeller <- function(x){
  
  x <- st_transform(x, 5070)
  x <- dplyr::mutate(x, 
        Latitude = jitter(sf::st_coordinates(x)[,2])/1000, 
        Longitude = jitter(sf::st_coordinates(x)[,1])/1000,
        .before = 'geometry') |>
    st_transform(4326)
  
  taxon <- sf::st_drop_geometry(x$scntfcnm)[1]
  
  # we will determine which features have utility in predicting whether a 
  # population has individuals in a phenophase. 
   ctrl <- caret::rfeControl(functions = caret::gamFuncs,
                             method = "repeatedcv", number = 10,  repeats = 5,
                             verbose = FALSE, rerank = TRUE, 
                             allowParallel = TRUE, seeds = NA)
  
  col_range <- (grep('flowering', colnames(x)) + 1):(grep('geometry',  colnames(x)) - 3)
  inde <- as.matrix(sf::st_drop_geometry(x[,col_range]))
  de <- sf::st_drop_geometry(x$flowering)
  
  cl <- parallel::makeCluster(parallel::detectCores(), type='PSOCK')
  doParallel::registerDoParallel(cl)
  
  ob <- myTryCatch(
     caret::rfe(
     as.matrix(sf::st_drop_geometry(x[,col_range])), # independent variables
     sf::st_drop_geometry(x$flowering), # dependent variables
     sizes  = c(2:5),
     rank = TRUE, 
     rfeControl = ctrl)
  )
  
  if(is.null(ob$error) == TRUE){lmProfile <- ob[['value']]} else{
    'unable to perform recursive feature elimination with this species. sample size
    is likely too small to accomodate the number of features and knots, peforming
    selection with generalized linear models.'
    ctrl <- caret::rfeControl(functions = caret::lmFuncs,
                              method = "repeatedcv", number = 10,  repeats = 5,
                              verbose = FALSE, rerank = TRUE, 
                              allowParallel = TRUE, seeds = NA)
    lmProfile <-  caret::rfe(
       as.matrix(sf::st_drop_geometry(x[,col_range])), # independent variables
       sf::st_drop_geometry(x$flowering), # dependent variables
       sizes  = c(2:5),
       rank = TRUE, 
       rfeControl = ctrl)
  }
  invisible(ParallelLogger::stopCluster(cl))
  rm(cl)
  
  terms <- lmProfile[['optVariables']][grep('doy', lmProfile[['optVariables']], invert = TRUE)]
  if( sum(grepl('gddlgd', terms)) > 1) {
    pos <- grep('gddlgd', terms)
    pos <- pos[2:length(pos)]
    terms <- terms[-pos]
  }
  if(length(terms) > 3){terms <- terms[1:3]}
  if(nrow(x) < 50){kval <- 15} else {kval <- 25}
  m_form <- as.formula(
      paste(
        "flowering ~ ",  "s(doy, bs = 'cc', k = ", kval, ") + ", # smooth for cyclic data
        paste0(terms, collapse = " + ")
      )
  )

  message('feature selection complete, fitting model ', m_form)
  mm <- mods(x = x, m_form, round = 1)
  msel_tab <- mm[[1]] ;  mod.final <- mm[[2]] 

 #  determine whether all terms were needed. - if not remove them. 
  
  new_form <- new_form_fn(mod.final, terms, m_form)
  if(m_form != new_form){
    message("Refitting a final model ", new_form, " without a covariate(s) for better prediction")
    top_mod <- row.names(msel_tab[1,])
    mm <- mods(x, m_form = new_form, round = 2, model = top_mod); 
    msel_tab <- mm[[1]] ;  mod.final <- mm[[2]] 
    }
  
  # write out results to local 
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
  
}

#' fit gamms with different spatial terms. 
mods <- function(x, m_form, round, model){
  
  cor_form <- as.formula("~ Latitude + Longitude")
  urGamm <- MuMIn::uGamm ; formals(urGamm)$na.action <- 'na.omit' 
  formals(urGamm)$family <- 'binomial'; formals(urGamm)$method <- 'REML'
  
  if(round == 1){
  mod.aspatial <- mgcv::gam(m_form, data = x, family = 'binomial')
  if(morans_wrapper(x, mod.aspatial) == TRUE){
    rm(mod.aspatial)
    message('Fitting spatial-autocorrelation error structures.')
   mod.corExp <- conv_ob(myTryCatch(urGamm(
      m_form, data = x, correlation = nlme::corExp(form = cor_form, nugget=T), family = 'binomial')))
    mod.corGaus <- conv_ob(myTryCatch(urGamm(
      m_form, data = x, correlation = nlme::corGaus(form = cor_form, nugget=T), family = 'binomial')))
    mod.corSpher <- conv_ob(myTryCatch(urGamm(
      m_form, data = x, correlation = nlme::corSpher(form = cor_form, nugget=T), family = 'binomial')))
    mod.corRatio <- conv_ob(myTryCatch(urGamm(
      m_form, data = x, correlation = nlme::corRatio(form = cor_form, nugget=T), family = 'binomial')))
    mod.corLin <- conv_ob(myTryCatch(urGamm(
      m_form, data = x, correlation = nlme::corLin(form = cor_form, nugget=T), family = 'binomial')))
  }
  
  tf <- unlist(lapply(X = mget(ls(pattern = 'mod[.]')), FUN = is.null))
  rm(list = c(names(tf)[which(tf == TRUE)]))

  # select and refit the top model with method 'ML'
  msel_tab <- MuMIn::model.sel(mget(ls(pattern = 'mod[.]')))
  top_mod <- row.names(msel_tab[1,])
  mod.final <- f_modeller(top_mod, m_form = m_form, data = x)
  } else {mod.final <- f_modeller(model, m_form = m_form, data = x)
  msel_tab <- MuMIn::model.sel(mget(ls(pattern = 'mod[.]')))
  rownames(msel_tab) <- model}
  
  return(list(msel_tab, mod.final))
  
}

#' this function flags records which warrant visual review as they are at extreme 
#' ends of the distributions
visual_review_flagger <- function(x, probs){
  
  ob_r <- max(x$PC1, na.rm = T) - min(x$PC1, na.rm = T) # range of PC1 values. 
  
  flggr <- function(x){
    qu <- quantile(x$doy,  probs = probs)
    x$phen_flag <- ifelse(x$doy < qu[1] | x$doy > qu[2], 'Broad', NA)
    dplyr::relocate(x, phen_flag, .before = geometry)
    return(x)
  }
  
  if(ob_r < 0.3){
    
    x <- flggr(x)
    
  } else {
    
    x <- x |>
      dplyr::mutate(
        cuts = cut(PC1, breaks = 3, labels = F), 
        qu_b_l = quantile(doy, probs = probs[1]), 
        qu_b_h = quantile(doy, probs = probs[2])
      ) |> 
      dplyr::group_by(cuts) |> 
      dplyr::mutate(
        qu_f_l = quantile(doy, probs = probs[1]), 
        qu_f_h = quantile(doy, probs = probs[2]), 
      ) |> 
      dplyr::rowwise() |> 
      dplyr::mutate(phen_flag = case_when(
        doy < qu_b_l | doy > qu_b_h ~ 'Broad',
        doy < qu_f_l | doy > qu_f_h ~ 'Fine'
      ), .before = geometry) |>
      dplyr::select(-cuts, -qu_b_l, -qu_b_h, -qu_f_l, -qu_f_h)
    
  }
  return(x)
}

#' quick and dirty drop far off records and close dupes - see other repos for proper drops. 
rec_selec <- function(x){
  
  x <- st_transform(x, 5070)
  nf <-  sf::st_nearest_feature(x)
  y <- x[ as.numeric(sf::st_distance(x, x[nf,], by_element = T) ) < 80000 , ]
  # remove far off records.
  
  # remove records right on top of each other. 
  nf <- sf::st_nearest_feature(y)
  z <- y[ as.numeric(sf::st_distance(y, y[nf,], by_element = T)) < 250, ] %>% 
    sf::st_buffer(250) 
  y <- y[ as.numeric(sf::st_distance(y, y[nf,], by_element = T)) > 250, ]
  
  z$grps <- sf::st_intersects(z) 
  z <- z %>% 
    dplyr::group_by(grps) %>% 
    dplyr::slice_sample(n = 1) %>% 
    sf::st_centroid() %>% 
    dplyr::select(-grps)
  
  out <- dplyr::bind_rows(y, z) %>% 
    sf::st_transform(4326) 
  
  return(out)
}



# identify local minima
l_min <- function(x){
  
  localMinima <- function(x) { # @ tommy on SO. 
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
      y <- y[-1]
    }
    y
  }
  
  ident_min <- x$doy[localMinima(x$valu) ] 
  v_low <- x[x$doy %in% ident_min, ] |>
    dplyr::group_by(doy) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::pull(valu)
  
  return(list(v_low, ident_min))
  
}



# identify local maxima
l_max <- function(x){
  
  localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
      y <- y[-1]
    }
    y
  }
  
  ident_max <- x$doy[localMaxima(x$valu) ] 
  v_high <- x[x$doy %in% ident_max, ] |>
    dplyr::group_by(doy) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::pull(valu)
  
  return(list(v_high, ident_max))

}


#' find splits for modal distributions at a set threshold and plot all data for visual review

modal_finder <- function(x, path){
  
  dd <- density(x$doy)
  
  vals <- data.frame(
    valu = dd$y[dd$x >= 0 & dd$x <= 366], 
    doy = round(dd$x[dd$x >= 0 & dd$x <= 366], 0)
  )
  
  low <- l_min(vals); high <- l_max(vals)
  v_low <- low[[1]]; v_high <- high[[1]] 
  ident_min <- low[[2]]; ident_max <- high[[2]]
  
  out <- data.frame( # no longer returned. 
    'DOY' = c(ident_min, ident_max),
    'Value' = c(v_low, v_high), 
    'Variable' =  c(rep('Min', length(ident_min)), rep('Max', length(ident_max))), 
    'prcnt.max' = (c(v_low, v_high) / max(v_high)) * 100
  )  
  
  if(nrow(out[out$Variable == 'Max' & out$prcnt.max >= 10.0, ]) > 1){
    
    doys <- out[out$Variable == 'Max' & out$prcnt.max >= 10.0, 'DOY']
    min <- out[out$Variable == 'Min' & out$DOY > doys[1] & out$DOY < doys[2],]
    min <- min[ which.min(min$Value), ]
    
    if(min$prcnt.max >= 70){min <- 999} else{min <- min$DOY}
  } else {min <- 999}
  
  
  png(filename = paste0(path, 
                        gsub(' ', '_', sf::st_drop_geometry(x$scntfcnm[1])), '.png'),
      width = 480, height = 320)
  hist(x = x$doy, main = x$scntfcnm[1], xlab = 'Day of Year (DOY)', prob = TRUE, 
       ylim = c(0, max(vals$valu)+0.005), col = '#F3E9DC', xlim = c(0, 365))
  axis(side = 1, at = seq(0, 350, by = 50), labels = T)
  lines(y = vals$valu, x = vals$doy, lwd = 2)
  abline(v = ident_min , col="#007FFF", lty = 2) 
  abline(v = ident_max, col = '#590925', lty = 2)
  abline(v = min, col = 'red', lwd = 3)
  legend("topright", legend=c("Maxima", "Minima", 'Split'),
         col=c("#590925", "#007FFF", 'red'), lty = c(2, 2, 1), cex=0.8)
  mtext(paste0('n = ', length(x$doy)), side=1, line=3.5, at=350)
  dev.off()
  
  return(min)
  
}


#' subset SPEI data to domain and resample to grains
speidR <- function(x, template, pout){
  
  x <- rast(x)
  
  dttm <- c(paste0(rep(1950:2022, each = 12), '-', month.abb), 
            paste0(rep(2023, each = 9), '-', month.abb[1:9]))
  lyr_name <- gsub('[.]nc', '', basename(sources(x)))
  L <- dim(x)[3]
  
  names(x)  <- paste0(lyr_name, '-', dttm)
  x <- x[[361:L]]
  
  x <- terra::crop(template, x)
  x <- terra::resample(template, x)
  
  terra::writeCDF(x, 
                  filename = file.path(pout, paste0(lyr_name, '.nc')), 
                  varname = lyr_name) 
}



ince_writer <- function(x, bs){
  
  nf <- x[ st_nearest_feature(x), ]
  recs2clust <- x[ as.numeric( st_distance(x, nf, by_element = T) ) < 80000 , ] 
  clusts <- geoClust_wrap(recs2clust, 'PC1') # identify clusters here.
  
  ince <- initiation_cessation(clusts, n_cores = no_cores, iter = 250)
  taxon <- gsub(' ', '_', sf::st_drop_geometry(x$scntfcnm[1]))
  write.csv(ince, paste0(
    '../results/initation_cessation_tables/', taxon, '.csv'), row.names = F)
  
  message(taxon, ' complete ' , Sys.time(), ' cooling down for 60 seconds before starting next species')
  Sys.sleep(60) # let the computer cool down for 60 seconds before re-engaging.
}


#' determine whether a co-variate should be dropped from a final model 
new_form_fn <- function(x, terms, m_form){
  
  ob <- broom::tidy(x, parametric = TRUE) |>
    dplyr::filter(term != '(Intercept)')
  
  if(nrow(ob) > 1){
    
    terms2remove <- dplyr::filter(ob, p.value > 0.1) |>  
      dplyr::pull(term) # identify terms to remove.
    
    if(length(terms2remove)>0){
      terms_sub <- terms[- grep(paste(terms2remove, collapse = "|"), terms) ] 
      if(length(terms_sub) > 0){
        m_form <- as.formula(
          paste(
            "flowering ~ ",  "s(doy, bs = 'cc', k = ", kval, ") + ", # smooth for cyclic data
            paste0(terms_sub, collapse = " + ")
          )
        )
      } else {
        
        terms2remove <- dplyr::arrange(ob, p.value) |>
          dplyr::pull(term)
        terms2remove <- terms2remove[2:length(terms2remove)]
        terms_sub <- terms[- grep(paste(terms2remove, collapse = "|"), terms) ] 
        m_form <- as.formula(
          paste(
            "flowering ~ ",  "s(doy, bs = 'cc', k = ", kval, ") + ", # smooth for cyclic data
            paste0(terms_sub, collapse = " + ")
          )
        )
      }
      return(m_form)
    } else {return(m_form)}
  } else {return(m_form)}
}


#' determine whether we should model spatial autocorrelation
morans_wrapper <- function(x, model){
  
  mat <- x
  mat[,c('Grp', 'x', 'y')] <- rep(1:(nrow(x)/3), each = 3)
  simulationOutput <- DHARMa::simulateResiduals(model)
  simulationOutput1 <- DHARMa::recalculateResiduals(simulationOutput, group = mat$Grp)
  matty <- dplyr::distinct(mat, Grp, .keep_all = T) |>
    sf::st_drop_geometry()
  
  morans <- DHARMa::testSpatialAutocorrelation(
    simulationOutput1, matty$x, matty$y, plot = FALSE)
  return(morans[["p.value"]] < 0.01) # if we reject the null hypothesis, that no spatial pattern exists, 
  # than we will create the spatial models. 
}

# reload the rasters into a stack - ensure they are in ascending DOY Order. 
orderLoad <- function(path){
  doy_f <- list.files(path)
  doy_stack <- paste0(path, 
                      sort(as.numeric(gsub('.tif', '', 
                                           basename(doy_f)))), '.tif')
  doy_stack <- terra::rast(lapply(doy_stack, terra::rast))
  return(doy_stack)
}


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
