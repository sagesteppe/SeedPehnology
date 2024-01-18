library(sf)
library(tidyverse)
library(terra)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
achy <- read.csv('../data/processed/high_priority_sheets.csv') %>% 
  filter(scientificname == 'achnatherum hymenoides')

# subset to scored sheets, and add '0' when a phenophase was not observed. 
achy <- achy %>% 
  dplyr::select(-accessuri) |>
  filter(! if_all(Pct_Bud:Pct_Dropped, ~ is.na(.))) |>
  mutate(across(Pct_Bud:Pct_Dropped, ~ replace_na(.x, 0)))

# now add spatial attributes to data. These will be used to look up the 
# relevant independent variables
p <- '../data/processed/high_priority_sheets'
tab_data <- st_read(file.path(p, (paste0( basename(p), '.shp'))), quiet = T) %>% 
  filter(scntfcn == 'achnatherum hymenoides') %>% 
  dplyr::select(siteID, scientificname = scntfcn, doy)

achy <- left_join(achy, tab_data) |>
  st_as_sf() |>
  st_transform(4326)

rm(tab_data)
# we now have Day of year, and latitude. DOY will be a fixed predictor, and I anticipate
# will hold the most explanatory power, where other terms will be covariates for it



# now we can import variables which we believe correlate with flowering. 

p <- '../data/spatial/processed'
preds <- rast(file.path(p, list.files(p)))

achy <- terra::extract(preds, achy, bind = TRUE) |>
  st_as_sf()

## we don't have much scored data, so we will impute NA bulk density data ##
achy <- achy |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
         across(Pct_Bud:Pct_Dropped, ~ ifelse(.x > 0, 1, .x)))


## now we will include the longitude and latitude as a correlation term to avoid
# inflating the residuals


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
  
  inde <- as.matrix(sf::st_drop_geometry(x[,8:23]))
  de <- sf::st_drop_geometry(x$Pct_Anthesis)
  
  cl <- parallel::makeCluster(parallel::detectCores(), type='PSOCK')
  doParallel::registerDoParallel(cl)
  lmProfile <- caret::rfe(
    inde, # independent variables
    de, # dependent variables
    sizes  = seq.int(from = 2, to = 10, by = 2),
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
  
  urGamm <- MuMIn::uGamm
  formals(urGamm)$na.action <- 'na.omit' 
  formals(urGamm)$family <- 'binomial'
  formals(urGamm)$method <- 'REML'
  
  mod.aspatial <- urGamm(formula, data = x)
  mod.corExp <- urGamm(
    formula, data = x, correlation = nlme::corExp(form = cor_form, nugget=T))
  mod.corGaus <- urGamm(
    formula, data = x, correlation = nlme::corGaus(form = cor_form, nugget=T))
  mod.corSpher <- urGamm(
    formula, data = x, correlation = nlme::corSpher(form = cor_form, nugget=T))
  mod.corRatio <- urGamm(
    formula, data = x, correlation = nlme::corRatio(form = cor_form, nugget=T))
  mod.corLin <- urGamm(
    formula, data = x, correlation = nlme::corLin(form = cor_form, nugget=T))
  
  tf <- unlist(lapply(X = mget(ls(pattern = 'mod[.]')), FUN = is.null))
  rm(list = c(names(tf)[which(tf == TRUE)]))

  msel_tab <- MuMIn::model.sel(mget(ls(pattern = 'mod[.]')))
  top_mod <- row.names(msel_tab[1,])
  
  if(top_mod == 'mod.aspatial'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial')
  } else if(top_mod == 'mod.corExp'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial',
                            correlation = nlme::corExp(form = cor_form, nugget=T))
  } else if(top_mod == 'mod.corGaus'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial',
                            correlation = nlme::corGaus(form = cor_form, nugget=T))
  } else if(top_mod == 'mod.corSpher'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial',
                            correlation = nlme::corSpher(form = cor_form, nugget=T))
  } else if(top_mod == 'mod.corRatio'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial',
                            correlation = nlme::corRatio(form = cor_form, nugget=T))
  } else if(top_mod == 'mod.corLin'){
    mod.final <- mgcv::gamm(formula, data = x, method = 'ML', family = 'binomial',
                            correlation = nlme::corLin(form = cor_form, nugget=T))
  }
  
  # write out results to local #
  msel_tab <- data.frame(msel_tab) |>
    tibble::rownames_to_column('model')
  
  mod.final <- mod.final[['gam']]
  
  saveRDS(mod.final, 
          file.path('../results/models', paste0(gsub(' ', '_', taxon), '.rds')))
  
  colnames(msel_tab) <- gsub('^s[.]', '', colnames(msel_tab))
  colnames(msel_tab) <- gsub('\\..*$', '', colnames(msel_tab))
  msel_tab <- data.frame(lapply(msel_tab, function(y) if(is.numeric(y)) round(y, 5) else y) ) 
  write.csv(msel_tab, row.names = FALSE,
    file.path('../results/selection_tables', paste0(gsub(' ', '_', taxon), '.csv')))
  
  return(mod.final)
  
}

set.seed(28)
modeller(achy)


# play with prediction

pred_df <- expand.grid(doy = seq(min(achy$doy), max(achy$doy)), 
                       bulk_density = seq(min(achy$bulk_density), max(achy$bulk_density), by = 1))

pred_df <- expand.grid(doy = seq(50, 200), 
                       bulk_density = mean(achy$bulk_density)) # toy around to see results

pred <- predict(mod.aspatial, pred_df, type = 'response', se = TRUE)

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
