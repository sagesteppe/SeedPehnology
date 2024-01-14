library(sf)
library(tidyverse)
library(mgcv)
library(terra)
library(spThin)

# thin(loc.data = ,  lat.col = , long.col = , reps = 100, thin.par = 3, )
setwd('~/Documents/SeedPhenology/scripts')
achy <- read.csv('../data/processed/high_priority_sheets.csv') %>% 
  filter(scientificname == 'achnatherum hymenoides')

# subset to scored sheets, and add '0' when a phenophasewas not observed. 
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



# the gams require absences of flowering both before and after the observed records. 
# we will generate an equal amount of non-flowering records to the flowering data set. 




# now we can import variables which we believe correlate with flowering. 

p <- '../data/spatial/processed'
preds <- rast(file.path(p, list.files(p)))

achy <- extract(preds, achy, bind = TRUE) |>
  st_as_sf()

## we don't have much scored data, so we will impute NA bulk density data ##
achy <- achy |>
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
         across(Pct_Bud:Pct_Dropped, ~ ifelse(.x > 0, 1, .x)))


## now we will include the longitude and latitude as a correlation term to avoid
# inflating the residuals

achy <- achy |>
  st_transform(5070) %>% 
  mutate(
    Latitude = sf::st_coordinates(.)[,2], 
    Longitude = sf::st_coordinates(.)[,1],
     .before = 'geometry') |>
  st_transform(4326)

  
# we will determine which features have utiltity in predicting whether a 
# population in in a phenophase. 
ctrl <- caret::rfeControl(functions = gamFuncs,
                   method = "repeatedcv", number = 10,  repeats = 5,
                 s  verbose = FALSE, rerank = TRUE, 
                   allowParallel = TRUE, seeds = NA)

inde <- as.matrix(sf::st_drop_geometry(achy[,8:23]))
de <- sf::st_drop_geometry(achy$Pct_Anthesis)

cl <- parallel::makeCluster(parallel::detectCores(), type='PSOCK')
doParallel::registerDoParallel(cl)
lmProfile <- caret::rfe(
  inde, # independent variables
  de, # dependent variables
  sizes  = c(1:10),
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

g_model <- mgcv::gam(formula,
            na.action = 'na.omit', family = 'binomial', 
            correlation = corExp(form = ~ Latitude + Longitude),
            data = achy, select = TRUE, method = 'ML')

# now we add a spatial autocorrelation term 

# now refit the top model using method = 'ML' 

par(mar=c(4,4,2,2))
plot(g_model)
summary(g_model)
gam.check(g_model)
concurvity(g_model, full=TRUE) 
anova(g_model) 

g_model <- update(g_model, method='ML', na.action='na.fail')
summary(g_model)

rm(terms, formula)








# play with prediction

pred_df <- expand.grid(doy = seq(min(achy$doy), max(achy$doy)), 
                       gdd10 = seq(min(achy$gdd10), max(achy$gdd10), by = 1))

pred_df <- expand.grid(doy = seq(50, 200), 
                       gdd10 = mean(achy$gdd10)) # toy around to see results

pred <- predict(g_model, pred_df, type = 'response', se = TRUE)


plot(pred_df$doy, pred$fit)

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
