library(MuMIn)
library(sf)
library(tidyverse)
library(mgcv)
library(terra)

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

# now we can import variables which we believe correlate with flowering. 

p <- '../data/spatial/processed'
preds <- rast(file.path(p, list.files(p)))

achy <- extract(preds, achy, bind = TRUE) |>
  st_as_sf()

## we don't have much scored data, so we will impute NA bulk density data ##
achy <- achy |>
  mutate(bulk_density = ifelse(is.na(bulk_density), 
                               mean(bulk_density, na.rm = TRUE), bulk_density), 
         across(Pct_Bud:Pct_Dropped, ~ ifelse(.x > 0, 1, .x)))


# models terms: s() means smoothing, we will not smooth a term like latitude

library(caret)

ctrl <- rfeControl(functions = gamFuncs,
                   method = "repeatedcv", number = 10,  repeats = 5,
                   verbose = FALSE, rerank = TRUE, 
                   allowParallel = TRUE)

inde <- as.matrix(st_drop_geometry(achy[,8:23]))
dep <- st_drop_geometry(achy$Pct_Anthesis)

cl <- parallel::makeCluster(parallel::detectCores(), type='PSOCK')
doParallel::registerDoParallel(cl)
lmProfile <- rfe(inde, dep, # y are outcomes
                 sizes  = c(1:10),
                 rank = TRUE, 
                 rfeControl = ctrl)
ParallelLogger::stopCluster(cl)


lmProfile[['optVariables']]


g_model <- gam(Pct_Anthesis ~ 
                 s(doy, bs = 'tp') + 
               #  s(gdgfgd10, bs = 'tp') + 
                 s(bio10, bs = 'tp'), 
            #     s(bio14, bs = 'tp') + 
            #     s(gdd0, bs = 'tp') + 
            #     s(cti, bs = 'tp') + 
            #     s(bulk_density, bs = 'tp') + 
             #    s(gdgfgd5, bs = 'tp'),
            na.action = 'na.omit', family = 'binomial', 
            data = achy, select=TRUE, method = 'ML', gamma = 7)

par(mar=c(4,4,2,2))
plot(g_model)
summary(g_model)
gam.check(g_model)
concurvity(g_model, full=TRUE) 
anova(g_model) 

g_model <- update(g_model, method='ML', na.action='na.fail')
summary(g_model)

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!

oz <- read.csv('https://raw.githubusercontent.com/selva86/datasets/master/ozone.csv')
oz.gam <- mgcv::gam(ozone_reading ~ s(Pressure_gradient, k=5, bs='tp') +
                s(Wind_speed, k=5, bs='tp') +
                s(Temperature_Sandburg, k=5, bs='tp'), data=oz,
              method='GCV', select=TRUE)

# METHODS GCV, AND AIC ARE BEST FOR ATTEMPTING TO PREDICT RESULTS. 

par(mar=c(4,4,2,2))
plot(oz.gam)
summary(oz.gam)

gam.check(oz.gam) # edfs approaching k' are problematic

concurvity(oz.gam, full=TRUE) # values approaching 0 are problematic

anova(oz.gam) # p value not as valuable as ML or REML. 

oz.ml <- update(oz.gam, method='ML', na.action='na.fail') # ML is used to compare models with different fixed effects - WILL NOT WORK FOR RANDOM EFFECTS. 

head(dredge(oz.ml, rank='AIC'),2)
