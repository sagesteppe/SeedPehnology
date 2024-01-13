library(sf)
library(tidyverse)
library(mgcv)
library(terra)

setwd('~/Documents/SeedPhenology/scripts')
achy <- read.csv('../data/processed/high_priority_sheets.csv') %>% 
  filter(scientificname == 'achnatherum hymenoides')

names(preds)
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
  mutate(
    across(
      doy:cti, ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
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
rm(cl)

terms <- unique(c('doy', lmProfile[['optVariables']]))

formula <- as.formula(
  paste(
    "Pct_Anthesis ~ ",  paste0("s(", terms, ", bs = 'tp')", collapse = " + ")
    )
)

g_model <- gam(formula,
            na.action = 'na.omit', family = 'binomial', 
            data = achy, select = TRUE, method = 'ML')

par(mar=c(4,4,2,2))
plot(g_model)
summary(g_model)
gam.check(g_model)
concurvity(g_model, full=TRUE) 
anova(g_model) 

g_model <- update(g_model, method='ML', na.action='na.fail')
summary(g_model)


rm(terms, formula)
# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!
