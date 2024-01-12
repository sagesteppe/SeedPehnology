library(MuMIn)
library(tidyverse)
library(mgcv)


setwd('~/Documents/SeedPhenology/scripts')
achy <- read.csv('../data/processed/high_priority_sheets.csv') %>% 
  filter(scientificname == 'achnatherum hymenoides')


# subset to scored sheets, and add '0' when a phenophasewas not observed. 
achy <- achy %>% 
  filter(! if_all(Pct_Bud:Pct_Dropped, ~ is.na(.))) |>
  mutate(across(Pct_Bud:Pct_Dropped, ~ replace_na(.x, 0)))

# now add spatial attributes to data. These will be used to look up the 
# relevant independent variables


# we now have Day of year, and latitude. DOY will be a fixed predictor, and I anticipate
# will hold the most explanatory power, where other terms will be covariates for it

# extract spatial terms: 

# models terms: s() means smoothing, we will not smooth a term like latitude, which 

(gam(flowering ∼s(doy, bs = 'tp'), family='binomial'))

# https://stacyderuiter.github.io/s245-notes-bookdown/gams-generalized-additive-models.html
# notes on GAMS!!

oz <- read.csv('https://raw.githubusercontent.com/selva86/datasets/master/ozone.csv')
oz.gam <- mgcv::gam(ozone_reading ~ s(Pressure_gradient, k=5, bs='tp') +
                s(Wind_speed, k=5, bs='tp') +
                s(Temperature_Sandburg, k=5, bs='tp'), data=oz,
              method='REML', select=TRUE)

par(mar=c(4,4,2,2))
plot(oz.gam)
summary(oz.gam)

gam.check(oz.gam) # edfs approaching k' are problematic

concurvity(oz.gam, full=TRUE) # values approaching 0 are problematic

anova(oz.gam) # p value not as valuable as ML or REML. 

oz.ml <- update(oz.gam, method='ML', na.action='na.fail') # ML is used to compare models with different fixed effects - WILL NOT WORK FOR RANDOM EFFECTS. 

head(dredge(oz.ml, rank='AIC'),2)
