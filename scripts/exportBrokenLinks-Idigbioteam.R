library(tidyverse)
setwd('/home/sagesteppe/Documents/SeedPhenology/scripts')

scored <- bind_rows(
  read.csv('../data/processed/high_priority_sheets-scored.csv'),
  read.csv('../data/processed/second_phenology_review-scored.csv'), 
  read.csv('../data/processed/low_priority_sheets-scored.csv'), 
  read.csv('../data/processed/troublesomeSpecies-scored.csv')
) %>% 
  filter(str_detect(accessr, 'https://www.pnwherbaria.org/images/j'))  


scored <- scored[is.na(scored$Anthesis), c('scntfcnm', 'accessr')]


write.csv(scored, '../data/brokenCPNWH.csv', row.names = F)
