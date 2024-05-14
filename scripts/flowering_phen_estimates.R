library(tidyverse)
library(sf)
library(lubridate)

setwd('/media/steppe/hdd/SeedPhenology/scripts')

scored <- bind_rows(
  read.csv('../data/processed/high_priority_sheets-scored.csv'),
  read.csv('../data/processed/second_phenology_review-scored.csv'), 
  read.csv('../data/processed/low_priority_sheets-scored.csv'), 
  read.csv('../data/processed/troublesomeSpecies-scored.csv')
) %>%            
  mutate(
    Anthesis = if_else(is.na(Anthesis), 999, Anthesis)
  ) %>% 
  select(-comments)

specimens <- st_read('../data/processed/training/presences.shp', quiet = TRUE ) %>% 
  filter(! scntfcnm %in% scored[scored$Anthesis %in% c(0, 999), 'accessr'] )
lp <- st_read('../data/processed/lowest_priority_sheets/lowest_priority_sheets.shp', quiet = TRUE) %>% 
  group_by(scntfcnm) %>% # these records were never cleaned
  arrange(doy) %>% 
  slice_head(prop = 0.9) %>% 
  slice_tail(prop = 0.9)

records <- bind_rows(specimens, lp)
rm(scored, lp, specimens)

# reduce window to when more than a few populations are flowering

breaks <- c(1, 32, 61, 92, 121, 152, 182, 213, 244, 275, 306, 336)
labels <- c('Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'June', 'July', 'Aug.', 'Sept.', 'Oct.', 'Nov.', 'Dec.')

spliecies <- split(records, records$scntfcnm)

phen_maker_fn <- function(x){
  
  species <- gsub('_', ' ', x$scntfcnm[1])
  
  mindoy <- (min(x$doy) - 28) ; maxdoy <- (max(x$doy) + 28)
  mindoy <- if(mindoy <= 0){mindoy <- 1} else{mindoy}
  maxdoy <- if(maxdoy >= 365){maxdoy <- 365} else {maxdoy}
  
  ggplot(x) +
    geom_density(aes(x = doy), fill = '#FF7B9C', color = '#FFC759', alpha = 0.5) +
    theme_bw() + 
    labs(title = 'Estimated Flowering') + 
    theme(aspect.ratio = 6/16, 
          plot.title = element_text(
            hjust = 0.5, colour = "black", size = 7, face = "bold"),
          axis.text.x= element_text(
            family = "Tahoma", face = "bold", colour = "black", size=5),
          panel.background = element_rect(fill='#607196'),
          plot.background = element_rect(fill='#607196'),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(breaks = breaks, labels=labels, limits = c(mindoy, maxdoy) )
  
  ggsave(filename = paste0('../results/phen/', species, '.png'), plot = last_plot(), 
         dpi = 150, width = 520, height = 300,  units = "px",  bg = 'transparent')
}

lapply(spliecies, phen_maker_fn)
