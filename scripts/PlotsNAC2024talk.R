library(tidyverse)
library(sf)

setwd('~/Documents/SeedPhenology/scripts')

modelled <- list.files('../data/processed/timestamps/')

t <- read.csv('../data/SOS/Scouting.csv') %>%
  select(taxa, scoutDate, percentDormant:percentPost, LATITUDE_DECIMAL, LONGITUDE_DECIMAL) |>
  st_as_sf(coords = c(x = 'LONGITUDE_DECIMAL', y = 'LATITUDE_DECIMAL'), crs = 4326) |>
  st_transform(5070) |>
  mutate(taxa = str_replace_all(taxa, ' ', '_')) |>
  mutate(
    taxa = case_when(
    taxa == 'Machaeranthera_canescens' ~ 'Dieteria_canescens',
    taxa != 'Machaeranthera_canescens' ~ taxa
  ),
  scoutDate = lubridate::yday( as.Date(scoutDate, tryFormats = c('%m/%d/%Y'))),
  sDate_join = as.numeric(scoutDate), .before = geometry) |>
  filter(taxa %in% modelled) %>%
  # we will remove Phacelia crenulata and Eriogonum fusiforme which did not flower this year #
  filter(!taxa %in% c('Phacelia_crenulata', 'Eriogonum_fusiforme', 'Eriogonum_shockleyi'))

# load and extract raster values to ground verified point.

extractPredictions <- function(x){

  taxon <- sf::st_drop_geometry(x$taxa[1])
  fp <- paste0('../data/processed/timestamps/', taxon, '/doy_preds')
  doys <- data.frame(
    sDate_join = as.numeric(gsub('.tif', '', list.files(fp))),
    path = paste0(fp, '/', list.files(fp))
  )

  # each species had its own biweekly intervals generated, these do not align.
  # we will add the nearest DOY prediction to each file. If the nearest date is #
  # before the start, or after the end. We will convert these into predictions of
  # 0 for the phenophase.

  data.table::setDT(doys); data.table::setkey(doys, sDate_join)
  data.table::setDT(x); data.table::setkey(x, sDate_join)
  set_merged <- doys[x, roll = "nearest"]

  set_merged <- set_merged |>
    dplyr::mutate(
      sDate_join = as.numeric(gsub('.tif', '', basename(path))),
      path = dplyr::if_else(abs(scoutDate - sDate_join) >= 14, NA, path)) |>
    sf::st_as_sf()

  # now split each of these objects by the raster for the date we will be comparing them to.
  # we can load the raster once, and compare it to several scouting points simultaneously.

  NAindexes <- which(is.na(set_merged$path))

  if(length(NAindexes)>0){
    noRasters <- set_merged[NAindexes,]
    set_merged <- set_merged[-NAindexes,]
  }

  # and now run the process for the populated cells.
  obbies <- split(set_merged, f = set_merged$sDate_join)

  extractByDate <- function(x){

    Path <- sf::st_drop_geometry(x) |>
      dplyr::pull(path) |>
      unique()

    surface <- terra::rast(Path)
    extracted_values <- terra::extract(surface, terra::vect(x)) |>
      dplyr::select(-ID)

    extracted_values <- data.frame(
      doy = colnames(extracted_values),
      prob = extracted_values[,1]
    )


    # combbine the extracted values with the features of the input data set.
    x <- dplyr::bind_cols(x, extracted_values)

    # if any of the rows are NA, use a feature engineered probability of flowering for them

    if(any(is.na(x$prob))){

      # now proceed row wise
      need_imputed <- x[which(is.na(x$prob)),]

      for(i in 1:nrow(need_imputed)){
      # mask any raster sites > 10k from population #
      mask <- sf::st_buffer(need_imputed[i,], 10000)
      surface_m <- terra::mask(surface, mask, inverse = FALSE)
      m <- terra::global(surface_m, mean, na.rm = TRUE)
      need_imputed[i,'prob'] <- m
      }

      x <- dplyr::bind_rows(
        x[which(!is.na(x$prob)),],
        need_imputed)

    }
      return(x)
  }


  ebd <- lapply(obbies, extractByDate)
  ebd <- dplyr::bind_rows(ebd)
  if(exists('noRasters')){

    if(nrow(ebd)> 0){
      noRasters <- dplyr::bind_cols(
        noRasters,
        data.frame(
          doy = as.character(noRasters$sDate_join),
          prob = 999))
    ebd <- dplyr::bind_rows(ebd, noRasters) |>
      dplyr::arrange(scoutDate)} else {
      ebd <- noRasters
      }}

  # now make values outside of the range 0-1 - 0 these are our not-flowering dates.
  return(ebd)

}

dat <- split(t, f = t$taxa)
ob <- lapply(dat, extractPredictions) |>
  bind_rows()

ob1 <- select(ob, -path, -sDate_join) |>
  sf::st_drop_geometry() |>
  rename_with( ~ str_replace(., 'percent', 'pct')) |>
  mutate(
    prob = if_else(prob > 1| is.na(prob), 0, prob),
    Class = case_when(
      prob > 0.5 & pctFlower > 0 ~ TRUE, 
      prob < 0.5 & pctFlower == 0 ~ TRUE,
      prob > 0.5 & pctFlower == 0 ~ FALSE
    ),
    Class2 = case_when(
      prob > 0.5 & pctFlower > 0 | prob > 0.5 & pctBud > 0 | prob > 0.5 & pctPre > 0 ~ TRUE, 
      prob > 0.5 & pctBud > 0 & pctPre > 0  ~ TRUE, 
      prob > 0.5 & pctBud > 0 & pctSeed > 0 ~ TRUE, 
      prob < 0.5 & pctFlower == 0 ~ TRUE,
      prob > 0.5 & pctFlower == 0 ~ FALSE
    )
  )

#  write.csv(ob1, '../data/SOS/GroundVerified.csv')
ob1 <- read.csv('../data/SOS/GroundVerified.csv') |>
  mutate(
    ExtendedClass = case_when(
      prob > 0.5 & pctFlower > 0 ~ 'TP', 
      prob < 0.5 & pctFlower == 0 ~ 'TN',
      prob > 0.5 & pctFlower == 0 ~ 'FP',
      prob < 0.5 & pctFlower > 0 ~ 'FN',
    ),
    ExtendedClass2 = case_when(
      prob > 0.5 & pctFlower > 0 ~ 'TP', 
      prob > 0.5 & pctFlower > 0 | prob > 0.5 & pctBud > 0 | prob > 0.5 & pctPre > 0 ~ 'TP', 
      prob < 0.5 & pctFlower == 0 ~ 'TN',
      prob > 0.5 & pctFlower == 0 ~ 'FP',
      prob < 0.5 & pctFlower > 0 ~ 'FN',
    ))

table(ob1$ExtendedClass)
table(ob1$ExtendedClass2)

#' @param x a dataframe containing the results
#' @param y column name containing the classification of results
metrics <- function(x, y){
  
  ct <- data.frame(table(x[,y]))
  colnames(ct) <- c('Var', 'Freq')
  
  Sensitivity <- ct[ct$Var=='TP','Freq'] / (ct[ct$Var=='TP','Freq'] + ct[ct$Var=='FN','Freq'])
  Specificity <- ct[ct$Var=='TN','Freq'] /  (ct[ct$Var=='TN','Freq'] + ct[ct$Var=='FP','Freq'])
  Accuracy <- (ct[ct$Var=='TN','Freq'] + ct[ct$Var=='TP','Freq']) / sum(ct$Freq)
  
  metrics <- data.frame(
    Var = c('Sensitivity', 'Specificity', 'Accuracy'), 
    Freq = round(c(Sensitivity, Specificity, Accuracy), 3)
  )
  
  ct <- rbind(ct, metrics)
  lvls <- c('TP', 'FN', 'Sensitivity', 'FP', 'TN', 'Specificity', 'Accuracy')
  ct <- dplyr::arrange(ct, factor(Var, levels = lvls))
  
  positions <- data.frame(
    x = c(rep(c(1, 2), each = 3), 3), 
    y = c(rep(c(3, 2, 1), times = 2), 1)
  )
  ct <- cbind(ct, positions)
  
  return(ct)
  
}

ec_metrics <- metrics(ob1, 'ExtendedClass') 
ec_metrics2 <- metrics(ob1, 'ExtendedClass2') 

contig_example <- data.frame(
  x = rep(c(1, 2, 3), each = 3), 
  y = rep(c(3, 2, 1), times = 3), 
  values = c(
    'True Positive\n(TP)', 'False Negative\n(FN)', 'Sensitivity', 
    'False Positive\n(FN)', 'True Negative\n(TN)', 'Specificity',
    'Positive\nPredictive\nValue (PPV)', 'Negative\nPredictive\nValue (NPV)', 'Accuracy'
    ),
  fill = c(
    'Green', 'Red', 'grey90',
    'Yellow', 'Green', 'grey90',
    'grey90', 'grey90', 'grey90'
    )
)

fill_lkp <- c('Green' = '#093824', 'Red' = 'firebrick1',
              'Yellow' = '#d95f02', 'grey90' = 'transparent')


  
p <- ggplot(contig_example, aes(x, y, label = values, fill = fill)) + 
  geom_tile(colour = "grey50") + 
  scale_fill_manual(values = fill_lkp) + 
  theme_void() + 
  scale_x_continuous(
    breaks = 1:3,
    position = 'top',
    labels = c('Positive', 'Negative', 'Metrics')) + 
  scale_y_continuous(
    breaks = 3:1,
    labels = c('Positive', 'Negative', 'Metrics')) +
  labs(x = 'Observed Class', y = 'Predicted Class', title = 'Contingency Table') + 
  theme( 
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10, angle = 90), 
    axis.text.x = element_text(size = 10), 
    axis.title.x = element_text(
      size = 12, hjust = 0.25,
      margin = margin(t = 10, r = 0, b = 5, l = 0)
      ), 
    axis.title.y = element_text(
      size = 12, angle = 90, hjust = 0.75,
      margin = margin(t = 0, r = 5, b = 0, l = 5)
      ), 
    legend.position = 'none', 
    text = element_text(color = 'white'),
    panel.border = element_blank(),
    panel.background = element_rect(fill='transparent', color =NA), 
    plot.background = element_rect(fill='transparent', color=NA)
    ) + 
  geom_segment(aes(0.5, 3.5, xend = 2.5, yend = 3.5), lwd = 1.1, color = '#7570b3') + 
  geom_segment(aes(0.5, 2.5, xend = 2.5, yend = 2.5), lwd = 1.1, color = '#7570b3') +
  geom_segment(aes(0.5, 1.5, xend = 2.5, yend = 1.5), lwd = 1.1, color = '#7570b3') +
  
  geom_segment(aes(0.5, 3.5, xend = 0.5, yend = 1.5), lwd = 1.1, color = '#7570b3') + 
  geom_segment(aes(1.5, 3.5, xend = 1.5, yend = 1.5), lwd = 1.1, color = '#7570b3') + 
  geom_segment(aes(2.5, 3.5, xend = 2.5, yend = 1.5), lwd = 1.1, color = '#7570b3') 
  

explain <- p + geom_text(color = 'white') + 
  
  annotate(
    "text", x = 1, y = 0.7, size = 2.5, parse = TRUE, color = 'white', # Sensitivity 
    label = latex2exp::TeX('\\frac{TP}{TP+FN}', output = "character")) + 
  annotate(
    "text",  x = 2, y = 0.7, size = 2.5, parse = TRUE, color = 'white',# Specificity 
    label = latex2exp::TeX('\\frac{TN}{TN+FP}', output = "character")) + 
  annotate(
    "text", x = 3, y = 0.7, size = 2.5, parse = TRUE, color = 'white', # Accuracy 
    label = latex2exp::TeX('\\frac{TP + TN}{TP+TN+FP+FN}', output = "character")) 

exact <- p + 
  geom_text(data = ec_metrics, aes(label = Freq, x = x, y = y),
            color = 'white', inherit.aes = F) + 
  labs(title = 'Observed') 

realistic <- p + 
  geom_text(data = ec_metrics2, aes(label = Freq, x = x, y = y), 
            color = 'white', inherit.aes = F)+ 
  labs(title = 'Realistic') 

ggsave(filename = '../plots/explanatory.png', 
       plot = explain, units = 'in', width = 5, height = 5)
ggsave(filename = '../plots/exact.png', 
       plot = exact, units = 'in', width = 5, height = 5)
ggsave(filename = '../plots/realistic.png', 
       plot = realistic, units = 'in', width = 5, height = 5)

