library(terra)

setwd('~/Documents/SeedPhenology/scripts')


p <- '../data/spatial/processed'
preds <- rast(file.path(p, list.files(p)))

set.seed(27)
vars <- spatSample(preds, 10000)
vars <- vars[,2:14]
pca <- prcomp(vars)

plot(pca)

factoextra::fviz_eig(pca) # very simple ! all variables so correlated 
# the first dimension contains all of the info !
factoextra::fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

predict(preds, pca, cores = parallel::detectCores(),
             filename = file.path('../results/spatial/gddPCA-AB.tif'))

comp1 <- rast('../results/spatial/gddPCA-AB.tif')
writeRaster(comp1, '../results/spatial/gddPCA.tif')

plot(comp1)
