library(terra)

setwd('~/Documents/SeedPhenology/scripts')

p <- '../data/spatial/processed'
preds <- rast(file.path(p, list.files(p, pattern = 'tif'))[1:3])

set.seed(27)
vars <- spatSample(preds, 15000)
vars <- vars[,2:17]
vars <- vars[ complete.cases(vars), ]
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

comp1 <- rast('../results/spatial/gddPCA-AB.tif')[[1]]

v <- values(comp1)
values(comp1) <- scales::rescale(v, to = c(0, 1))
plot(comp1[[1]])

writeRaster(comp1, '../results/spatial/gddPCA.tif', overwrite = TRUE)

rm(v, p, preds, comp1)

# first raster is a little big... 8.5m cells. 
# now bring down the resolution so we can interpolate more quickly. 

comp1 <- rast( '../results/spatial/gddPCA.tif')
cp1_c <- dim(comp1)[1] * dim(comp1)[2] 

comp2 <- aggregate(comp1, factor = 2, fun = 'mean')  
dim(comp2)[1] * dim(comp2)[2] / cp1_c # one quarter the size. ;-)

writeRaster(comp2, '../results/spatial/gddPCA.tif', overwrite = TRUE)

rm(comp1, comp2, cp1_c)

