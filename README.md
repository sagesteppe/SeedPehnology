# SeedPhenology
Estimating the initiation of flowering, peak flower, and the cessation of flowering using herbarium specimens to inform wildland seed collection in the Western United States. 

Considerable effort has been directed towards understanding how climate change affects the phenology of species.
These studies have shown causal links between a number of climate metrics, and generally, advanced initiation of phenological events, which has resulted in the generation of a variety of theory and spatial data which serve as independent variables. 
However, these observational studies and manipulative experiments have generally been limited to a few dozen species in only one to a couple populations, or when using herbarium sheets from many populations than only a few species. 
Hence their findings have seldom been applied to the documentation of phenological traits for species in general. 

A component missing from most of these studies has been the utilization of spatial covariates in the workflows and statistics used to model phenological processes. 
We believe that the inability to incorporate spatial terms in the modelling process limits the documentation of phenology to only a handful of populations with similar environmental conditions, rather than reflecting the species as a whole. 
Generalized Additive Models (GAM's) are often used to document a phenophase because a single model can have their splines fit to both initiation, peak, and cessation of an phase, using a single or multiple independent variables - a limitation of several other methods of estimation.
The use of independent variable(s) alongside GAM's ability to incorporate an error-correlation structure which accommodates spatial autocorrelation allows them to model the phenological parameters of a species across it's geographic and concomitant environmental range. 

Hererin we use GAM's to model phenophases, inferred from herbarium specimens and using environmental predictors identified as casual cues of phenology, in space. 
Our necessity to more accurately understand the phenology of species arose from our goal of native seed collection for both native plant germplasm development, and *ex-situ* conservation. 
The identification of putative populations with enough individuals to warrant germplasm development is a time consuming process, because most plant species can only be identified when they have reproductive organs, and a populations ability to support these collections varies wildly with the years weather, pathogen load, and various stochastic processes. 
The collection of seeds, which is generally occurring for both many species and many populations each year, is challenged by both the need for crews to collect from other species and the natural dispersal of seeds - simply put *timing is everything*.  

# Methods

## Data Sources

Species records were derived from the Symbiota herbarium portal for all years from 1981-2021, these years reflected the climate means used as independent variables. 
All records were downloaded, and the records in the 2.5% Day of Year (DOY) quantile were manually reviewed. 
These early records were reviewed because, especially with graminoids, collectors may actually collect material without reproductive organs yet reaching anthesis.
The later records were reviewed because collectors may have collections of individuals entirely post-anthesis - a situation very common with certain clades where species are distinguished by characteristics of their fruits (e.g. Leguminosae). 
In both scenarios analysts proceeded towards the mean of the distribution until they encountered 5 consecutive sheets with the desired phenophase. 

CHELSA climate variables for Growing Degree Days (GDD) heat sums (at 0C, 5C, 10C), first (gdgfgd) and last (gddlgd) growing day of year, vapor pressure deficits (vpd), Bio10 (mean daily mean air temperatures of the warmest quarter), and Bio14 (precipitation amount of the driest month). 
Soil bulk density, which is shown to reflect the amount of air/water space in soil, was downloaded from SoilGrids... 
Compound Topographic index, which describes the potential of areas throughout a landscape to accumulate soil moisture via a combination of its landform position, slope, aspect, and upslope catchment area, was downloaded from geomorpho90m and resampled from 90m to the 250m resolution of the previous data sets. 

GAM's require data on when a species was **not** flowering in order to develop splines for the onset of flowering. 
Pseudo-floral absences were created using known sites, and their observed phenology. 
All of the CHELSA climate variables were decomposed using PCA, and the first axis (explaining 97.6% of the variation; 750m x 750M CELL RESOLUTION?) was used as a feature space in a Ward-like hierarchical clustering algorithm which seeks to maximize homogeneity of both the feature and constraint (here geography) space (hclustgeo).  
A suitable number of clusters from the independent variable were automatically selected using kgs (maptree), these clusters were then reanalyzed in light of the constraint space using automatic selection of an alpha parameter which blends the feature and constraint space and re-clustered using hclustgeo (Clustgeo).  
Each cluster had weibull estimates of flowering initiation and cessation modelled, and any doy within 28 days preceding onset or following cessation were drawn for each group (phenesse). 
These values were arranged by ascending day of year and joined to the members of the group via decreasing warm to cool values along the PCA axis. 
Points in clusters which had a nearest geographic neighbor in another cluster had their randomly generated pseudo-absences wiped, and thin plate spline regression using the PCA axis as an independent variable, and interpolation was used to repopulate the floral pseudo-absence dates (fields, terra). 

## Modelling

All independent variables were extracted to the dependent variables, and if a value for an independent variable was missing - which was not uncommon for Soil Bulk Density, where the modellers excluded the fringes of several vernally wet playas - it was imputed as the mean of the variable for the species. 
All independent variables then underwent feature selection using the Recursive Feature Elimination (rfe) with 10 Cross-Validations (CV) folds, 5 replicates, and from 1-10 variables using caret (@kuhn).
The remaining variable(s) were used as covariates with day of year always included in the models. 
A gamm was fit using presence/absence of flowering as a response, as well as gamm's with error structure of gaussian, spherical, and exponential variograms, with REML (@package). 
All models were subjected to model selection, and the top model determined via AUC scores (MuMIn).. 

## Surfaces

Convex hulls were drawn around the geographic extent of the species, and a prediction was generated for each DOY with greater than 5% probability of flowering... 
