# SeedPehnology
Estimating the initiation of flower, peak flower, fruit maturation, and fruit dispersal using herbarium specimens to inform wildland seed collection in the Western United States. 


## Workflow for Neural Net Training

Mask R-CNN models have been successfully applied to  determining and counting the phenological status of reproductive organs on herbarium sheets.
Previous efforts have largely focused on Mid and Upper Atlantic flowering herbs, generally monocots, with showy flowers, and very large samples. 
These efforts have been published in the peer-reviewed literature, and we believe can be operationalized on a larger scale, albeit with less accurate results. 

Given the phylogenetic and accompanying morphological diversity of the target taxa, new mask R-CNN models will be generated for several clades. 
...

The workflow for mask r-cnn training requires the use of many humans to label image data. 
Best practice involves sampling a diverse set of images which feature a diversity of reproductive states across the range of the species (for statistical modelling, see OTHER SECTION). 


Generation of a training image set for MASK R-CNN includes:

1) Use SEINET to access all occurrence records, and associated image data in text format.
2) Using all occurrence records extract the latitude, longitude, and Day of Year (DOY) the specimen was collected.
3) Rescale the above variables to the ranges of 0-1
4) From an orthogonal latin hypercube design, gather 100 points. 
5) From each of the 100 points from the LHS, identify the scanned herbarium record, with the smallest distance.
6) Identify the next two closest herbarium records, and randomly select 50 records for over-sample records. 
7) use wget to download all relevant images for distribution and annotation
8) Split the 100 specimens and their associated over-samples into 3 30 specimen data sets, each for a separate analyst. 
9) A Senior botanist will annotate 10 records as a guide which the analysts can refer to as they work. 

