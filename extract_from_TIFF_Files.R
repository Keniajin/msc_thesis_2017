#########################################################################################################################################################
################## Script for extracting MODIS vegetation indices for analysis##########################################################################
#######################   processing MODIS EVI to get EVI data   ###########################################################################################
##***********************     remotely-sensed data in R.    #############################################################################################
##**************************    Author: KMwai Feb 2016  ################################################################################################
set.seed(1221223)
library(MODIS)
library(dplyr)
library(foreign)
library(tidyr)

##read in the data with cordinates from KHDSS
##  1km by 1km data
dat <- read.dta("mogeni/KHDSS.DTA")

## read in the admissions data with admission cordinates
dataAdm <- read.csv("data/admision_lat_long_sec.csv")

#'         \\\_Extract data process_\\\         #
#'#'\\_________Starting the extract of data_______________________\\
#'(this process takes some time and needs 30 GB of space free)

##pull from
# # Custom command. You could build your own using a loop and paste
# f <- paste0("modis_download.py -r -p MOD13Q1.005 -t ",tile," -e ","2014-05-01"," -f ","2014-05-30"," tmp/","h21v09")
# # Call the python script
# system(f)
# # Then go into your defined download folder (tmp/ in our case)


##identify the folder with the monthly  tiff data for EVI
vi <- preStack(path = "modis/monthly_data/", pattern = "*.tif$")

### stack the data to data frame
s <- stack(vi)
s <-  s * 0.0001 # Rescale the downloaded Files with the scaling factor (from modis) EVI

#'#'\\_________extracting for the 1km by 1km data_________\\
# And extract the mean value for our point from before.
# First Transform our coordinates from lat-long to to the MODIS sinus Projection
## method='bilinear' used for extraction 
## If 'simple' values for the cell a point falls in are returned.
## If 'bilinear' the returned values are interpolated from the values of the four nearest raster cells.

sp <-  SpatialPoints(coords = cbind(dat[2], dat[1]),
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") )
sp <- spTransform(sp, CRS(proj4string(s)))
dataExtract <- raster::extract(s, sp , df=T, method="bilinear") # Extract the EVI
write.csv(dataExtract , "data/modis_data_2001_2015.csv")

#'#'\\_________extracting for the admissions data_________\\
sp2 <-  SpatialPoints(coords = cbind(dataAdm[2], dataAdm[1]),
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
sp2 <- spTransform(sp2, CRS(proj4string(s)))
dataExtractAdm <- raster::extract(s, sp2 , df=T, method="bilinear") # Extract the EVI
write.csv(dataExtractAdm , "data/dataExtractADM.csv")
#'@_______________________________________________________________
#'         \\\_End Extract data process_\\\         #
#'         

#'#'\\_________extracting for the sublocation centroids_________\\
#'## read in the sublocations centroids
pointsKilifi <- read.csv("data/satScan/points_kilifi.csv")
pointsKilifi <- pointsKilifi %>%  select(lat , long , SLNAME)
sp2 <-  SpatialPoints(coords = cbind(pointsKilifi[2], pointsKilifi[1]),
                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
sp2 <- spTransform(sp2, CRS(proj4string(s2)))
dataExtractPointsKilifi <- raster::extract(s2, sp2 , df=T, method="bilinear") # Extract the EVI
write.csv(dataExtractPointsKilifi , "data/pointsKilifi.csv")
#'@_______________________________________________________________
#'         \\\_End Extract data process_\\\         #         
       

####
#'#'\\_________Extracting the rainfall data_________\\

###Extracting the rainfall data

##identify the folder with the monthly  tiff data for EVI
vi2 <- preStack(path = "rainfall/", pattern = "resampledchirps-v2.0.20*")

### stack the data to data frame
si2 <- stack(vi2)

#'#'\\_________extracting for the admissions data_________\\
sp2_b <-  SpatialPoints(coords = cbind(dataAdm[2], dataAdm[1]),
                        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
sp2_b <- spTransform(sp2_b, CRS(proj4string(si2)))
dataExtractAdm2 <- raster::extract(si2, sp2_b , df=T, method="bilinear") # Extract the rainfall
write.csv(dataExtractAdm2 , "dataExtractADM_rainfall.csv")
#'@_______________________________________________________________
###change the path name for temporary file
#load("I:/Project/admissionExtract.Rdata")
#s2@file@name <- "G:\\Rtmp8WMcwe\\raster\\r_tmp_2017-03-07_094800_1668_02103.gr"

# library("multidplyr")
# cluster <- create_cluster(2)
# #> Initialising 2 core cluster.
# cluster



