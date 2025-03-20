#============= Running Bfast Model ===================#
#  Location  Venejuala                                #
#  Code improved by Mukti Subedi                      #
#   Date 1/29/2019                                    #
#=====================================================#

####- Setting options -####
# ------------------------------------------------------------------------------
options(timeout   = 15,     prompt    = "< ",
        continue  = "...",  digits    = 5,
        width     =  80,    papersize = "letter" )

####- Clear workspace -####
rm(list = ls())

## loading bfast spatial packages
## requires devtools because the bfast spatial is not available in CRAN

####- Install devtools -####
# if you have already installed bfast spatial following line of code will
# check and verify that you have upto date package. If there is a update then
# it will install the updated packages
# devtools::install_github('loicdtx/bfastSpatial', ref = 'develop')

####- Loading Packages -####
library(bfastSpatial)
library(parallel)      ## does not require in windows
library(doParallel)    ## does not require in windows
library(raster)

####- Setting up Directory -####
(setwd("C:\\NepImage")) # Change this in your computer
#------------------------------------------------------------------------------#
(path         <- getwd())
(tmpDir       <- rasterOptions()$tmpdir) # you can create your own temp directory
(inDir        <- file.path(path,    'Data'))     # path to Data directory
(stepDir      <- file.path(inDir,   'Datastep')) # path to Datasetp
(landsat7Dir  <- file.path(stepDir, 'Landsat7')) # path to Landsat7 Data
(landsat8Dir  <- file.path(stepDir, 'Landsat8')) # path to landsat8 Data
(outDir       <- file.path(path,    'Out'))       # path to Output   Data

# These folders were created in windows explorer and if you decided to give
# different folder name it is oaky to do so however make sure it is consistent
# through our your code. :-)
# if you want to create a folder from R you can use dir.create(path)
### Setting directory for nbr2 data
(nbr2Dir      <- file.path(stepDir, 'nbr2'))
#------------------------------------------------------------------------------#
rasterOptions(progress  = "text") # it will show update in the raster processing
rasterOptions(overwrite = TRUE)   # to save space it it is  a good idea to
                                  # delete and overwrite files
#------------- this part can be ignored in windows ----------------------------#
####- Raster Processing setting -####
# look at rasterOptions
# There are two parameters that have a lot of influence on the performance of
# the {raster} package: chunksize and maxmemory.
# chunksize: integer. Maximum number of cells to read/write in a single
# chunk while processing (chunk by chunk) disk based Raster* objects
# maxmemory: integer. Maximum number of cells to read into memory.

(rasterOptions())
### lets increase the chunkzie and maxmemory by 10 times
(rasterOptions(chunksize = 5e+09)) # default 1e+07
(rasterOptions(maxmemory = 9e+10)) # default 1e+08

### detect how many cores your computer has
(detectCores())

# by default R will use 1 core we can change it to higher number (max -1)

(use.core<- detectCores()-1) # However, multicore is not allowed in windows
# in my machine that I wrote code on has 8 cores
# normal laptops will have 4 cores
# check your computer how many cores you have
## now to perfom parallel processing we have to define how many cores we want
# to use (use. core) object showed we wanted to use total -1 = 7 in this
# computer # The function registerDoparallel()  below to register your cores
# for the parallelization using registering CoreCluster

(clster<- makeCluster(use.core))

registerDoParallel(clster)
## to end the cluster use call stopCluster()
## In my computer (Windows) the above code doesnt make any difference and only
## one core will be uzed while processing data
#------------------------------------------------------------------------------#
####- Batch Processing Landsat7 -####
#### Discription:
# Processes Landsat data (a single scene for processLandsat and multiple
# scenes for processLandsatBatch), from tarball or directory containing surface
# reflectance bands to vegetation index. The data must be surface reflectance
# obtained from espa (https://espa.cr.usgs.gov/). These data may already
# contained pre-processed indices layers, in which case they are directly used.
# The batcher (processLandsatBatch) allows to process multiple scenes with one
# command; sequentially or in parallel (parallel processing only work on unix
#                                       systems (linux and mac)).

#### Usage:
#
# processLandsat(x, outdir, vi = "ndvi", srdir = NULL, delete = FALSE,
#                mask = NULL, keep = c(0), e = NULL, fileExt = "grd",
#                overwrite = FALSE)

# printing time before processing landsat7 images
#==============================================================================#
#==============================================================================#
(land7Time1<- Sys.time())
#------------------------------------------------------------------------------#
processLandsatBatch(x         = landsat7Dir, # directory to load data
                    outdir    = stepDir,     # output directory
                    vi        = 'nbr2',      # vegetation index
                    delete    = TRUE,        # deletes temporary files
                    mask      = 'pixel_qa',  # selects those pixles only
                    keep      = c(66, 130),  # cloudless pixels values
                    overwrite = TRUE,       # overwrite files if they already
                    mc.cores  = 1)
                                             # exists in the directory
#------------------------------------------------------------------------------#
# printing time after processing landsat7 images

(land7Time2<- Sys.time())

## printing time differrence
## This let us understand how much time it took to process
## landsat 7 images
##
(landsatProcessTime<- land7Time2 - land7Time1)

####- Batch Processing Landsat8 -####
# printing time before processing landsat8 images
(land8Time1<- Sys.time())
#------------------------------------------------------------------------------#
processLandsatBatch(x         = landsat8Dir,
                    outdir    = stepDir,
                    vi        = 'nbr2',
                    delete    = TRUE,
                    mask      = 'pixel_qa',
                    keep      = c(322, 386,834,898,1346), # see link below
                    overwrite = TRUE,
                    mc.cores = 1)

# https://landsat.usgs.gov/documents/lasrc_product_guide.pdf
(land8Time2<- Sys.time())
#------------------------------------------------------------------------------#
(landsatProcessTime8<- land8Time2 - land8Time1)
## printing time differrence
## This let us understand how much time it took to process
## landsat 7 images
##

#-----------------------SKIP, if TIMESTACK Works ------------------------------#
#------------------------------------------------------------------------------#
# If landsat full scenes are used from several landsat products or even with
# same sensor due to radiometric/geometric error the exact path may not align,
# however, timeStack() doesn't work if the areal extent is not same for all the
# layers you want to stck. The trouble starts here.
#
#                 HOW DO WE SOLVE THIS PROBLEM?
#------------------------------------------------------------------------------#
# Actually there could be several options. I find it easier to use extents of
# already processed images and find the common extent and use this extent as a
# "e" function in processLandsatBatch() command. below is the example i used for
# five scenes. DONT FORGET TO CHANGE THE CODE AS PER YOUR REQUIREMENT
#------------------------------------------------------------------------------#
####- Taking care of common extent -####
#------------------------------------------------------------------------------#
#    reading nbr2 data from nbr2 folder as AllData
###  The end
###  Next improvement on "Spatial BFASTMonitor" :-)
###  Import shapefile.

#------------------------------------------------------------------------------#
####- Taking care of common extent1 -####
#------------------------------------------------------------------------------#
# reading nbr2 data from nbr2 folder as AllData
AllData <- list.files(file.path(stepDir, 'nbr2'),
                      full.names = TRUE,
                      pattern = ".gri$")
(AllData)
# creating list of names for extent
xmincol<-    c("MinX", "MinY", "MaxX", "MaxY")
# creating data frame for extent values
data<-        data.frame(matrix(ncol = 4, nrow = 0))
# applying names to dataframe
names(data)<- xmincol

# looping through gri filies from nbr2 folder and reading min and max of x and y
# and storing them in data frame created earlier
for(i in AllData) {
  x<- assign(unlist(strsplit(i, "[.]"))[1],  raster(i))
  data[i,1]<- xmin(x)
  data[i,2]<- ymin(x)
  data[i,3]<- xmax(x)
  data[i,4]<- ymax(x)
}

# we need a maximum of minX and everyother  data should be minimum for all the
# layers to have a common extent
(x.min<- max(data$MinX))
(y.min<- max(data$MinY))
(x.max<- min(data$MaxX))
(y.max<- min(data$MaxY))

# extent to be applied in batch processing
(exten<- c(x.min,x.max, y.min, y.max)) # x.max = 1968285

# define extent for batch processing
(extentForBatch<- extent(exten))
# making a raster object  based on extent with spatial resolution
# of 30 by 30 m

 extenras<-raster(extentForBatch)
 #crs(extenras)<-crs(exten)
# applying projection/cooridnate system to exten objected based on
# original data

# projection(extenras)<-crs(x) # remember we have lists of data on x list

########## NOW STOP HERE AND APPLY exten OBJECT ON BATCH PROCESS
#==============================================================================#
#         ####- Running Batch Process with Extent -####                        #
#==============================================================================#
# running batch processing again. Remember once you know the common extention of
# your data you don't have to bother to run the process twice.
#------------------------------------------------------------------------------#

(land7Time1<- Sys.time())
#------------------------------------------------------------------------------#
processLandsatBatch(x         = landsat7Dir, # directory to load data
                    outdir    = stepDir,     # output directory
                    vi        = 'nbr2',      # vegetation index
                    delete    = TRUE,        # deletes temporary files
                    mask      = 'pixel_qa',  # selects those pixles only
                    keep      = c(66, 130),  # cloudless pixels values
                    e         = extenras,
                    overwrite = TRUE,
                    mc.cores  = 1)        # overwrite files if they already
# exists in the directory
#------------------------------------------------------------------------------#
# printing time after processing landsat7 images

(land7Time2<- Sys.time())

## printing time differrence
## This let us understand how much time it took to process
## landsat 7 images
##
(landsatProcessTime<- land7Time2 - land7Time1)

####- Batch Processing Landsat8 -####
# printing time before processing landsat8 images
(land8Time1<- Sys.time())
#------------------------------------------------------------------------------#
processLandsatBatch(x         = landsat8Dir,
                    outdir    = stepDir,
                    vi        = 'nbr2',
                    delete    = TRUE,
                    mask      = 'pixel_qa',
                    keep      = c(322, 386,834,898,1346), # see link below
                    e         = extenras,
                    overwrite = TRUE,
                    mc.cores  = 1)
# you can look at the file size in nbr2 directory if they all are of same size
# https://landsat.usgs.gov/documents/lasrc_product_guide.pdf
(land8Time2<- Sys.time())
#------------------------------------------------------------------------------#
(landsatProcessTime8<- land8Time2 - land8Time1)
#------------------------------------------------------------------------------#

crop_extent <- readOGR("../NepAnalysis/Data/Bara.shp")

##Plot shapefile.
plot(crop_extent,
     main   = "Shapefile imported into R - crop extent",
     axes   = TRUE,
     border = "blue",
     add    = TRUE)

##First locate the directory where the raster are stored.
stack_list <- list.files(file.path(stepDir, 'nbr2'),
                      full.names = TRUE,
                      pattern = ".gri$")
## Looping over every stack in filelist
for(i in 1:length(stack_list)) {

  # Stack Images in list as IMG

  im  <- stack(stack_list[i])

  img<- projectRaster(im, crs = proj4string(crop_extent))
  # CRS(proj4string(test1)))
  # Calculate NDVI for each list[i]
  cropped_raster <- crop(img, extent(crop_extent))
  masked.raster  <- mask(cropped_raster,crop_extent)

  # naming outfile name
  outname <- sub(pattern     = ".grd",
                 replacement = "",
                 x           = stack_list[i])
  # export raster
  writeRaster(masked.raster, filename  = outname, overwrite = TRUE)
}
# After extentras projection

stack_list <- list.files(file.path(stepDir, 'nbr2'),
                         full.names = TRUE,
                         pattern = ".gri$")
## Looping over every stack in filelist
for(i in 1:length(stack_list)) {

  # Stack Images in list as IMG

  im  <- stack(stack_list[i])

  img<- projectRaster(im, crs = proj4string(extenras))
  # CRS(proj4string(test1)))
  # Calculate NDVI for each list[i]
  cropped_raster <- crop(img, extent(extenras))
  #masked.raster  <- mask(cropped_raster,crop_extent)

  # naming outfile name
  outname <- sub(pattern     = ".grd",
                 replacement = "",
                 x           = stack_list[i])
  # export raster
  writeRaster(cropped_raster, filename  = outname, overwrite = TRUE)
}



#------------------------------------------------------------------------------#
####- Creating Stack -####
#
### use timeStack function and Creates a time stack of Landsat layers
### USE:
### timeStack(x, pattern = NULL, orderChrono = TRUE, ...)
(nbrStackTime1<- Sys.time())
#------------------------------------------------------------------------------#
nbr2Stack <- timeStack(x         = nbr2Dir,
                       pattern   = glob2rx('*.gri'),
                       filename  = file.path(inDir, 'nbr2_stack.grd'),
                       datatype  = 'INT2S',
                     orderChrono = TRUE,
                       overwrite = TRUE)
#------------------------------------------------------------------------------#
# plotting stacked nbr images Example we are plotting #30
(nbrStackTime2<- Sys.time())

# nbr2Stack<-raster(x = file.path(inDir, 'nbr2_stack.grd'))

(nbrStackTimeDiff = nbrStackTime2 - nbrStackTime1)

#nbr2Stack<- brick(raster(file.path(inDir,'nbr2_stack.grd')))
####- Running BfastSpatial -####

#### Function (bfmSpatial) to run bfastmonitor on any kind of raster brick,
#### with parallel support
#
#### Description
#
#### Implements bfastmonitor function, from the bfast package on any kind of
#### rasterBrick object. Time information is provided as an extra object and
#### the time series can be regular or irregular.
#
#### Usage
#
# bfmSpatial(x, dates = NULL, pptype = "irregular", start, monend = NULL,
#            formula = response ~ trend + harmon, order = 3, lag = NULL,
#            slag = NULL, history = c("ROC", "BP", "all"), type = "OLS-MOSUM",
#            h = 0.25, end = 10, level = 0.05, mc.cores = 1,
#            returnLayers = c("breakpoint", "magnitude", "error"),
#            sensor = NULL, ...)
## variables after which it is commented can be changed
#  The start function specifies the start of the ##monitoring period.
#  The h value refers to the window size for the calculation of MOSUM.
#------------------------------------------------------------------------------#
#------------------------------ Pixel level -----------------------------------#
plot(nbr2Stack, 1)

# add a column for years and plot # of scenes per year
nbr2Stack$year <- as.numeric(substr(nbr2Stack$date, 1, 4))
# show the layer names
names(nbr2Stack)
s <- getSceneinfo(names(nbr2Stack))
s
# add a column for years and plot # of scenes per year
s$year <- as.numeric(substr(s$date, 1, 4))
hist(s$year, breaks=c(2004:2019), main="p170r55: Scenes per Year",
     xlab="year", ylab="# of scenes")


bfm <- bfmPixel(nbr2Stack, start=c(2009,1), interactive=TRUE,h = 0.25, plot = TRUE)

tarCell<- 18377380
bfmPixel(nbr2Stack, start=c(2007,1),
         sensor="ETM+", plot=TRUE, interactive = TRUE)
bfm$bfm
bfm1 <- bfmPixel(nbr2Stack, cell=tarCell, start=c(2009, 1),
                 formula=response~harmon, plot=TRUE)
bfm1$bfm

bfm2<- bfmPixel(nbr2Stack, cell=tarCell, start=c(2007, 1),
                formula=response~trend+harmon, order = 3, lag = NULL, plot=TRUE)
bfm2$bfm

#------------------------------------------------------------------------------#
## Prining time before running Bfast Spatial
(BfastTime1<- Sys.time())

parallel::stopCluster(clster)
# Run BfastSpatial
#nbr2Stack<-stack((file.path(inDir, 'nbr2_stack.grd')))

if (!file.exists(fn <- file.path(outDir, 'bfm_harmon_NBR2.grd'))) {
  bfm_NBR2 <- bfmSpatial(x      = crop_stack_brick,
                         pptype = 'irregular',
                         start  = c(2007,1), ## numeric data. The starting date of the
                         #monend = c(2019,11),
                         # monitoring period. Can either be given as a float
                         #(e.g.,2000.5) or a vector giving period/cycle
                         #(e.g., c(2000, 7)).
                         history = "all",
                         type    = "OLS-MOSUM",
                         h       = 0.25, ## numeric scalar from interval (0,1)
                         # specifying the bandwidth relative to the sample size
                         # in MOSUM/ME monitoring processes
                         #end = 10,
                         formula  = response ~ trend + harmon,
                         order    = 3,
                         lag = NULL,
                         #level = 0.05,
                         overwrite = TRUE,
                         filename = file.path(inDir, 'bfm_harmon_NBR2.grd')
  )
} else {
  bfm_NBR2 <- brick(fn)
}

# Printing time after running BfastSpatial
(BfastTime2<- Sys.time())

(TimeTakenBfast<- BfastTime2 - BfastTime1)

#------------------------------------------------------------------------------#
#--- post processing of Bfast spatial output ----------------------------------#
#------------------------------------------------------------------------------#
####-   Create Breakpoint, Month, and Magnitude -####
## Create objects to create months and magnitude products

months_NBR2 <- changeMonth(change_NBR2)

monthlabs   <- c("jan", "feb", "mar", "apr", "may", "jun",
                 "jul", "aug", "sep", "oct", "nov", "dec")
## applying color scheme
cols <- rainbow(12) # if required Rcolorbrewer can be used for different colors

magn_NBR2     <- raster(bfm_NBR2, 2) / 10000

magn_bkp_NBR2 <- magn_NBR2

# applying NA's where there is no change or has NAs
magn_bkp_NBR2[is.na(change_NBR2)] <- NA

# creating space for two 1 figure per column
# number of column == 2
opar <- par(mfrow = c(1, 2))

##Plot breakpoints
plot(nbr2Stack[[80]], col = grey.colors(255), legend = F)

plot(bfm_NBR2[[1]], add = TRUE)
months_NBR2 <- changeMonth(raster(bfm_NBR2,1))
monthlabs   <- c("jan", "feb", "mar", "apr", "may", "jun",
               "jul", "aug", "sep", "oct", "nov", "dec")
## applying color scheme
cols <- rainbow(12) # if required Rcolorbrewer can be used for different colors

magn_NBR2     <- raster(bfm_NBR2, 2) / 10000

magn_bkp_NBR2 <- magn_NBR2

# applying NA's where there is no change or has NAs
magn_bkp_NBR2[is.na(change_NBR2)] <- NA

# creating space for two 1 figure per column
# number of column == 2
opar <- par(mfrow = c(1, 2))

##Plot breakpoints
plot(nbr2Stack[[80]], col = grey.colors(255), legend = F)

plot(bfm_NBR2[[1]], add = TRUE)

##Plot months
plot(months_NBR2, col = cols, breaks = c(1:12), legend = FALSE)
legend("bottomright", legend = monthlabs, cex = 0.5, fill = cols,
       ncol = 2)

##Plot magnitudes
plot(magn_bkp_NBR2, main = "Magnitude of a breakpoint")
plot(magn_NBR2, main = "Magnitude: all pixels")

#------------------------------------------------------------------------------#

####- Writing Raster Files -####
### Create *.tiff files for breakpoints, months and magnitude products
##  Write raster files using write raster command
##  Write layer 1 of bfm_NBR2 layer
writeRaster(bfm_NBR2[[1]], filename = "NBR2_breaks.tif",
            format = "GTiff", outDir = outDir, overwrite = TRUE)

##  Write changeMonth layers
writeRaster(months_NBR2, filename = "NBR2_breaksmos.tif",
            format = "GTiff", outDir = outDir,  overwrite = TRUE)

## Write magnitude of breakpoints
writeRaster(magn_bkp_NBR2, filename = "NBR2_magbreaks.tif",
            format = "GTiff", outDir = outDir, overwrite = TRUE)

#------------------------------------------------------------------------------#

### The end
### Next improvement on "Spatial BFASTMonitor" :-)

change <- raster(bfm_NBR2, 1)
Cmonths <- changeMonth(change)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun",
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)


magn <- raster(bfm_NBR2, 2)
# make a version showing only breakpoing pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")
#------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------#

magn09thresh <- magn_bkp_NBR2
magn09thresh[magn_bkp_NBR2 > -0.05] <- NA
# compare
op <- par(mfrow=c(1, 2))
plot(magn_bkp_NBR2, main="magnitude")
plot(magn09thresh, main="magnitude < -0.05")

#------------------------------------------------------------------------------#
magn09_sieve <- areaSieve(magn09thresh, thresh=1800)
install.packages("igraph")
library(igraph)
magn09_areasieve <- areaSieve(magn09thresh)
magn09_as_rook <- areaSieve(magn09thresh, directions=4)
# compare all magn rasters
op <- par(mfrow=c(2, 2))
plot(magn09thresh, main="magnitude")
plot(magn09_sieve, main="pixel sieve")
plot(magn09_areasieve, main="0.5ha sieve")
plot(magn09_as_rook, main="0.5ha sieve, rook's case")
