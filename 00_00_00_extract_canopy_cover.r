############# --------------- TREE CROWN based on lidar
################## ---------- Lubbock County, Texas, USA
####- 1. Global Option -####
# -----------------------------------------------------------------------------#
options(timeout   = 15,     prompt    = "< ",
        continue  = "...",  digits    = 5,
        width     =  80,    papersize = "letter" )

####- 2. Clear workspace -####
# -----------------------------------------------------------------------------#
rm(list = ls())

####- 3. Load Package -####
# -----------------------------------------------------------------------------#
library(sf)         # sf
library(imager)     # for image data
library(lidR)       # lpc data 
library(RMCC)       # require for mcc algorithm
library(terra)      # shade/raster data


####- 4. Read Data -####
# -----------------------------------------------------------------------------#
drive = "K:/"
tif_path <- paste0(drive, "Mukti_TTU_CC/LBB_Project/downloadLidar/LB_2018_Out2/")
out_path <- paste0(drive, "Mukti_TTU_CC/LBB_Project/downloadLidar/")
rast_vect <- list.files(tif_path,pattern = "*.tif$",full.names = FALSE)
#------------------------------------------------------------------------------#
#### - 4. Source Functions  - ####
#### - 4.1 Maxima Detection - ####
local_maxima_detection <- function(chm, max_width = 12, id_field = "tree_id",
                                   filter_type = "median", dmin = 1, ddrop = 0.05,
                                   min_height = 2.0, drop_zero = TRUE) {
  # Convert raster to cimg object if necessary
  if (inherits(chm, "SpatRaster")) {
    #chm_gs <- rast_image_cimg(chm)
    chm_gs <- as.cimg(terra::as.matrix(chm, wide = TRUE), 
                      na.rm = TRUE, na.val = 0)
    chm_res <- terra::res(chm)[1]
    max_width <- max_width / chm_res
  } else {
    chm_gs <- chm
  }
  
  
  #Apply condition to set values below zero to zero
  if (drop_zero) {
    chm_gs[chm_gs < 0] <- 0
  }
  if (!is.null(min_height)) {
    chm_gs[chm_gs < min_height] <- 0
  }
  
  
  # # Apply median or other filters
  # if (tolower(filter_type) == "median") {
  #   chm_gs <- imager::medianblur(chm_gs, n = 3)
  # } else if (tolower(filter_type) == "mean") {
  #   chm_gs <- imager::imfilter(chm_gs, 
  #                              filter = matrix(1, nrow = 3, ncol = 3) / 9)
  # }
  # 
  max_radius <- max_width %/% 2
  
  # Create a square structuring element for dilation
  strel <- imager::imfill(3, 3, val = 1)
  
  # Initialize result
  chm_dilate <- imager::dilate(chm_gs, strel)
  maxi <- imager::as.cimg(chm_gs == chm_dilate)
  
  # Loop on increasing radius
  for (i in 2:max_radius) {
    chm_dilate <- imager::dilate(chm_dilate, strel)
    maxi <- imager::parmax(list(maxi, 
                                imager::as.cimg(chm_gs == chm_dilate) * i), na.rm = TRUE)
  }
  
  ### initial chm for maxima detection = chm_gs
  max_copy<- maxi
  
  max_copy[chm_gs< min_height]<- 0
  
  max_copy[max_copy<(dmin+chm_gs *ddrop)]<0
  
  # Convert cimg objects to raster if necessary
  if (inherits(chm, "SpatRaster")) {
    maxi <- cimg_raster(max_copy, chm)
  }
  plot(maxi)
  #-------------------------- convert to polygon ###############
  # Convert to points
  local_max_pts <- sf::st_as_sf(as.points(maxi, na.rm = TRUE, values = TRUE))
  names(local_max_pts)[1] <- 'height'
  #local_max_pts<- local_max_pts %>% dplyr::filter(height>min_height)
  local_max_pts[[id_field]] <- 1:nrow(local_max_pts)
  local_max_pts <- local_max_pts[, c(id_field, "height", "geometry")]
  
  # Return output
  return(local_max_pts)
  #return(maxi)
}

#------------------------------ Local Maxima End ------------------------------#

#### -4.2 watershed segmentation -####
#------------------------------------------------------------------------------#
marker_watershed <- function(tree_tops, chm,ndvi = ndvi_resp, min_height = 2,
                             format = "raster", tree_id_field = "treeID",
                             min_slope = NULL,
                             slope_neighbors = 8,
                             slope_unit = "degrees", ndvi_resp = ndvi_resp) {
  # Apply watershed segmentation
  chm_mask <- chm < min_height
  chm_mask[is.na(chm)] <- TRUE
  # Convert tree_tops to a raster
  tree_tops_ras <- terra::rasterize(tree_tops, chm, field = tree_id_field, 
                                    background = 0)
  
  # Convert data to 'cimg' format
  chm_img   <- imager::as.cimg(terra::as.matrix(chm, wide = TRUE))
  ttops_img <- imager::as.cimg(terra::as.matrix(tree_tops_ras, wide = TRUE))
  
  # Apply watershed function
  ws_img <- imager::watershed(ttops_img, chm_img)
  ws_ras <- terra::rast(as.matrix(ws_img), extent = terra::rast(chm), crs = terra::crs(chm))
  ws_ras[chm_mask] <- NA
  
  if (toupper(format) %in% c("POLYGONS", "POLYGON", "POLY")) {
    # Convert raster to polygons
    polys <- sf::st_as_sf(terra::as.polygons(ws_ras, dissolve = TRUE))
    polys <- sf::st_cast(polys, "MULTIPOLYGON")
    polys <- sf::st_cast(polys, "POLYGON", warn = FALSE)
    
    if (nrow(polys) == 0) {
      stop("No segments created")
    }
    polys <- sf::st_as_sf(fillHoles(terra::vect(polys)))
    #----------------------------------------  # merge
    # polys_union <- st_combine(polys) %>%st_union(.,by_feature = FALSE)
    # polys <- st_cast(polys_union, "POLYGON") %>% st_as_sf()
    # polys <- polys %>% dplyr::rename(geometry = x)
    # names(polygon_sf)[names(polygon_sf) == "x"] <- "geometry"
    st_geometry(polys) <- "geometry"
    names(polys)[1] <- tree_id_field
    row.names(polys) <- 1:nrow(polys)
    # polys$tree_id_field <- 1:nrow(polys)
    # Calculate slope raster
    slope_raster <- terra::terrain(chm,v = "slope", 
                                   neighbors = slope_neighbors,
                                   unit = slope_unit)
    # Extract the mean slope within each polygon
    mean_slopes <- terra::extract(slope_raster, polys, fun = mean, na.rm = TRUE)
    mean_height <- terra::extract(chm, polys, fun = mean, na.rm = TRUE)
    mean_ndvi   <- terra::extract(ndvi_resp, polys, fun = mean, na.rm = TRUE)
    min_ndvi    <- terra::extract(ndvi_resp, polys, fun = min, na.rm = TRUE)
    min_ht      <- terra::extract(chm, polys, fun = min, na.rm = TRUE)
    names(mean_ndvi)[2]  <- "mean_ndvi"
    names(mean_height)[2]<- "mean_ht"
    names(min_ndvi)[2]   <- "min_ndvi"
    names(min_ht)[2]<- "min_ht"
    attbr = mean_slopes %>% dplyr::left_join(mean_height, by = c("ID" = "ID")) %>%
            dplyr::left_join(min_ht, by = c("ID" = "ID")) %>% 
            dplyr::left_join(mean_ndvi, by = c("ID" = "ID")) %>% 
            dplyr::left_join(min_ndvi, by = c("ID" = "ID"))
     
    # Add a "slope" column to the polygons
    #polys$slope <- mean_slopes
    polys<- cbind(polys, attbr[,2:6])
    # Filter polygons based on the mean slope threshold
    if (!is.null(min_slope)) {
      polys <- polys[polys$slope >= min_slope, ]
      if (nrow(polys) == 0) {
        stop("No polygons meet the slope threshold")
      }
    }
    
    return(polys)
  } else {
    return(ws_ras)
  }
}

####------------------------- Watershed Segmentation  End ------------------####
#------------------------------------------------------------------------------#

# if(!dir.exists(paste0(out_path, "LB_2011_shp"))){
#   dir.create(paste0(out_path,   "LB_2011_shp"))
# }

out_folder = "LB_2018_shp_2"

if(!dir.exists(paste0(out_path, out_folder))){
  dir.create(paste0(out_path,   out_folder))
}

ndvi_path  = paste0(drive,"Mukti_TTU_CC/chm_2018/ndvi_2019.img")
##### - 5. Create Function - ####

rast_image_cimg <- function(raster, NA_replace = 0, maxpixels = 1e+10) {
  # replace NA values
  raster[is.na(raster)] <- NA_replace
  # convert to cimg object
  if (ncol(raster) * nrow(raster) > maxpixels) {
    warning("Too many raster pixels: conversion to cimg partial; try with higher maxpixels arguments")
  }
  imager::as.cimg(matrix(raster, ncol = nrow(raster)), maxpixels = maxpixels)
}
cimg_raster <- function(cimg, reference = NULL) {
  # Convert to SpatRaster
  chm <- terra::rast(t(as.matrix(cimg)))
  
  # If a reference raster is provided
  if (!is.null(reference)) {
    # Set extent and CRS from the reference
    terra::ext(chm) <- terra::ext(reference)
    terra::crs(chm) <- terra::crs(reference)
  }
  chm[chm ==0]<- NA
  return(chm)
}
#------------------------------------------------------------------------------#
process_all_data <- function(rast_vect = rast_vect[1], min_height = 2,
                             drop_zero = TRUE, out_folder = out_folder,
                             ndvi_path = ndvi_path)
{ # main function
  #------------------------------------------------------------------------------#  
  ndvi <- terra::rast(ndvi_path)
  # ndvi_resp <- terra::resample(ndvi,y = chm, method = "bilinear")
  # 
  # # Create a mask based on the conditions
  # mask <- (chm > 2) & (ndvi_resp > 0.30)
  # 
  # # Apply the mask to the chm raster
  # #filtered_chm <- mask(chm, mask)
  # filtered_chm <- chm*terra::values(mask, ifelse = NA)
  
  for(file in rast_vect){ # start for loop
    ind = which(rast_vect == file)
    cat(sprintf("processing file: %s:remaining file are:%d", file,length(rast_vect)-ind),"\n")
    img <- terra::rast(paste0(tif_path,file))
    ## ----------------------------------------------------------------------##
    ndvi_resp = terra::resample(ndvi,y = img, method = "bilinear")
    # Create a mask based on the conditions
    mask <- (img > 2) & (ndvi_resp > 0.30)
    # Apply the mask to the chm raster
    #filtered_chm <- mask(chm, mask)
    filtered_chm <- img*terra::values(mask, ifelse = NA)
    ## ----------------------------------------------------------------------##  
    ttops = local_maxima_detection(chm = filtered_chm,max_width = 12,
                                   id_field = "tree_id",
                                   filter_type = "median",min_height = 2.0,
                                   drop_zero = TRUE, dmin = 1, ddrop = 0.05)
    
    cat(sprintf("Tree tops is calculated for: %s\nnrow: = %d", file, nrow(ttops)),"\n")
    
    poly = marker_watershed(tree_tops = ttops,
                            chm = filtered_chm,
                            min_height = 2,format = "POLY",
                            tree_id_field = "tree_id",
                            min_slope = NULL,slope_neighbors = 8,
                            slope_unit = "degrees", ndvi_resp = ndvi_resp)
    name = sub("\\.[^.]+$", "", file) 
    cat(sprintf("Writing file: %s", file))  
    
    if (nrow(poly) == 0) {
      next}
    sf::st_write(poly,paste0(out_path,out_folder,"/",name,"_",ind,".shp"),driver = "ESRI Shapefile")
    
  } # close for loop
}# close main function

process_all_data(rast_vect = rast_vect[1],min_height = 2,drop_zero = TRUE,out_folder = out_folder,ndvi_path = ndvi_path)