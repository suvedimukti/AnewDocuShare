# ============================================================================ #
# Site index modelling for GA :- Loblolly Pine
#------------------------------------------------------------------------------#
#### - 1. set global option - ####
#------------------------------------------------------------------------------#

options( continue = "...", digits = 5, width = 80, timeout = 30, prompt = "> ")

#-------------- remove everything from environment

rm(list = ls())

####- 2. Load library - ####
#------------------------------------------------------------------------------#
library(nls2)
library(tidyverse)
library(broom)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(patchwork)
library(nlme)
library(lme4)
library(robustbase)
library(ModelMetrics)
library(rsample)
library(magrittr)
library(data.table)
library(NISTnls)
library(gslnls)
library(glue)
library(sf)
library(meteo)
library(sp)
library(sftime)
library(terra)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)
library(ranger)

#------------------------------------------------------------------------------#
#### -1. Load Data - ####
#------------------------------------------------------------------------------#
# 1. Read SI data
alldat_loc <- readRDS("../02_01_Output_Prj03/00_01_Combined_Outputs/combined_RSPPPLPYP.rds")
####- based on normalization -####

#------------------------------------------------------------------------------#
si_data <- alldat_loc %>% st_as_sf(coords = c("rlong", "rlat"),
                                   crs = "EPSG:4269") %>% st_transform(., crs = "EPSG:5070")

# 2. Read soil with grid data

ga_grid_clim_soil <- readRDS("../02_01_Output_Prj03/01_01_SoilStacRaster/ga_grid_clim_soil.rds")


si_dat_sel <- si_data %>% dplyr::select(all_of(c("si_b0", "rsi")))

si_dat_extract <- st_intersection(x = ga_grid_clim_soil,y = si_dat_sel)

#### - Variogram - ####
sf_dat1 <- si_dat_extract %>% st_drop_geometry(si_dat_extract) %>%
            mutate(x = sf::st_coordinates(si_dat_extract)[,1],
                   y = sf::st_coordinates(si_dat_extract)[,2])

fit_multiple_variogram<- function (data, coords = c("X", "Y"), length = 99)

mul_var <- stdcab::fit_variogram(data = si_dat_extract,response = "rsi",coords = NULL)

# multi_var
vif_var = c("PPT_wt", "DD_0_sp", "DD_0_at", "DD_18_sm", "NFFD_wt",
            "NFFD_at", "PAS_wt", "CMD_sm", "CMD_at", "RH_wt", "RH_sp", "RH_sm",
            "RH_at", "EXT", "RH", "awc_r", "claytotal_r", "sandtotal_r",
            "silttotal_r", "dbovendry_r", "om_r", "ksat_r", "kffact", "pi_r",
            "MAT", "MAP", "MSP", "SHM", "CMI", "Tave_sp",
            "Tave_sm", "Tave_wt"
)

mul_pred <- stdcab::fit_multiple_variogram(data = si_dat_extract[, vif_var],coords =NULL,length = 100 )

plt_mulvar <- stdcab::plot_variogram(object = mul_var,length = 100, show_range = TRUE)

(plt_predvar<- stdcab::plot_multiple_variogram(object = mul_pred,plot_type = "bar",show_range = TRUE))

plt_predvar$plot+
  # scale_y_continuous(trans = 'log10',
  #                                     breaks = scales::trans_breaks('log10', function(x) 10^x),
  #                                     labels = function(x)round(x,3), limits = c(-0.001, 10000))+
  coord_flip()+scale_y_log10()+guides(color = guide_legend(override.aes = aes(label = "")))
cluster_data = stdcab::spatial_quadgrid_sample(data = si_dat_extract,
                                               cellsize = c(mean(plt_mulvar$data$range, na.rm = TRUE), mean(plt_mulvar$data$range, na.rm = TRUE)),
                                               show_grid = TRUE,
                                               rotation_angle = 0,
                                               fold_selection = "random",
                                               k = 10)

block<- cluster_data$blocks
block$obs<- NULL
block$gname <- NULL


#### - 5. RFSI - ####
#------------------------------------------------------------------------------#

si_dat <- st_intersection(x = block,y = si_dat_extract)

st_geometry(si_dat)<- "geometry"

sidf = cbind(st_drop_geometry(si_dat), data.frame(st_coordinates(si_dat)))

cluster_data = stdcab::spatial_cluster_sample(data = sidf,coords =c("X", "Y"),
                                              v = 10,spatial = FALSE,
                                              clust_method = "hclust",
                                              dist_clust = "ward.D")
# library(rsample)

df_train_test <- rsample::rsample2caret(cluster_data, data = c("analysis", "assessment"))


si_dat<- si_dat %>% select_if(~ !is.numeric(.) || sum(.) != 0)
# preparing data
dep_var =  c("Tmax_wt",
             "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmin_wt", "Tmin_sp", "Tmin_sm",
             "Tmin_at", "Tave_wt", "Tave_sp", "Tave_sm", "Tave_at", "PPT_wt",
             "PPT_sp", "PPT_sm", "PPT_at", "DD_0_wt", "DD_0_sp", "DD_0_at",
             "DD5_wt", "DD5_sp", "DD5_sm", "DD5_at", "DD_18_wt", "DD_18_sp",
             "DD_18_sm", "DD_18_at", "DD18_wt", "DD18_sp", "DD18_sm", "DD18_at",
             "NFFD_wt", "NFFD_sp", "NFFD_sm", "NFFD_at", "PAS_wt", "Eref_wt",
             "Eref_sp", "Eref_sm", "Eref_at", "CMD_sp", "CMD_sm", "CMD_at",
             "RH_wt", "RH_sp", "RH_sm", "RH_at", "CMI_wt", "CMI_sp", "CMI_sm",
             "CMI_at", "MAT", "MWMT", "MCMT", "TD", "MAP", "MSP", "AHM", "SHM",
             "DD_0", "DD5", "DD_18", "DD18", "NFFD", "bFFP", "eFFP", "FFP",
             "PAS", "EMT", "EXT", "Eref", "CMD", "RH", "CMI", "DD1040", "awc_r",
             "claytotal_r", "sandtotal_r", "silttotal_r", "dbovendry_r", "om_r",
             "ksat_r", "kffact", "pi_r")


equn_form <- as.formula(paste0("rsi ~ ", paste0(dep_var, collapse = " + ")))

class(si_dat)

rfsi_model <- meteo::rfsi(formula = equn_form,
                          data = si_dat,
                          #zero.tol = 0,
                          n.obs = 5, # number of nearest observations
                          # s.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
                          # p.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
                          cpus = detectCores()-16,
                          progress = TRUE,
                          # ranger parameters
                          importance = "impurity",
                          seed = 42,
                          num.trees = 2000,
                          mtry = 25,
                          splitrule = "variance",
                          min.node.size = 10,
                          sample.fraction = 0.95,
                          quantreg = TRUE) # quantile regression model

rfsi_model
var_imp <- rfsi_model$variable.importance
var_imp <- importance(rfsi_model,  scale = TRUE)
var_imp <- data.frame(variables = names(var_imp), importance = var_imp)

# readxl::read_excel("D:/FPA_LIDAR_2023/10_Project_03/08_ClimateData/ClimateNA_v742/Variables.xlsx",sheet = "Climate Data",col_names = TRUE)

var_imp <- var_imp %>% dplyr::filter(!str_detect(variables, c("obs"))) %>%
  dplyr::filter(!str_detect(variables, c("dist")))

## Create a plot of variable importance
var_imp %>%

  ## Sort the data by importance
  arrange(importance) %>% dplyr::slice_max(order_by = importance, n = 30) %>%

  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) +

  ## Plot the bar graph
  geom_bar(stat='identity') +

  ## Flip the graph to make a horizontal bar plot
  coord_flip() +

  ## Add x-axis label
  xlab('Variables\n') +

  ## Add a title
  labs(y = "\nImportance") +

  ## Some layout for the plot
  # theme_bw() +
  # theme(axis.text = element_text(size = 10),
  #       axis.title = element_text(size = 15),
  #       plot.title = element_text(size = 20),
  # )
  theme_bw()+
  theme(
    axis.text.y.left    = element_text(face = "bold", size = 12),
    axis.title.y.left   = element_text(face = "bold", size = 14),
    axis.text.x.bottom  = element_text(face = "bold", size = 12),
    axis.title.x.bottom = element_text(face = "bold", size = 14),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.ticks          = element_line(linewidth = 1),
    axis.ticks.length.x = unit(-0.1, "cm"),
    panel.grid.major    = element_line(linetype = "longdash", color = "grey95"),
    panel.grid.minor    = element_line(linetype = "dotted", color = "grey90")
  )


# Set desired resolution (1km x 1km in this case)
ga_grid_clim_soil %<>% mutate_at(vars(awc_r, claytotal_r, sandtotal_r,
                                      silttotal_r, dbovendry_r,
                                      om_r, ksat_r, kffact, pi_r),
                                 ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

# ga_grid_clim_soil
# res <- 1000
# # Create a raster with the same extent and resolution as the polygons
# template_raster <- rast(ext(ga_grid_clim_soil), resolution = res)
# # Convert each variable in the polygons to a raster layer, using the template raster
# raster_stack <- rast(lapply(1:89, function(i) rasterize(ga_grid_clim_soil,
#                                   template_raster, field = names(ga_grid_clim_soil)[-c(1:4)][i])))
si_dat_vif <- si_dat %>%
  dplyr::select(any_of(c("PPT_wt", "DD_0_sp", "DD_0_at", "DD_18_sm", "NFFD_wt",
                         "NFFD_at", "PAS_wt", "CMD_sm", "CMD_at", "RH_wt", "RH_sp", "RH_sm",
                         "RH_at", "EXT", "RH", "awc_r", "claytotal_r", "sandtotal_r",
                         "silttotal_r", "dbovendry_r", "om_r", "ksat_r", "kffact", "pi_r",
                         "rsi", "MAT", "MAP", "MSP", "SHM", "CMI", "Tave_wt", "Tave_sp", "Tave_sm", "Tave_wt")))

equn_rfsi = as.formula(paste0("rsi ~ ",
                              paste0(c("PPT_wt", "DD_0_sp", "DD_0_at", "DD_18_sm", "NFFD_wt",
                                       "NFFD_at", "PAS_wt", "CMD_sm", "CMD_at", "RH_wt", "RH_sp", "RH_sm",
                                       "RH_at", "EXT", "RH", "awc_r", "claytotal_r", "sandtotal_r",
                                       "silttotal_r", "dbovendry_r", "om_r", "ksat_r", "kffact", "pi_r",
                                       "MAT", "MAP", "MSP", "SHM", "CMI", "Tave_wt", "Tave_sp",
                                       "Tave_sm", "Tave_wt"),
                                     collapse = "+")))


##----------------------------------------------------------------------------##
# si_dat_fil = cbind(st_drop_geometry(si_dat_vif),
#                    data.frame(st_coordinates(si_dat_vif)))
vif_var = c("PPT_wt", "DD_0_sp", "DD_0_at", "DD_18_sm", "NFFD_wt",
            "NFFD_at", "PAS_wt", "CMD_sm", "CMD_at", "RH_wt", "RH_sp", "RH_sm",
            "RH_at", "EXT", "RH", "awc_r", "claytotal_r", "sandtotal_r",
            "silttotal_r", "dbovendry_r", "om_r", "ksat_r", "kffact", "pi_r",
            "MAT", "MAP", "MSP", "SHM", "CMI", "Tave_sp",
            "Tave_sm", "Tave_wt"
)
tune_for = as.formula(paste0("rsi ~ ", paste0(vif_var, collapse = "+")))
si_dat_fil1 = si_dat[, c(vif_var,"rsi", "sgname", "id", "tile")]
si_dat_fil <-  si_dat %>%select_if(~ !is.numeric(.) || sum(.) != 0)
si_dat_fil1 = cbind(st_drop_geometry(si_dat_fil), data.frame(st_coordinates(si_dat_fil)))

cluster_data = stdcab::spatial_cluster_sample(data = si_dat_fil1,coords =c("X", "Y"),
                                              v = 10,spatial = FALSE,
                                              clust_method = "hclust",
                                              dist_clust = "ward.D")

cluster_data1 = stdcab::spatial_quadgrid_sample(data = si_dat_fil,
                                                cellsize = c(mean(plt_mulvar$data$range, na.rm = TRUE), mean(plt_mulvar$data$range, na.rm = TRUE)),
                                                show_grid = TRUE,fold_selection = "random",k = 10
)
# library(rsample)
dep_var = c("sgname", "id", "tile", "Latitude", "Longitude", "Tmax_wt",
            "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmin_wt", "Tmin_sp", "Tmin_sm",
            "Tmin_at", "Tave_wt", "Tave_sp", "Tave_sm", "Tave_at", "PPT_wt",
            "PPT_sp", "PPT_sm", "PPT_at", "DD_0_wt", "DD_0_sp", "DD_0_at",
            "DD5_wt", "DD5_sp", "DD5_sm", "DD5_at", "DD_18_wt", "DD_18_sp",
            "DD_18_sm", "DD_18_at", "DD18_wt", "DD18_sp", "DD18_sm", "DD18_at",
            "NFFD_wt", "NFFD_sp", "NFFD_sm", "NFFD_at", "PAS_wt", "Eref_wt",
            "Eref_sp", "Eref_sm", "Eref_at", "CMD_sp", "CMD_sm", "CMD_at",
            "RH_wt", "RH_sp", "RH_sm", "RH_at", "CMI_wt", "CMI_sp", "CMI_sm",
            "CMI_at", "MAT", "MWMT", "MCMT", "TD", "MAP", "MSP", "AHM", "SHM",
            "DD_0", "DD5", "DD_18", "DD18", "NFFD", "bFFP", "eFFP", "FFP",
            "PAS", "EMT", "EXT", "Eref", "CMD", "RH", "CMI", "DD1040", "awc_r",
            "claytotal_r", "sandtotal_r", "silttotal_r", "dbovendry_r", "om_r",
            "ksat_r", "kffact", "pi_r")
dev_var = c("Tmax_wt",
            "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmin_wt", "Tmin_sp", "Tmin_sm",
            "Tmin_at", "Tave_wt", "Tave_sp", "Tave_sm", "Tave_at", "PPT_wt",
            "PPT_sp", "PPT_sm", "PPT_at", "DD_0_wt", "DD_0_sp", "DD_0_at",
            "DD5_wt", "DD5_sp", "DD5_sm", "DD5_at", "DD_18_wt", "DD_18_sp",
            "DD_18_sm", "DD_18_at", "DD18_wt", "DD18_sp", "DD18_sm", "DD18_at",
            "NFFD_wt", "NFFD_sp", "NFFD_sm", "NFFD_at", "PAS_wt", "Eref_wt",
            "Eref_sp", "Eref_sm", "Eref_at", "CMD_sp", "CMD_sm", "CMD_at",
            "RH_wt", "RH_sp", "RH_sm", "RH_at", "CMI_wt", "CMI_sp", "CMI_sm",
            "CMI_at", "MAT", "MWMT", "MCMT", "TD", "MAP", "MSP", "AHM", "SHM",
            "DD_0", "DD5", "DD_18", "DD18", "NFFD", "bFFP", "eFFP", "FFP",
            "PAS", "EMT", "EXT", "Eref", "CMD", "RH", "CMI", "DD1040", "awc_r",
            "claytotal_r", "sandtotal_r", "silttotal_r", "dbovendry_r", "om_r",
            "ksat_r", "kffact", "pi_r")
dev_var<- c("awc_r", "claytotal_r", "CMD_at", "CMD_sm", "dbovendry_r",
            "DD_0_at", "DD_0_sp", "EXT", "kffact", "ksat_r", "NFFD_sm", "om_r",
            "pi_r", "PPT_wt", "RH", "RH_at", "RH_sm", "RH_sp", "RH_wt",
            "silttotal_r","CMI")

equn_f <- as.formula(paste0("rsi ~ ", paste0(dep_var, collapse = "+")))
equn_f <- as.formula(paste0("rsi ~ ", paste0(dev_var, collapse = "+")))
df_train_test <- rsample::rsample2caret(cluster_data, data = c("analysis", "assessment"))

# Define the new class names
new_classes <- c("spatial_clustering_cv", "rset", "tbl_df", "tbl", "data.frame")

# Assign the new class names to your object
class(cluster_data1$splits) <- new_classes

# Verify the changes
class(cluster_data1$splits)

df_train_test1<- rsample::rsample2caret(cluster_data1$splits,data = c("analysis", "assessment"))

#--------------- fit statistics ------------------------------#


estimate_fit_stats <- function(fold_obs, fold_pred){
  # if(length(dtype)>1){
  #   stop('data type should of one of "train" or "test" ')
  # }
  # #----------------------------------------------------------------------------#
  # if (any(grepl(pattern = "log", x = summary(model)$formula) == TRUE) && dtype == "test") {
  #   pr <- exp(as.vector(predict(model, df)))
  # } else if (any(grepl(pattern = "log", x = summary(model)$formula) == FALSE) && dtype == "test") {
  #   pr <- as.vector(predict(model, df))
  # } else if (any(grepl(pattern = "log", x = summary(model)$formula) == TRUE) && dtype == "train") {
  #   pr <- exp(as.vector(fitted(model, df)))
  # } else {
  #   pr <- as.vector(fitted(model, df))
  # }
  #-------------------------------------Mean Bias -----------------------------#
  obs <- fold_obs
  pr <- fold_pred
  n  <- length(fold_obs)
  mb <- sum(pr-obs)/n
  #return(mb)
  #-------------------------------------    MAE -------------------------------#
  # mae
  absval = abs(obs-pr)
  mae = sum(absval)/n
  #-----------------------------   RMSE ---------------------------------------#
  paras  <- length(coef(model))
  sumsq  <- sum((pr-obs)^2)
  rmse     <- sqrt(sumsq/(n-paras))
  #---------------------------------   RSQ adj --------------------------------#
  rss  <-sum((pr-obs)^2)
  tss  <- sum((obs - mean(pr))^2)
  par1 <- length(coef(model))-1
  rs = cor(obs,pr)^2
  r_sq = rs
  # r_sq <- 1- (rss/tss)
  r_ad <- 1-(((1-r_sq)*(n-1))/(n-par1))
  #------------------------------------ info criteria
  AIC <- n * log(rss / n) + 2 * paras
  BIC <- n * log(rss / n) + paras * log(n)

  ret_df <- round(data.frame("mean_bias" = mb, "MAE" = mae, "RMSE" = rmse,
                             "adj_Rsq" = r_ad,"AIC" = AIC, "BIC" = BIC),3)
  ret_df$pars <- paste0("par_comb_", h)
  ret_df$fold <- paste0("fold_",f)
  ret_df<- subset(ret_df, select = c("pars" ,"fold", "mean_bias", "MAE", "RMSE", "adj_Rsq", "AIC", "BIC"))
  return(ret_df)
}

#-----------------------------------------------------------------------------#
# it was df_train_test before
# indices = df_train_test
indices = df_train_test1
# making tgrid
n.obs <- 1:10
min.node.size <- 2:30
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 500 #c(200, 500, 700, 1000, 1500) # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

hp <- expand.grid(min.node.size = min.node.size,
                  mtry = mtry, sf = sample.fraction)
set.seed(1318)
hp <- hp[sample(nrow(hp), 50),]
hp <- hp[order(hp$mtry),]
rmse_hp <- list()#rep(NA, nrow(hp))

for (h in 1:nrow(hp)) {
  comb <- hp[h, ]
  print(paste("combination: ", h, sep=""))
  print(comb)
  fold_obs <- c()
  fold_pred <- c()

  for (f in 1:length(indices$index)) {
    print(paste("fold: ", f, sep=""))
    dev_df1 <- si_dat_fil1[indices$index[[f]], ]
    val_df1 <- si_dat_fil1[indices$indexOut[[f]],  ] #c("rsi",dep_var)
    # comb <- rfsi_tuned$tuned.parameters
    model <- ranger(tune_for, data = si_dat_fil1, importance = "none", seed = 42,
                    num.trees = ntree, mtry = comb$mtry,
                    splitrule = "variance",
                    min.node.size = comb$min.node.size,
                    sample.fraction = comb$sf,
                    oob.error = FALSE)
    # model = rfsi_model2
    fold_obs  <- c(fold_obs, val_df1$rsi)
    fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
    #df_stat   <- data.frame( obs = fold_obs, pred = fold_pred)
    #rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
    fstat = estimate_fit_stats(fold_obs = fold_obs, fold_pred = fold_pred)
  }
  # rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  rmse_hp[[h]] <- fstat
  print(paste("rmse: ", rmse_hp[h], sep=""))

}
ga_shp = sf::st_read("../04_GAshpfile/ga_inlad_shp.shp")

newdata <- terra::rast(si_dat_extract)
class(newdata)
df<- st_centroid(ga_grid_clim_soil)
#v <- c(awc_r, claytotal_r, sandtotal_r, silttotal_r, dbovendry_r, om_r, ksat_r, kffact, pi_r)
df %<>% mutate_at(vars(awc_r, claytotal_r, sandtotal_r, silttotal_r, dbovendry_r, om_r, ksat_r, kffact, pi_r),~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

# df1<- df
# df1<- as_Spatial(df1)
# sp::gridded(df1)<- TRUE
# testspat <- terra::rast(df1)

ranp = predict(model,data = st_drop_geometry(df))
df$rsi = ranp$predictions

#tp = terra::rasterizeGeom(df,res =1000,unit = "m")
df1<- df
df1<- as_Spatial(df1)
sp::gridded(df1)<- TRUE
tp <- terra::rast(df1)

tp1<- tp$rsi
#----------------------------------------#
rfp_clip = terra::crop(x = tp1,y = terra::vect(ga_shp), mask = TRUE)

plot(rfp_clip)

#-----------------------------------------------------------------------------#


df_hp = do.call(rbind, rmse_hp)

which(df_hp[, "adj_Rsq"] == max(df_hp[, "adj_Rsq"]))

# best is hp[49, ]
comb =  hp[which(df_hp[, "adj_Rsq"] == max(df_hp[, "adj_Rsq"]))[1], ]
# now estimate for 10 fold_obs

estimate_fit_stats_folds <- function(fold_obs, fold_pred){
  # if(length(dtype)>1){
  #   stop('data type should of one of "train" or "test" ')
  # }
  # #----------------------------------------------------------------------------#
  # if (any(grepl(pattern = "log", x = summary(model)$formula) == TRUE) && dtype == "test") {
  #   pr <- exp(as.vector(predict(model, df)))
  # } else if (any(grepl(pattern = "log", x = summary(model)$formula) == FALSE) && dtype == "test") {
  #   pr <- as.vector(predict(model, df))
  # } else if (any(grepl(pattern = "log", x = summary(model)$formula) == TRUE) && dtype == "train") {
  #   pr <- exp(as.vector(fitted(model, df)))
  # } else {
  #   pr <- as.vector(fitted(model, df))
  # }
  #-------------------------------------Mean Bias -----------------------------#
  obs <- fold_obs
  pr <- fold_pred
  n  <- length(fold_obs)
  mb <- sum(pr-obs)/n
  #return(mb)
  #-------------------------------------    MAE -------------------------------#
  # mae
  absval = abs(obs-pr)
  mae = sum(absval)/n
  #-----------------------------   RMSE ---------------------------------------#
  paras  <- length(coef(model))
  sumsq  <- sum((pr-obs)^2)
  rmse     <- sqrt(sumsq/(n-paras))
  #---------------------------------   RSQ adj --------------------------------#
  rss  <-sum((pr-obs)^2)
  tss  <- sum((obs - mean(pr))^2)
  par1 <- length(coef(model))-1
  r_sq = cor(obs, pr)^2
  #r_sq <- 1- (rss/tss)
  r_ad <- 1-(((1-r_sq)*(n-1))/(n-par1))
  #------------------------------------ info criteria
  AIC <- n * log(rss / n) + 2 * paras
  BIC <- n * log(rss / n) + paras * log(n)
  CCC <- DescTools::CCC(x = obs, y = pr, ci = "z-transform",
                        conf.level = 0.95, na.rm=TRUE)$rho.c[1,1]

  ret_df <- round(data.frame("mean_bias" = mb, "MAE" = mae, "RMSE" = rmse,
                             "adj_Rsq" = r_ad, "CCC" = CCC, "AIC" = AIC, "BIC" = BIC),3)
  # ret_df$pars <- paste0("par_comb_", h)
  ret_df$fold <- paste0("fold_",f)
  ret_df<- subset(ret_df, select = c( "fold", "mean_bias", "MAE", "RMSE", "adj_Rsq", "CCC", "AIC", "BIC"))
  return(ret_df)
}
rmse_hp1<- rmse_hp
rmse_hp <- list()
# model   <- ranger(equn_rfsi, data = si_dat_fil1, importance = "none", seed = 42,
#                   num.trees = ntree, mtry = comb$mtry,
#                   splitrule = "variance",
#                   min.node.size = comb$min.node.size,
#                   sample.fraction = comb$sf,
#                   oob.error = FALSE)

model <- ranger(equn_f, data = si_dat_fil1, importance = "none", seed = 42,
                num.trees = ntree, mtry = max(comb$mtry),
                splitrule = "variance",
                min.node.size = max(comb$min.node.size),
                sample.fraction = max(comb$sf),
                oob.error = FALSE)

# model <- meteo::rfsi(equn_f, data = si_dat_fil, importance = "none", seed = 42,
#                      num.trees = ntree, mtry = comb$mtry,
#                      splitrule = "variance",
#                      min.node.size = comb$min.node.size,
#                      sample.fraction = comb$sf,
#                      oob.error = FALSE
# )

for (f in 1:length(indices$index)) {
  print(paste("fold: ", f, sep=""))
  dev_df1 <- si_dat_fil1[indices$index[[f]], ]
  val_df1 <- si_dat_fil1[indices$indexOut[[f]],  ] #c("rsi",dep_var)

  model = model
  fold_obs  <- c(fold_obs, val_df1$rsi)
  fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
  #df_stat   <- data.frame( obs = fold_obs, pred = fold_pred)
  #rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  fstat = estimate_fit_stats_folds(fold_obs = fold_obs, fold_pred = fold_pred)
  # rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  rmse_hp[[f]] <- fstat
  print(paste("rmse: ", rmse_hp[f], sep = ""))
}

ran_pred = predict(object = model, data = st_drop_geometry(df))
df$ran_pred = ranp$predictions

#tp = terra::rasterizeGeom(df,res =1000,unit = "m")
df1<- df
df1<- as_Spatial(df1)
sp::gridded(df1)<- TRUE
tp <- terra::rast(df1)

tp1<- tp$ran_pred
#----------------------------------------#
rfp_clip = terra::crop(x = tp1,y = terra::vect(ga_shp), mask = TRUE)

plot(rfp_clip)

hp_df = do.call(rbind, rmse_hp)
gp_df <- hp_df %>% pivot_longer(., cols = -fold, names_to = "Variables",values_to = "value") %>% filter (Variables != c("AIC", "BIC"))
gp_df <- gp_df %>% mutate(var = case_when(Variables == "mean_bias" ~ "Mean Bias (MB)",
                                          Variables  == "MAE"~ "Mean Absolute Error (MAE)",
                                          Variables  == "RMSE" ~ "Root Mean Square Error (RMSE)",
                                          .default  =  "R\U00B2 Adjusted")
)
acc_mb <- gp_df %>% dplyr::filter(Variables == "mean_bias") %>% ggplot(aes (factor(var), value))+geom_violin(trim = FALSE,draw_quantiles = TRUE)+
  # geom_jitter(height = 0, width = 0.1)+
  geom_boxplot(width = 0.5)+
  facet_wrap(vars(var), scales = "free")+
  ylab(" ")+
  xlab (" ")+
  #scale_y_continuous(breaks = c(-2,0,2))+
  theme_bw(10)+
  theme(panel.background = element_rect( fill = "white" , color = "white"),
        panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "white", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        strip.text = element_text(size = 8),
        # strip.background = element_rect(colour = "black", fill = "#F2F2F2", size = 0.2),
        axis.title = element_text(size = rel(1.05) , color = "black"),
        axis.text = element_text(size = rel(1.05) , color = "black"),
        axis.text.x = element_text(size = 14, angle = 360, face = "plain"),
        axis.text.y = element_text(size = 14,face = "plain",angle = 360),
        axis.title.y = element_text(size = 14, face = "plain"),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

acc_mae<- gp_df %>% dplyr::filter(Variables == "MAE") %>% ggplot(aes (factor(var), value))+geom_violin(trim = FALSE,draw_quantiles = TRUE)+
  # geom_jitter(height = 0, width = 0.1)+
  geom_boxplot(width = 0.5)+
  facet_wrap(vars(var), scales = "free")+
  ylab(" ")+
  xlab (" ")+
  #scale_y_continuous(breaks = c(-2,0,2))+
  theme_bw(10)+
  theme(panel.background = element_rect( fill = "white" , color = "white"),
        panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "white", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        strip.text = element_text(size = 8),
        # strip.background = element_rect(colour = "black", fill = "#F2F2F2", size = 0.2),
        axis.title = element_text(size = rel(1.05) , color = "black"),
        axis.text = element_text(size = rel(1.05) , color = "black"),
        axis.text.x = element_text(size = 14, angle = 360, face = "plain"),
        axis.text.y = element_text(size = 14,face = "plain",angle = 360),
        axis.title.y = element_text(size = 14, face = "plain"),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )
acc_rmse<- gp_df %>% dplyr::filter(Variables == "RMSE") %>% ggplot(aes (factor(var), value))+geom_violin(trim = FALSE,draw_quantiles = TRUE)+
  # geom_jitter(height = 0, width = 0.1)+
  geom_boxplot(width = 0.5)+
  facet_wrap(vars(var), scales = "free")+
  ylab(" ")+
  xlab (" ")+
  #scale_y_continuous(breaks = c(-2,0,2))+
  theme_bw(10)+
  theme(panel.background = element_rect( fill = "white" , color = "white"),
        panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "white", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        strip.text = element_text(size = 8),
        # strip.background = element_rect(colour = "black", fill = "#F2F2F2", size = 0.2),
        axis.title = element_text(size = rel(1.05) , color = "black"),
        axis.text = element_text(size = rel(1.05) , color = "black"),
        axis.text.x = element_text(size = 14, angle = 360, face = "plain"),
        axis.text.y = element_text(size = 14,face = "plain",angle = 360),
        axis.title.y = element_text(size = 14, face = "plain"),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

acc_adj_Rsq<- gp_df %>% dplyr::filter(Variables == "adj_Rsq") %>% ggplot(aes (factor(var), value))+geom_violin(trim = FALSE,draw_quantiles = TRUE)+
  # geom_jitter(height = 0, width = 0.1)+
  geom_boxplot(width = 0.5)+
  facet_wrap(vars(var), scales = "free")+
  ylab(" ")+
  xlab (" ")+
  #scale_y_continuous(breaks = c(-2,0,2))+
  theme_bw(10)+
  theme(panel.background = element_rect( fill = "white" , color = "white"),
        panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "white", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        strip.text = element_text(size = 8),
        # strip.background = element_rect(colour = "black", fill = "#F2F2F2", size = 0.2),
        axis.title = element_text(size = rel(1.05) , color = "black"),
        axis.text = element_text(size = rel(1.05) , color = "black"),
        axis.text.x = element_text(size = 14, angle = 360, face = "plain"),
        axis.text.y = element_text(size = 14,face = "plain",angle = 360),
        axis.title.y = element_text(size = 14, face = "plain"),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

# tiff(filename = '../02_01_Output_Prj03/00_01_Combined_Outputs/accuracy_met_rfsi4.tif',width = 12,height = 12,units = "in",res = 600)
# # plot_layout(sp_11|lp_31)/(yp21|pp_21)+
# #   theme(plot.tag.position  = c(0.1, .90))+
# #   plot_annotation(tag_levels = c("A"))
# # dev.off()
#
# ((acc_mb /acc_mae) &
#     theme(plot.tag.position  = c(.22, .96))) -
#
#   ((acc_rmse / acc_adj_Rsq) &
#      theme(plot.tag.position  = c(.22, .96))) +
#
#   plot_annotation(tag_levels = "A",tag_prefix = "(",tag_suffix = ")") &
#   theme(plot.tag = element_text(face = 'bold'))
# dev.off()
#
png(filename = '../02_01_Output_Prj03/00_01_Combined_Outputs/accuracy_met_rfsi4_SINorm_10fold_1vif.png',width = 12,height = 12,units = "in",res = 600)
# plot_layout(sp_11|lp_31)/(yp21|pp_21)+
#   theme(plot.tag.position  = c(0.1, .90))+
#   plot_annotation(tag_levels = c("A"))
# dev.off()

((acc_mb /acc_mae) &
    theme(plot.tag.position  = c(.22, .96))) -

  ((acc_rmse / acc_adj_Rsq) &
     theme(plot.tag.position  = c(.22, .96))) +

  plot_annotation(tag_levels = "A",tag_prefix = "(",tag_suffix = ")") &
  theme(plot.tag = element_text(face = 'bold'))

dev.off()

write.csv(hp_df, "../02_01_Output_Prj03/00_01_Combined_Outputs/Valid_stat_10_fold1vif.csv")

