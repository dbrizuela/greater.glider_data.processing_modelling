### Script to fit Maxent models of the greater glider (Petauroides volans) across its entire distribution range.

### Models are fitted using four different ways of background sampling (see main text):
# random (rndm)
# bias model (bmodel) 
# bias file (bfile)
# target group background (tgb)

### Models are fitted with the 3 splits of fitting - external evaluation data and then final models are fitted using all data.

# load libraries
library(raster)
library(tidyverse)
library(rJava)
# options(java.parameters = "-Xmx1g") # this has to be done before loading dismo library
library(dismo)
library(maptools)
library(blockCV)
library(rgdal)
library(sf)
library(mapview)
library(rgeos)

##################################################################################################
######### Load occ data and background points ####################################################

### Load data if not already in workspace:

# ##### presences (1 per 500m cell) and removed from eval blocks
# ### Sets for external evaluation
# # Set 1
# load("./outputs/processed_data/PB_fit_1.RData")
# # Set 2
# load("./outputs/processed_data/PB_fit_2.RData")
# # Set 3
# load("./outputs/processed_data/PB_fit_3.RData")
# ### Set for final models
# load("./outputs/processed_data/PB_fmodel.RData")
# 
# 
# ##### rndm bckgr points. 1 per 500m cell, removed from eval blocks. sample of 100000
# ### Sets for external evaluation
# # Set 1
# load("./outputs/processed_data/rndm_bckgr_1.RData")
# # Set 2
# load("./outputs/processed_data/rndm_bckgr_2.RData")
# # Set 3
# load("./outputs/processed_data/rndm_bckgr_3.RData")
# ### Set for final models
# load("./outputs/processed_data/rndm_bckgr_fmodel.RData") 
# 
# 
# ##### bmodel bckgr points. 1 per 500m cell, removed from eval blocks. sample of 100000
# ### Sets for external evaluation
# # Set 1
# load("./outputs/processed_data/bmodel_bckgr_1.RData")
# # Set 2
# load("./outputs/processed_data/bmodel_bckgr_2.RData")
# # Set 3
# load("./outputs/processed_data/bmodel_bckgr_3.RData") 
# ### Set for final models
# load("./outputs/processed_data/bmodel_bckgr_fmodel.RData") 
# 
# 
# ##### bfile bckgr points. 1 per 500m cell, removed from eval blocks. sample of 100000
# ### Sets for external evaluation
# # Set 1
# load("./outputs/processed_data/bfile_bckgr_1.RData")
# # Set 2
# load("./outputs/processed_data/bfile_bckgr_2.RData")
# # Set 3
# load("./outputs/processed_data/bfile_bckgr_3.RData")
# ### Set for final models
# load("./outputs/processed_data/bfile_bckgr_fmodel.RData") 
# 
# 
# ##### tgb bckgr points. 1 per 500m cell, removed from eval blocks.
# ### Sets for external evaluation
# # Set 1
# load("./outputs/processed_data/tgb_bckgr_1.RData")
# # Set 2
# load("./outputs/processed_data/tgb_bckgr_2.RData")
# # Set 3
# load("./outputs/processed_data/tgb_bckgr_3.RData")
# ### Set for final models
# load("./outputs/processed_data/tgb_bckgr_fmodel.RData")
# 
# 
# ##### eval dataset all
# ### evaluation sets for external evaluation
# # Set 1
# load("./outputs/processed_data/PA_PIA_eval_1.RData")
# # Set 2
# load("./outputs/processed_data/PA_PIA_eval_2.RData")
# # Set 3
# load("./outputs/processed_data/PA_PIA_eval_3.RData")

# ### External evaluation sets Southern half
# # Set 1
# load("./outputs/processed_data/PA_PIA_eval_1_south.RData")
# # Set 2
# load("./outputs/processed_data/PA_PIA_eval_2_south.RData")
# # Set 3
# load("./outputs/processed_data/PA_PIA_eval_3_south.RData")

# # for some reason everytime I upload these again, PA column becomes factor, so change to numeric:
# PA_PIA_eval_1$PA <- as.numeric(levels(PA_PIA_eval_1$PA))[PA_PIA_eval_1$PA]
# PA_PIA_eval_2$PA <- as.numeric(levels(PA_PIA_eval_2$PA))[PA_PIA_eval_2$PA]
# PA_PIA_eval_3$PA <- as.numeric(levels(PA_PIA_eval_3$PA))[PA_PIA_eval_3$PA]
# 
# PA_PIA_eval_1_south$PA <- as.numeric(levels(PA_PIA_eval_1_south$PA))[PA_PIA_eval_1_south$PA]
# PA_PIA_eval_2_south$PA <- as.numeric(levels(PA_PIA_eval_2_south$PA))[PA_PIA_eval_2_south$PA]
# PA_PIA_eval_3_south$PA <- as.numeric(levels(PA_PIA_eval_3_south$PA))[PA_PIA_eval_3_south$PA]




##################################################################################################
######### Load and prepare env var  ##############################################################

### This is done in the 'data_processing' script as well:

# # All variables have been resampled (500 m), projected (epsg:3577) and masked to modelling extent
# vars_stack_all <- raster::stack(list.files(paste0(getwd(),"./spatial_data/environmental_variables"), pattern = '.asc', full.names=TRUE))
# crs(vars_stack_all) <- ("+init=epsg:3577")
# 
# names(vars_stack_all)
# 
# # Change names to variables in stack
# vars.names <- c("temp_mean", 
#                 "temp_warm", 
#                 "temp_cold",
#                 "pp_annual", 
#                 "pp_driest", 
#                 "temp_season", 
#                 "pp_season", 
#                 "mean_drought", 
#                 "fPAR_mean", 
#                 "fPAR_var", 
#                 "GPP_mean", 
#                 "GPP_var", 
#                 "mean_EDI", 
#                 "tsf", 
#                 "c_height",
#                 "mean_xveg")
# 
# names(vars_stack_all) <- vars.names
# names(vars_stack_all)
# 
# # pairs(vars_stack_all, maxpixels = 10000) # check correlation
# 
# See 'data_processing' script for reasoning on variables' removal
#
#
# # drop variables
# names(vars_stack_all)
# vars_stack <- dropLayer(vars_stack_all, c(1,3,8,11))
# names(vars_stack)
# rm(vars_stack_all)




##################################################################################################
######### prepare 'sites with data' (SWD) objects for model fitting   ############################

##### random
### Sets for external evaluation
# Set 1
SWD_rndm_1 <- rbind(PB_fit_1, rndm_bckgr_1)
id_rndm_1 <- c(rep(1, nrow(PB_fit_1)), rep(0, nrow(rndm_bckgr_1)))
# Set 2
SWD_rndm_2 <- rbind(PB_fit_2, rndm_bckgr_2)
id_rndm_2 <- c(rep(1, nrow(PB_fit_2)), rep(0, nrow(rndm_bckgr_2)))
# Set 3
SWD_rndm_3 <- rbind(PB_fit_3, rndm_bckgr_3)
id_rndm_3 <- c(rep(1, nrow(PB_fit_3)), rep(0, nrow(rndm_bckgr_3)))
### Set for final model
SWD_rndm_fmodel <- rbind(PB_fmodel, rndm_bckgr_fmodel)
id_rndm_fmodel <- c(rep(1, nrow(PB_fmodel)), rep(0, nrow(rndm_bckgr_fmodel)))


##### bias model
### Sets for external evaluation
# Set 1
SWD_bmodel_1 <- rbind(PB_fit_1, bmodel_bckgr_1)
id_bmodel_1 <- c(rep(1, nrow(PB_fit_1)), rep(0, nrow(bmodel_bckgr_1)))
# Set 2
SWD_bmodel_2 <- rbind(PB_fit_2, bmodel_bckgr_2)
id_bmodel_2 <- c(rep(1, nrow(PB_fit_2)), rep(0, nrow(bmodel_bckgr_2)))
# Set 3
SWD_bmodel_3 <- rbind(PB_fit_3, bmodel_bckgr_3)
id_bmodel_3 <- c(rep(1, nrow(PB_fit_3)), rep(0, nrow(bmodel_bckgr_3)))
### Set for final model
SWD_bmodel_fmodel <- rbind(PB_fmodel, bmodel_bckgr_fmodel)
id_bmodel_fmodel <- c(rep(1, nrow(PB_fmodel)), rep(0, nrow(bmodel_bckgr_fmodel)))


##### bfile
### Sets for external evaluation
# Set 1
SWD_bfile_1 <- rbind(PB_fit_1, bfile_bckgr_1)
id_bfile_1 <- c(rep(1, nrow(PB_fit_1)), rep(0, nrow(bfile_bckgr_1)))
# Set 2
SWD_bfile_2 <- rbind(PB_fit_2, bfile_bckgr_2)
id_bfile_2 <- c(rep(1, nrow(PB_fit_2)), rep(0, nrow(bfile_bckgr_2)))
# Set 3
SWD_bfile_3 <- rbind(PB_fit_3, bfile_bckgr_3)
id_bfile_3 <- c(rep(1, nrow(PB_fit_3)), rep(0, nrow(bfile_bckgr_3)))
### Set for final model
SWD_bfile_fmodel <- rbind(PB_fmodel, bfile_bckgr_fmodel)
id_bfile_fmodel <- c(rep(1, nrow(PB_fmodel)), rep(0, nrow(bfile_bckgr_fmodel)))


##### tgb
### Sets for external evaluation
# Set 1
SWD_tgb_1 <- rbind(PB_fit_1, tgb_bckgr_1)
id_tgb_1 <- c(rep(1, nrow(PB_fit_1)), rep(0, nrow(tgb_bckgr_1)))
# Set 2
SWD_tgb_2 <- rbind(PB_fit_2, tgb_bckgr_2)
id_tgb_2 <- c(rep(1, nrow(PB_fit_2)), rep(0, nrow(tgb_bckgr_2)))
# Set 3
SWD_tgb_3 <- rbind(PB_fit_3, tgb_bckgr_3)
id_tgb_3 <- c(rep(1, nrow(PB_fit_3)), rep(0, nrow(tgb_bckgr_3)))
### Set for final model
SWD_tgb_fmodel <- rbind(PB_fmodel, tgb_bckgr_fmodel)
id_tgb_fmodel <- c(rep(1, nrow(PB_fmodel)), rep(0, nrow(tgb_bckgr_fmodel)))





##################################################################################################
######### create blocks for internal evaluation with blockCV #####################################

# Internal evaluation is doing using all data (i.e. the data for the final models)

### random background
locations_rndm <- as.data.frame(rbind(PB_fmodel[4:5], rndm_bckgr_fmodel[4:5])) # 
loc4block_rndm <- data.frame(id_rndm_fmodel, locations_rndm, SWD_rndm_fmodel)  #  ID for random background
loc.df_rndm <- SpatialPointsDataFrame(loc4block_rndm[,c("x", "y")], loc4block_rndm)
crs(loc.df_rndm) <- ("+init=epsg:3577")

### bmodel background
locations_bmodel <- as.data.frame(rbind(PB_fmodel[4:5], bmodel_bckgr_fmodel[4:5])) # 
loc4block_bmodel <- data.frame(id_bmodel_fmodel, locations_bmodel, SWD_bmodel_fmodel)  #  ID for bmodel background
loc.df_bmodel <- SpatialPointsDataFrame(loc4block_bmodel[,c("x", "y")], loc4block_bmodel)
crs(loc.df_bmodel) <- ("+init=epsg:3577")

### bfile background
locations_bfile <- as.data.frame(rbind(PB_fmodel[4:5], bfile_bckgr_fmodel[4:5])) # 
loc4block_bfile <- data.frame(id_bfile_fmodel, locations_bfile, SWD_bfile_fmodel)  #  ID for bfile background
loc.df_bfile <- SpatialPointsDataFrame(loc4block_bfile[,c("x", "y")], loc4block_bfile)
crs(loc.df_bfile) <- ("+init=epsg:3577")

### tgb background
locations_tgb <- as.data.frame(rbind(PB_fmodel[4:5], tgb_bckgr_fmodel[4:5])) # 
loc4block_tgb <- data.frame(id_tgb_fmodel, locations_tgb, SWD_tgb_fmodel)  #  ID for tgb
loc.df_tgb <- SpatialPointsDataFrame(loc4block_tgb[,c("x", "y")], loc4block_tgb)
crs(loc.df_tgb) <- ("+init=epsg:3577")


#### spatial blocking 

### random 
sb_rndm <- spatialBlock(speciesData = loc.df_rndm,
                        species = "id_rndm_fmodel",
                        rasterLayer = vars_stack[[1]],
                        theRange = 25000, # 25 km
                        k = 10,
                        selection = 'random',
                        maskBySpecies = T,
                        showBlocks = T,
                        biomod2Format = FALSE)
# save(sb_rndm, file = "./outputs/spatial_blocks/sb_rndm.RData") # advisable to save blocks

### bmodel 
sb_bmodel <- spatialBlock(speciesData = loc.df_bmodel,
                          species = "id_bmodel_fmodel",
                          rasterLayer = vars_stack[[1]],
                          theRange = 25000, # 25 km
                          k = 10,
                          selection = 'random',
                          maskBySpecies = T,
                          showBlocks = T,
                          biomod2Format = FALSE)
# save(sb_bmodel, file = "./outputs/spatial_blocks/sb_bmodel.RData") # advisable to save blocks

### bfile 
sb_bfile <- spatialBlock(speciesData = loc.df_bfile,
                         species = "id_bfile_fmodel",
                         rasterLayer = vars_stack[[1]],
                         theRange = 25000, # 25 km
                         k = 10,
                         selection = 'random',
                         maskBySpecies = T,
                         showBlocks = T,
                         biomod2Format = FALSE)
# save(sb_bfile, file = "./outputs/spatial_blocks/sb_bfile.RData") # advisable to save blocks

### tgb 
sb_tgb <- spatialBlock(speciesData = loc.df_tgb,
                       species = "id_tgb_fmodel",
                       rasterLayer = vars_stack[[1]],
                       theRange = 25000, # 25 km
                       k = 10,
                       selection = 'random',
                       maskBySpecies = T,
                       showBlocks = T,
                       biomod2Format = FALSE)
# save(sb_tgb, file = "./outputs/spatial_blocks/sb_tgb.RData") # advisable to save blocks





##################################################################################################
######### Internal evaluation ####################################################################


#### Random background

folds_rndm <- sb_rndm$folds # Get the fold information
# create the structures to hold the results for each model
k <- 10
rndm_blockCV_AllFC <- vector('list', k)
rndm_blockCV_LQ <- vector('list', k)
rndm_blockCV_LQH <- vector('list', k)

for (i in 1:k) {
  trainSet_rndm <- unlist(folds_rndm[[i]][1]) # extract the training set indices
  testSet_rndm <- unlist(folds_rndm[[i]][2]) # extract the testing set indices
  
  Maxent_rndm_blockCV_AllFC <- maxent(x=SWD_rndm_fmodel[7:18][trainSet_rndm,], p=id_rndm_fmodel[trainSet_rndm], removeDuplicates=FALSE)
  Maxent_rndm_blockCV_LQ <- maxent(x=SWD_rndm_fmodel[7:18][trainSet_rndm,], p=id_rndm_fmodel[trainSet_rndm], removeDuplicates=FALSE, args=c('noautofeature', 'nohinge', 'noproduct'))
  Maxent_rndm_blockCV_LQH <- maxent(x=SWD_rndm_fmodel[7:18][trainSet_rndm,], p=id_rndm_fmodel[trainSet_rndm], removeDuplicates=FALSE, args=c('noproduct'))
  testdat.env_rndm <- SWD_rndm_fmodel[7:18][testSet_rndm,]
  testdat.sp_rndm <- id_rndm_fmodel[testSet_rndm]
  locust.test.pres_rndm <- testdat.env_rndm[testdat.sp_rndm==1,]
  locust.test.bg_rndm  <- testdat.env_rndm[testdat.sp_rndm==0,]
  
  # evaluate model with held-out data
  rndm_blockCV_AllFC[[i]] <- evaluate(p=locust.test.pres_rndm, a=locust.test.bg_rndm, model=Maxent_rndm_blockCV_AllFC)
  rndm_blockCV_LQ[[i]] <- evaluate(p=locust.test.pres_rndm, a=locust.test.bg_rndm, model=Maxent_rndm_blockCV_LQ)
  rndm_blockCV_LQH[[i]] <- evaluate(p=locust.test.pres_rndm, a=locust.test.bg_rndm, model=Maxent_rndm_blockCV_LQH)
}


#### Bmodel background

folds_bmodel <- sb_bmodel$folds # Get the fold information
# create the structures to hold the results for each model
k <- 10
bmodel_blockCV_AllFC <- vector('list', k)
bmodel_blockCV_LQ <- vector('list', k)
bmodel_blockCV_LQH <- vector('list', k)

for (i in 1:k) {
  trainSet_bmodel <- unlist(folds_bmodel[[i]][1]) # extract the training set indices
  testSet_bmodel <- unlist(folds_bmodel[[i]][2]) # extract the testing set indices
  
  Maxent_bmodel_blockCV_AllFC <- maxent(x=SWD_bmodel_fmodel[7:18][trainSet_bmodel,], p=id_bmodel_fmodel[trainSet_bmodel], removeDuplicates=FALSE)
  Maxent_bmodel_blockCV_LQ <- maxent(x=SWD_bmodel_fmodel[7:18][trainSet_bmodel,], p=id_bmodel_fmodel[trainSet_bmodel], removeDuplicates=FALSE, args=c('noautofeature', 'nohinge', 'noproduct'))
  Maxent_bmodel_blockCV_LQH <- maxent(x=SWD_bmodel_fmodel[7:18][trainSet_bmodel,], p=id_bmodel_fmodel[trainSet_bmodel], removeDuplicates=FALSE, args=c('noproduct'))
  testdat.env_bmodel <- SWD_bmodel_fmodel[7:18][testSet_bmodel,]
  testdat.sp_bmodel <- id_bmodel_fmodel[testSet_bmodel]
  locust.test.pres_bmodel <- testdat.env_bmodel[testdat.sp_bmodel==1,]
  locust.test.bg_bmodel  <- testdat.env_bmodel[testdat.sp_bmodel==0,]
  
  # evaluate model with held-out data
  bmodel_blockCV_AllFC[[i]] <- evaluate(p=locust.test.pres_bmodel, a=locust.test.bg_bmodel, model=Maxent_bmodel_blockCV_AllFC)
  bmodel_blockCV_LQ[[i]] <- evaluate(p=locust.test.pres_bmodel, a=locust.test.bg_bmodel, model=Maxent_bmodel_blockCV_LQ)
  bmodel_blockCV_LQH[[i]] <- evaluate(p=locust.test.pres_bmodel, a=locust.test.bg_bmodel, model=Maxent_bmodel_blockCV_LQH)
}



#### Bfile background

folds_bfile <- sb_bfile$folds # Get the fold information
# create the structures to hold the results for each model
k <- 10
bfile_blockCV_AllFC <- vector('list', k)
bfile_blockCV_LQ <- vector('list', k)
bfile_blockCV_LQH <- vector('list', k)

for (i in 1:k) {
  trainSet_bfile <- unlist(folds_bfile[[i]][1]) # extract the training set indices
  testSet_bfile <- unlist(folds_bfile[[i]][2]) # extract the testing set indices
  
  Maxent_bfile_blockCV_AllFC <- maxent(x=SWD_bfile_fmodel[7:18][trainSet_bfile,], p=id_bfile_fmodel[trainSet_bfile], removeDuplicates=FALSE)
  Maxent_bfile_blockCV_LQ <- maxent(x=SWD_bfile_fmodel[7:18][trainSet_bfile,], p=id_bfile_fmodel[trainSet_bfile], removeDuplicates=FALSE, args=c('noautofeature', 'nohinge', 'noproduct'))
  Maxent_bfile_blockCV_LQH <- maxent(x=SWD_bfile_fmodel[7:18][trainSet_bfile,], p=id_bfile_fmodel[trainSet_bfile], removeDuplicates=FALSE, args=c('noproduct'))
  testdat.env_bfile <- SWD_bfile_fmodel[7:18][testSet_bfile,]
  testdat.sp_bfile <- id_bfile_fmodel[testSet_bfile]
  locust.test.pres_bfile <- testdat.env_bfile[testdat.sp_bfile==1,]
  locust.test.bg_bfile  <- testdat.env_bfile[testdat.sp_bfile==0,]
  
  # evaluate model with held-out data
  bfile_blockCV_AllFC[[i]] <- evaluate(p=locust.test.pres_bfile, a=locust.test.bg_bfile, model=Maxent_bfile_blockCV_AllFC)
  bfile_blockCV_LQ[[i]] <- evaluate(p=locust.test.pres_bfile, a=locust.test.bg_bfile, model=Maxent_bfile_blockCV_LQ)
  bfile_blockCV_LQH[[i]] <- evaluate(p=locust.test.pres_bfile, a=locust.test.bg_bfile, model=Maxent_bfile_blockCV_LQH)
}



#### Tgb background

folds_tgb <- sb_tgb$folds # Get the fold information
# create the structures to hold the results for each model
k <- 10
tgb_blockCV_AllFC <- vector('list', k)
tgb_blockCV_LQ <- vector('list', k)
tgb_blockCV_LQH <- vector('list', k)

for (i in 1:k) {
  trainSet_tgb <- unlist(folds_tgb[[i]][1]) # extract the training set indices
  testSet_tgb <- unlist(folds_tgb[[i]][2]) # extract the testing set indices
  
  Maxent_tgb_blockCV_AllFC <- maxent(x=SWD_tgb_fmodel[7:18][trainSet_tgb,], p=id_tgb_fmodel[trainSet_tgb], removeDuplicates=FALSE)
  Maxent_tgb_blockCV_LQ <- maxent(x=SWD_tgb_fmodel[7:18][trainSet_tgb,], p=id_tgb_fmodel[trainSet_tgb], removeDuplicates=FALSE, args=c('noautofeature', 'nohinge', 'noproduct'))
  Maxent_tgb_blockCV_LQH <- maxent(x=SWD_tgb_fmodel[7:18][trainSet_tgb,], p=id_tgb_fmodel[trainSet_tgb], removeDuplicates=FALSE, args=c('noproduct'))
  testdat.env_tgb <- SWD_tgb_fmodel[7:18][testSet_tgb,]
  testdat.sp_tgb <- id_tgb_fmodel[testSet_tgb]
  locust.test.pres_tgb <- testdat.env_tgb[testdat.sp_tgb==1,]
  locust.test.bg_tgb  <- testdat.env_tgb[testdat.sp_tgb==0,]
  
  # evaluate model with held-out data
  tgb_blockCV_AllFC[[i]] <- evaluate(p=locust.test.pres_tgb, a=locust.test.bg_tgb, model=Maxent_tgb_blockCV_AllFC)
  tgb_blockCV_LQ[[i]] <- evaluate(p=locust.test.pres_tgb, a=locust.test.bg_tgb, model=Maxent_tgb_blockCV_LQ)
  tgb_blockCV_LQH[[i]] <- evaluate(p=locust.test.pres_tgb, a=locust.test.bg_tgb, model=Maxent_tgb_blockCV_LQH)
}

save.image() # advisable to save workspace


#### View AUCs for all models 

# rndm
AUC.rndm_blockCV_AllFC <- sapply(rndm_blockCV_AllFC, slot, 'auc')
AUC.rndm_blockCV_LQH <- sapply(rndm_blockCV_LQH, slot, 'auc')
AUC.rndm_blockCV_LQ <- sapply(rndm_blockCV_LQ, slot, 'auc')

round(mean(AUC.rndm_blockCV_AllFC),3)
round(mean(AUC.rndm_blockCV_LQH),3)
round(mean(AUC.rndm_blockCV_LQ),3)

round(sd(AUC.rndm_blockCV_AllFC),2)
round(sd(AUC.rndm_blockCV_LQH),2)
round(sd(AUC.rndm_blockCV_LQ),2)


#  bmodel
AUC.bmodel_blockCV_AllFC <- sapply(bmodel_blockCV_AllFC, slot, 'auc')
AUC.bmodel_blockCV_LQH <- sapply(bmodel_blockCV_LQH, slot, 'auc')
AUC.bmodel_blockCV_LQ <- sapply(bmodel_blockCV_LQ, slot, 'auc')

round(mean(AUC.bmodel_blockCV_AllFC),3)
round(mean(AUC.bmodel_blockCV_LQH),3)
round(mean(AUC.bmodel_blockCV_LQ),3)

round(sd(AUC.bmodel_blockCV_AllFC),2)
round(sd(AUC.bmodel_blockCV_LQH),2)
round(sd(AUC.bmodel_blockCV_LQ),2)


# bfile
AUC.bfile_blockCV_AllFC <- sapply(bfile_blockCV_AllFC, slot, 'auc')
AUC.bfile_blockCV_LQH <- sapply(bfile_blockCV_LQH, slot, 'auc')
AUC.bfile_blockCV_LQ <- sapply(bfile_blockCV_LQ, slot, 'auc')

round(mean(AUC.bfile_blockCV_AllFC),3)
round(mean(AUC.bfile_blockCV_LQH),3)
round(mean(AUC.bfile_blockCV_LQ),3)

round(sd(AUC.bfile_blockCV_AllFC),2)
round(sd(AUC.bfile_blockCV_LQH),2)
round(sd(AUC.bfile_blockCV_LQ),2)


# tgb
AUC.tgb_blockCV_AllFC <- sapply(tgb_blockCV_AllFC, slot, 'auc')
AUC.tgb_blockCV_LQH <- sapply(tgb_blockCV_LQH, slot, 'auc')
AUC.tgb_blockCV_LQ <- sapply(tgb_blockCV_LQ, slot, 'auc')

round(mean(AUC.tgb_blockCV_AllFC),3)
round(mean(AUC.tgb_blockCV_LQH),3)
round(mean(AUC.tgb_blockCV_LQ),3)

round(sd(AUC.tgb_blockCV_AllFC),2)
round(sd(AUC.tgb_blockCV_LQH),2)
round(sd(AUC.tgb_blockCV_LQ),2)


### Given Maxent stochasticity and random blocks/folds selection, results of these metric will vary from the results reported in the main text.

# rm blockCV maxent obecjts to make project lighter. Once evaluation metrics are obtained, won't next the maxent objects anymore.
# rm(Maxent_rndm_blockCV_AllFC, Maxent_rndm_blockCV_LQ, Maxent_rndm_blockCV_LQH, Maxent_bmodel_blockCV_AllFC, Maxent_bmodel_blockCV_LQ, Maxent_bmodel_blockCV_LQH, Maxent_bfile_blockCV_AllFC, Maxent_bfile_blockCV_LQ, Maxent_bfile_blockCV_LQH, Maxent_tgb_blockCV_AllFC, Maxent_tgb_blockCV_LQ, Maxent_tgb_blockCV_LQH) 

# remove spatial blocks as well (these are saved already)
# rm(sb_bfile, sb_bmodel, sb_rndm, sb_tgb)





##################################################################################################
######### External evaluation ####################################################################

# 3 sets of models (one for each ext eval block) for each background sampling strategy



##### Random background
### Set_1
Maxent_rndm_allFC_1 <- maxent(x=SWD_rndm_1[7:18], p=id_rndm_1, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_rndm_LQ_1 <- maxent(x=SWD_rndm_1[7:18], p=id_rndm_1,
                           args=c('-P', '-J', 
                                  'noautofeature', 
                                  'nohinge', 
                                  'noproduct'), removeDuplicates=FALSE)
Maxent_rndm_LQH_1 <- maxent(x=SWD_rndm_1[7:18], p=id_rndm_1,
                            args=c('-P', '-J',
                                   'noproduct'), removeDuplicates=FALSE)

### Set_2
Maxent_rndm_allFC_2 <- maxent(x=SWD_rndm_2[7:18], p=id_rndm_2, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_rndm_LQ_2 <- maxent(x=SWD_rndm_2[7:18], p=id_rndm_2,
                           args=c('-P', '-J', 
                                  'noautofeature', 
                                  'nohinge', 
                                  'noproduct'), removeDuplicates=FALSE)
Maxent_rndm_LQH_2 <- maxent(x=SWD_rndm_2[7:18], p=id_rndm_2,
                            args=c('-P', '-J',
                                   'noproduct'), removeDuplicates=FALSE)

### Set_3
Maxent_rndm_allFC_3 <- maxent(x=SWD_rndm_3[7:18], p=id_rndm_3, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_rndm_LQ_3 <- maxent(x=SWD_rndm_3[7:18], p=id_rndm_3,
                           args=c('-P', '-J', 
                                  'noautofeature', 
                                  'nohinge', 
                                  'noproduct'), removeDuplicates=FALSE)
Maxent_rndm_LQH_3 <- maxent(x=SWD_rndm_3[7:18], p=id_rndm_3,
                            args=c('-P', '-J',
                                   'noproduct'), removeDuplicates=FALSE)



##### Bmodel background
### Set_1
Maxent_bmodel_allFC_1 <- maxent(x=SWD_bmodel_1[7:18], p=id_bmodel_1, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bmodel_LQ_1 <- maxent(x=SWD_bmodel_1[7:18], p=id_bmodel_1,
                             args=c('-P', '-J', 
                                    'noautofeature', 
                                    'nohinge', 
                                    'noproduct'), removeDuplicates=FALSE)
Maxent_bmodel_LQH_1 <- maxent(x=SWD_bmodel_1[7:18], p=id_bmodel_1,
                              args=c('-P', '-J',
                                     'noproduct'), removeDuplicates=FALSE)

### Set_2
Maxent_bmodel_allFC_2 <- maxent(x=SWD_bmodel_2[7:18], p=id_bmodel_2, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bmodel_LQ_2 <- maxent(x=SWD_bmodel_2[7:18], p=id_bmodel_2,
                             args=c('-P', '-J', 
                                    'noautofeature', 
                                    'nohinge', 
                                    'noproduct'), removeDuplicates=FALSE)
Maxent_bmodel_LQH_2 <- maxent(x=SWD_bmodel_2[7:18], p=id_bmodel_2,
                              args=c('-P', '-J',
                                     'noproduct'), removeDuplicates=FALSE)

### Set_3
Maxent_bmodel_allFC_3 <- maxent(x=SWD_bmodel_3[7:18], p=id_bmodel_3, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bmodel_LQ_3 <- maxent(x=SWD_bmodel_3[7:18], p=id_bmodel_3,
                             args=c('-P', '-J', 
                                    'noautofeature', 
                                    'nohinge', 
                                    'noproduct'), removeDuplicates=FALSE)
Maxent_bmodel_LQH_3 <- maxent(x=SWD_bmodel_3[7:18], p=id_bmodel_3,
                              args=c('-P', '-J',
                                     'noproduct'), removeDuplicates=FALSE)



##### Bfile background
### Set_1
Maxent_bfile_allFC_1 <- maxent(x=SWD_bfile_1[7:18], p=id_bfile_1, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bfile_LQ_1 <- maxent(x=SWD_bfile_1[7:18], p=id_bfile_1,
                            args=c('-P', '-J', 
                                   'noautofeature', 
                                   'nohinge', 
                                   'noproduct'), removeDuplicates=FALSE)
Maxent_bfile_LQH_1 <- maxent(x=SWD_bfile_1[7:18], p=id_bfile_1,
                             args=c('-P', '-J',
                                    'noproduct'), removeDuplicates=FALSE)

### Set_2
Maxent_bfile_allFC_2 <- maxent(x=SWD_bfile_2[7:18], p=id_bfile_2, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bfile_LQ_2 <- maxent(x=SWD_bfile_2[7:18], p=id_bfile_2,
                            args=c('-P', '-J', 
                                   'noautofeature', 
                                   'nohinge', 
                                   'noproduct'), removeDuplicates=FALSE)
Maxent_bfile_LQH_2 <- maxent(x=SWD_bfile_2[7:18], p=id_bfile_2,
                             args=c('-P', '-J',
                                    'noproduct'), removeDuplicates=FALSE)

### Set_3
Maxent_bfile_allFC_3 <- maxent(x=SWD_bfile_3[7:18], p=id_bfile_3, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_bfile_LQ_3 <- maxent(x=SWD_bfile_3[7:18], p=id_bfile_3,
                            args=c('-P', '-J', 
                                   'noautofeature', 
                                   'nohinge', 
                                   'noproduct'), removeDuplicates=FALSE)
Maxent_bfile_LQH_3 <- maxent(x=SWD_bfile_3[7:18], p=id_bfile_3,
                             args=c('-P', '-J',
                                    'noproduct'), removeDuplicates=FALSE)


##### Tgb background
### Set_1
Maxent_tgb_allFC_1 <- maxent(x=SWD_tgb_1[7:18], p=id_tgb_1, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_tgb_LQ_1 <- maxent(x=SWD_tgb_1[7:18], p=id_tgb_1,
                          args=c('-P', '-J', 
                                 'noautofeature', 
                                 'nohinge', 
                                 'noproduct'), removeDuplicates=FALSE)
Maxent_tgb_LQH_1 <- maxent(x=SWD_tgb_1[7:18], p=id_tgb_1,
                           args=c('-P', '-J',
                                  'noproduct'), removeDuplicates=FALSE)

### Set_2
Maxent_tgb_allFC_2 <- maxent(x=SWD_tgb_2[7:18], p=id_tgb_2, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_tgb_LQ_2 <- maxent(x=SWD_tgb_2[7:18], p=id_tgb_2,
                          args=c('-P', '-J', 
                                 'noautofeature', 
                                 'nohinge', 
                                 'noproduct'), removeDuplicates=FALSE)
Maxent_tgb_LQH_2 <- maxent(x=SWD_tgb_2[7:18], p=id_tgb_2,
                           args=c('-P', '-J',
                                  'noproduct'), removeDuplicates=FALSE)

### Set_3
Maxent_tgb_allFC_3 <- maxent(x=SWD_tgb_3[7:18], p=id_tgb_3, args=c('-P','-J'), removeDuplicates=FALSE)
Maxent_tgb_LQ_3 <- maxent(x=SWD_tgb_3[7:18], p=id_tgb_3,
                          args=c('-P', '-J', 
                                 'noautofeature', 
                                 'nohinge', 
                                 'noproduct'), removeDuplicates=FALSE)
Maxent_tgb_LQH_3 <- maxent(x=SWD_tgb_3[7:18], p=id_tgb_3,
                           args=c('-P', '-J',
                                  'noproduct'), removeDuplicates=FALSE)



######### Evaluation.

# Evaluation was done manually for each evaluation set of each background sampling strategy, but this could be writen as a function or a loop (making sure each model is evaluated with the corresponding evaluation set).



##### Random background

### Set_1
rndm_AllFC_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_rndm_allFC_1)

rndm_LQH_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_rndm_LQH_1)

rndm_LQ_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_rndm_LQ_1)


### Set_2
rndm_AllFC_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_rndm_allFC_2)

rndm_LQH_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_rndm_LQH_2)

rndm_LQ_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_rndm_LQ_2)


### Set_3
rndm_AllFC_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_rndm_allFC_3)

rndm_LQH_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_rndm_LQH_3)

rndm_LQ_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_rndm_LQ_3)



##### Bmodel background

### Set_1
bmodel_AllFC_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_1)

bmodel_LQH_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_1)

bmodel_LQ_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_1)


### Set_2
bmodel_AllFC_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_2)

bmodel_LQH_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_2)

bmodel_LQ_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_2)


### Set_3
bmodel_AllFC_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_3)

bmodel_LQH_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_3)

bmodel_LQ_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_3)



##### Bfile background

### Set_1
bfile_AllFC_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bfile_allFC_1)

bfile_LQH_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bfile_LQH_1)

bfile_LQ_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_bfile_LQ_1)


### Set_2
bfile_AllFC_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bfile_allFC_2)

bfile_LQH_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bfile_LQH_2)

bfile_LQ_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_bfile_LQ_2)


### Set_3
bfile_AllFC_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bfile_allFC_3)

bfile_LQH_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bfile_LQH_3)

bfile_LQ_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_bfile_LQ_3)



##### Tgb background
### Set_1
tgb_AllFC_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_tgb_allFC_1)

tgb_LQH_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_tgb_LQH_1)

tgb_LQ_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = Maxent_tgb_LQ_1)


### Set_2
tgb_AllFC_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_tgb_allFC_2)

tgb_LQH_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_tgb_LQH_2)

tgb_LQ_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = Maxent_tgb_LQ_2)


### Set_3
tgb_AllFC_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_tgb_allFC_3)

tgb_LQH_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_tgb_LQH_3)

tgb_LQ_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = Maxent_tgb_LQ_3)


##### results

### rndm
round(rndm_AllFC_eval_1@auc, 3)
round(rndm_LQH_eval_1@auc, 3)
round(rndm_LQ_eval_1@auc, 3)

round(rndm_AllFC_eval_2@auc, 3)
round(rndm_LQH_eval_2@auc, 3)
round(rndm_LQ_eval_2@auc, 3)

round(rndm_AllFC_eval_3@auc, 3)
round(rndm_LQH_eval_3@auc, 3)
round(rndm_LQ_eval_3@auc, 3)

### mean
round(mean(c(rndm_AllFC_eval_1@auc, rndm_AllFC_eval_2@auc, rndm_AllFC_eval_3@auc)),3) # All_FC
round(mean(c(rndm_LQH_eval_1@auc, rndm_LQH_eval_2@auc, rndm_LQH_eval_3@auc)),3) # LQH
round(mean(c(rndm_LQ_eval_1@auc, rndm_LQ_eval_2@auc, rndm_LQ_eval_3@auc)),3) # LQ


### bmodel
round(bmodel_AllFC_eval_1@auc, 3)
round(bmodel_LQH_eval_1@auc, 3)
round(bmodel_LQ_eval_1@auc, 3)

round(bmodel_AllFC_eval_2@auc, 3)
round(bmodel_LQH_eval_2@auc, 3)
round(bmodel_LQ_eval_2@auc, 3)

round(bmodel_AllFC_eval_3@auc, 3)
round(bmodel_LQH_eval_3@auc, 3)
round(bmodel_LQ_eval_3@auc, 3)

### mean
round(mean(c(bmodel_AllFC_eval_1@auc, bmodel_AllFC_eval_2@auc, bmodel_AllFC_eval_3@auc)),3) # All_FC
round(mean(c(bmodel_LQH_eval_1@auc, bmodel_LQH_eval_2@auc, bmodel_LQH_eval_3@auc)),3) # LQH
round(mean(c(bmodel_LQ_eval_1@auc, bmodel_LQ_eval_2@auc, bmodel_LQ_eval_3@auc)),3) # LQ

### bfile
round(bfile_AllFC_eval_1@auc, 3)
round(bfile_LQH_eval_1@auc, 3)
round(bfile_LQ_eval_1@auc, 3)

round(bfile_AllFC_eval_2@auc, 3)
round(bfile_LQH_eval_2@auc, 3)
round(bfile_LQ_eval_2@auc, 3)

round(bfile_AllFC_eval_3@auc, 3)
round(bfile_LQH_eval_3@auc, 3)
round(bfile_LQ_eval_3@auc, 3)

### mean
round(mean(c(bfile_AllFC_eval_1@auc, bfile_AllFC_eval_2@auc, bfile_AllFC_eval_3@auc)),3) # All_FC
round(mean(c(bfile_LQH_eval_1@auc, bfile_LQH_eval_2@auc, bfile_LQH_eval_3@auc)),3) # LQH
round(mean(c(bfile_LQ_eval_1@auc, bfile_LQ_eval_2@auc, bfile_LQ_eval_3@auc)),3) # LQ

### tgb
round(tgb_AllFC_eval_1@auc, 3)
round(tgb_LQH_eval_1@auc, 3)
round(tgb_LQ_eval_1@auc, 3)

round(tgb_AllFC_eval_2@auc, 3)
round(tgb_LQH_eval_2@auc, 3)
round(tgb_LQ_eval_2@auc, 3)

round(tgb_AllFC_eval_3@auc, 3)
round(tgb_LQH_eval_3@auc, 3)
round(tgb_LQ_eval_3@auc, 3)

### mean
round(mean(c(tgb_AllFC_eval_1@auc, tgb_AllFC_eval_2@auc, tgb_AllFC_eval_3@auc)),3) # All_FC
round(mean(c(tgb_LQH_eval_1@auc, tgb_LQH_eval_2@auc, tgb_LQH_eval_3@auc)),3) # LQH
round(mean(c(tgb_LQ_eval_1@auc, tgb_LQ_eval_2@auc, tgb_LQ_eval_3@auc)),3) # LQ

# save external eval models as list and remove objects from workspace to reduce project weight.
ext_eval_maxent <- list(Maxent_bfile_allFC_1, Maxent_bfile_allFC_2, Maxent_bfile_allFC_3, Maxent_bmodel_allFC_1, Maxent_bmodel_allFC_2, Maxent_bmodel_allFC_3, Maxent_rndm_allFC_1, Maxent_rndm_allFC_2, Maxent_rndm_allFC_3, Maxent_tgb_allFC_1, Maxent_tgb_allFC_2, Maxent_tgb_allFC_3, Maxent_bfile_LQH_1, Maxent_bfile_LQH_2, Maxent_bfile_LQH_3, Maxent_bmodel_LQH_1, Maxent_bmodel_LQH_2, Maxent_bmodel_LQH_3, Maxent_rndm_LQH_1, Maxent_rndm_LQH_2, Maxent_rndm_LQH_3, Maxent_tgb_LQH_1, Maxent_tgb_LQH_2, Maxent_tgb_LQH_3, Maxent_bfile_LQ_1, Maxent_bfile_LQ_2, Maxent_bfile_LQ_3, Maxent_bmodel_LQ_1, Maxent_bmodel_LQ_2, Maxent_bmodel_LQ_3, Maxent_rndm_LQ_1, Maxent_rndm_LQ_2, Maxent_rndm_LQ_3, Maxent_tgb_LQ_1, Maxent_tgb_LQ_2, Maxent_tgb_LQ_3)

# save list
save(ext_eval_maxent, file="./outputs/models/ext_eval_maxent.RData")

# remove list and all models
# 
# rm(ext_eval_models)
# rm(Maxent_bfile_allFC_1, Maxent_bfile_allFC_2, Maxent_bfile_allFC_3, Maxent_bmodel_allFC_1, Maxent_bmodel_allFC_2, Maxent_bmodel_allFC_3, Maxent_rndm_allFC_1, Maxent_rndm_allFC_2, Maxent_rndm_allFC_3, Maxent_tgb_allFC_1, Maxent_tgb_allFC_2, Maxent_tgb_allFC_3, Maxent_bfile_LQH_1, Maxent_bfile_LQH_2, Maxent_bfile_LQH_3, Maxent_bmodel_LQH_1, Maxent_bmodel_LQH_2, Maxent_bmodel_LQH_3, Maxent_rndm_LQH_1, Maxent_rndm_LQH_2, Maxent_rndm_LQH_3, Maxent_tgb_LQH_1, Maxent_tgb_LQH_2, Maxent_tgb_LQH_3, Maxent_bfile_LQ_1, Maxent_bfile_LQ_2, Maxent_bfile_LQ_3, Maxent_bmodel_LQ_1, Maxent_bmodel_LQ_2, Maxent_bmodel_LQ_3, Maxent_rndm_LQ_1, Maxent_rndm_LQ_2, Maxent_rndm_LQ_3, Maxent_tgb_LQ_1, Maxent_tgb_LQ_2, Maxent_tgb_LQ_3)





##################################################################################################
######### External evaluation south ##############################################################


# load eval south sets
# Set 1
load("./outputs/processed_data/PA_PIA_eval_1_south.RData")
# Set 2
load("./outputs/processed_data/PA_PIA_eval_2_south.RData")
# Set 3
load("./outputs/processed_data/PA_PIA_eval_3_south.RData")


##### Random background

### Set_1
rndm_AllFC_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_rndm_allFC_1)

rndm_LQH_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_rndm_LQH_1)

rndm_LQ_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_rndm_LQ_1)


### Set_2
rndm_AllFC_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_rndm_allFC_2)

rndm_LQH_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_rndm_LQH_2)

rndm_LQ_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_rndm_LQ_2)


### Set_3
rndm_AllFC_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_rndm_allFC_3)

rndm_LQH_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_rndm_LQH_3)

rndm_LQ_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_rndm_LQ_3)


##### Bmodel background

### Set_1
bmodel_AllFC_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_1)

bmodel_LQH_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_1)

bmodel_LQ_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_1)


### Set_2
bmodel_AllFC_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_2)

bmodel_LQH_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_2)

bmodel_LQ_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_2)


### Set_3
bmodel_AllFC_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bmodel_allFC_3)

bmodel_LQH_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQH_3)

bmodel_LQ_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bmodel_LQ_3)


##### Bfile background

### Set_1
bfile_AllFC_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bfile_allFC_1)

bfile_LQH_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bfile_LQH_1)

bfile_LQ_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_bfile_LQ_1)


### Set_2
bfile_AllFC_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bfile_allFC_2)

bfile_LQH_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bfile_LQH_2)

bfile_LQ_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_bfile_LQ_2)


### Set_3
bfile_AllFC_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bfile_allFC_3)

bfile_LQH_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bfile_LQH_3)

bfile_LQ_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_bfile_LQ_3)


##### Tgb background
### Set_1
tgb_AllFC_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_tgb_allFC_1)

tgb_LQH_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_tgb_LQH_1)

tgb_LQ_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = Maxent_tgb_LQ_1)


### Set_2
tgb_AllFC_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_tgb_allFC_2)

tgb_LQH_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_tgb_LQH_2)

tgb_LQ_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = Maxent_tgb_LQ_2)


### Set_3
tgb_AllFC_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_tgb_allFC_3)

tgb_LQH_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_tgb_LQH_3)

tgb_LQ_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = Maxent_tgb_LQ_3)


#### results

### rndm
round(rndm_AllFC_eval_1_south@auc, 3)
round(rndm_LQH_eval_1_south@auc, 3)
round(rndm_LQ_eval_1_south@auc, 3)

round(rndm_AllFC_eval_2_south@auc, 3)
round(rndm_LQH_eval_2_south@auc, 3)
round(rndm_LQ_eval_2_south@auc, 3)

round(rndm_AllFC_eval_3_south@auc, 3)
round(rndm_LQH_eval_3_south@auc, 3)
round(rndm_LQ_eval_3_south@auc, 3)


### bmodel
round(bmodel_AllFC_eval_1_south@auc, 3)
round(bmodel_LQH_eval_1_south@auc, 3)
round(bmodel_LQ_eval_1_south@auc, 3)

round(bmodel_AllFC_eval_2_south@auc, 3)
round(bmodel_LQH_eval_2_south@auc, 3)
round(bmodel_LQ_eval_2_south@auc, 3)

round(bmodel_AllFC_eval_3_south@auc, 3)
round(bmodel_LQH_eval_3_south@auc, 3)
round(bmodel_LQ_eval_3_south@auc, 3)


### bfile
round(bfile_AllFC_eval_1_south@auc, 3)
round(bfile_LQH_eval_1_south@auc, 3)
round(bfile_LQ_eval_1_south@auc, 3)

round(bfile_AllFC_eval_2_south@auc, 3)
round(bfile_LQH_eval_2_south@auc, 3)
round(bfile_LQ_eval_2_south@auc, 3)

round(bfile_AllFC_eval_3_south@auc, 3)
round(bfile_LQH_eval_3_south@auc, 3)
round(bfile_LQ_eval_3_south@auc, 3)


### tgb
round(tgb_AllFC_eval_1_south@auc, 3)
round(tgb_LQH_eval_1_south@auc, 3)
round(tgb_LQ_eval_1_south@auc, 3)

round(tgb_AllFC_eval_2_south@auc, 3)
round(tgb_LQH_eval_2_south@auc, 3)
round(tgb_LQ_eval_2_south@auc, 3)

round(tgb_AllFC_eval_3_south@auc, 3)
round(tgb_LQH_eval_3_south@auc, 3)
round(tgb_LQ_eval_3_south@auc, 3)




##################################################################################################
######### Final (full data) models ###############################################################



##### Random background
Maxent_rndm_allFC_final <- maxent(x=SWD_rndm_fmodel[7:18], p=id_rndm_fmodel, args=c('-P','-J'), removeDuplicates=FALSE, path='outputs/models/rndm_bckgr/all_FC')
Maxent_rndm_LQ_final <- maxent(x=SWD_rndm_fmodel[7:18], p=id_rndm_fmodel,
                               args=c('-P', '-J', 
                                      'noautofeature', 
                                      'nohinge', 
                                      'noproduct'), removeDuplicates=FALSE, path='outputs/models/rndm_bckgr/LQ')
Maxent_rndm_LQH_final <- maxent(x=SWD_rndm_fmodel[7:18], p=id_rndm_fmodel,
                                args=c('-P', '-J',
                                       'noproduct'), removeDuplicates=FALSE, path='outputs/models/rndm_bckgr/LQH')


##### Bmodel background
Maxent_bmodel_allFC_final <- maxent(x=SWD_bmodel_fmodel[7:18], p=id_bmodel_fmodel, args=c('-P','-J'), removeDuplicates=FALSE, path='outputs/models/bmodel_bckgr/all_FC')
Maxent_bmodel_LQ_final <- maxent(x=SWD_bmodel_fmodel[7:18], p=id_bmodel_fmodel,
                                 args=c('-P', '-J', 
                                        'noautofeature', 
                                        'nohinge', 
                                        'noproduct'), removeDuplicates=FALSE, path='outputs/models/bmodel_bckgr/LQ')
Maxent_bmodel_LQH_final <- maxent(x=SWD_bmodel_fmodel[7:18], p=id_bmodel_fmodel,
                                  args=c('-P', '-J',
                                         'noproduct'), removeDuplicates=FALSE, path='outputs/models/bmodel_bckgr/LQH')


##### Bfile background
Maxent_bfile_allFC_final <- maxent(x=SWD_bfile_fmodel[7:18], p=id_bfile_fmodel, args=c('-P','-J'), removeDuplicates=FALSE, path='outputs/models/bfile_bckgr/all_FC')
Maxent_bfile_LQ_final <- maxent(x=SWD_bfile_fmodel[7:18], p=id_bfile_fmodel,
                                args=c('-P', '-J', 
                                       'noautofeature', 
                                       'nohinge', 
                                       'noproduct'), removeDuplicates=FALSE, path='outputs/models/bfile_bckgr/LQ')
Maxent_bfile_LQH_final <- maxent(x=SWD_bfile_fmodel[7:18], p=id_bfile_fmodel,
                                 args=c('-P', '-J',
                                        'noproduct'), removeDuplicates=FALSE, path='outputs/models/bfile_bckgr/LQH')


##### Tgb background
Maxent_tgb_allFC_final <- maxent(x=SWD_tgb_fmodel[7:18], p=id_tgb_fmodel, args=c('-P','-J'), removeDuplicates=FALSE, path='outputs/models/tgb_bckgr/all_FC')
Maxent_tgb_LQ_final <- maxent(x=SWD_tgb_fmodel[7:18], p=id_tgb_fmodel,
                              args=c('-P', '-J', 
                                     'noautofeature', 
                                     'nohinge', 
                                     'noproduct'), removeDuplicates=FALSE, path='outputs/models/tgb_bckgr/LQ')
Maxent_tgb_LQH_final <- maxent(x=SWD_tgb_fmodel[7:18], p=id_tgb_fmodel,
                               args=c('-P', '-J',
                                      'noproduct'), removeDuplicates=FALSE, path='outputs/models/tgb_bckgr/LQH')




##### Predictions

##### rndm
Maxent_rndm_allFC_pred <- predict(Maxent_rndm_allFC_final, vars_stack, path='./outputs/raster/Maxent_rndm_allFC.tif', overwrite=TRUE)
Maxent_rndm_LQ_pred <- predict(Maxent_rndm_LQ_final, vars_stack, path='./outputs/raster/Maxent_rndm_LQ.tif', overwrite=TRUE)
Maxent_rndm_LQH_pred <- predict(Maxent_rndm_LQH_final, vars_stack, path='./outputs/raster/Maxent_rndm_LQH.tif', overwrite=TRUE) # sometimes the function does not write the predicitn raster into file, so save manually:
# write predictions
writeRaster(Maxent_rndm_allFC_pred, './outputs/raster/Maxent_rndm_allFC.tif', overwrite=TRUE)
writeRaster(Maxent_rndm_LQ_pred, './outputs/raster/Maxent_rndm_LQ.tif', overwrite=TRUE)
writeRaster(Maxent_rndm_LQH_pred, './outputs/raster/Maxent_rndm_LQH.tif', overwrite=TRUE)


##### bmodel
Maxent_bmodel_allFC_pred <- predict(Maxent_bmodel_allFC_final, vars_stack, path='./outputs/raster/Maxent_bmodel_allFC.tif', overwrite=TRUE)
Maxent_bmodel_LQ_pred <- predict(Maxent_bmodel_LQ_final, vars_stack, path='./outputs/raster/Maxent_bmodel_LQ.tif', overwrite=TRUE)
Maxent_bmodel_LQH_pred <- predict(Maxent_bmodel_LQH_final, vars_stack, path='./outputs/raster/Maxent_bmodel_LQH.tif', overwrite=TRUE)
# write predictions
writeRaster(Maxent_bmodel_allFC_pred, './outputs/raster/Maxent_bmodel_allFC.tif', overwrite=TRUE)
writeRaster(Maxent_bmodel_LQ_pred, './outputs/raster/Maxent_bmodel_LQ.tif', overwrite=TRUE)
writeRaster(Maxent_bmodel_LQH_pred, './outputs/raster/Maxent_bmodel_LQH.tif', overwrite=TRUE)


##### bfile
Maxent_bfile_allFC_pred <- predict(Maxent_bfile_allFC_final, vars_stack, path='./outputs/raster/Maxent_bfile_allFC.tif', overwrite=TRUE)
Maxent_bfile_LQ_pred <- predict(Maxent_bfile_LQ_final, vars_stack, path='./outputs/raster/Maxent_bfile_LQ.tif', overwrite=TRUE)
Maxent_bfile_LQH_pred <- predict(Maxent_bfile_LQH_final, vars_stack, path='./outputs/raster/Maxent_bfile_LQH.tif', overwrite=TRUE)
# write predictions
writeRaster(Maxent_bfile_allFC_pred, './outputs/raster/Maxent_bfile_allFC.tif', overwrite=TRUE)
writeRaster(Maxent_bfile_LQ_pred, './outputs/raster/Maxent_bfile_LQ.tif', overwrite=TRUE)
writeRaster(Maxent_bfile_LQH_pred, './outputs/raster/Maxent_bfile_LQH.tif', overwrite=TRUE)


##### tgb
Maxent_tgb_allFC_pred <- predict(Maxent_tgb_allFC_final, vars_stack, path='./outputs/raster/Maxent_tgb_allFC.tif', overwrite=TRUE)
Maxent_tgb_LQ_pred <- predict(Maxent_tgb_LQ_final, vars_stack, path='./outputs/raster/Maxent_tgb_LQ.tif', overwrite=TRUE)
Maxent_tgb_LQH_pred <- predict(Maxent_tgb_LQH_final, vars_stack, path='./outputs/raster/Maxent_tgb_LQH.tif', overwrite=TRUE)
# write predictions
writeRaster(Maxent_tgb_allFC_pred, './outputs/raster/Maxent_tgb_allFC.tif', overwrite=TRUE)
writeRaster(Maxent_tgb_LQ_pred, './outputs/raster/Maxent_tgb_LQ.tif', overwrite=TRUE)
writeRaster(Maxent_tgb_LQH_pred, './outputs/raster/Maxent_tgb_LQH.tif', overwrite=TRUE)



#### remove LQ and LQH predictions from workspace
rm(Maxent_rndm_LQ_pred, Maxent_rndm_LQH_pred, Maxent_bmodel_LQ_pred, Maxent_bmodel_LQH_pred, Maxent_bfile_LQ_pred, Maxent_bfile_LQH_pred, Maxent_tgb_LQ_pred, Maxent_tgb_LQH_pred)



