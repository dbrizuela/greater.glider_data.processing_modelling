### Script to fit generalized linear models across the entire distribution range of the greater glider

# load required libraries
library(raster)
library(dismo)
# install.packages("gbm")
library(gbm)
library(rgdal)
library(tidyverse)
library(blockCV)
library(sp)
library(car)
library(effects)
library(gam)
library(dominanceanalysis)
library(effects)
library(caret)



##################################################################################################
######### Load occ data and background points ####################################################


### Load data if not already in workspace:

# ##### Sets for external evaluation
# 
# ### PA
# # Set 1
# load("./outputs/PA_fit_1.RData")
# # Set 2
# load("./outputs/PA_fit_2.RData")
# # Set 3
# load("./outputs/PA_fit_3.RData")
# 
# ### PA_PIA
# # Set 1
# load("./outputs/PA_PIA_fit_1.RData")
# # Set 2
# load("./outputs/PA_PIA_fit_2.RData")
# # Set 3
# load("./outputs/PA_PIA_fit_3.RData")
# 
# ### External evaluation sets
# # Set 1
# load("./outputs/PA_PIA_eval_1.RData")
# # Set 2
# load("./outputs/PA_PIA_eval_2.RData")
# # Set 3
# load("./outputs/PA_PIA_eval_3.RData")
# 
# 
# ### External evaluation sets Southern half
# # Set 1
# load("./outputs/PA_PIA_eval_1_south.RData")
# # Set 2
# load("./outputs/PA_PIA_eval_2_south.RData")
# # Set 3
# load("./outputs/PA_PIA_eval_3_south.RData")
# 
# 
# ### Final models sets
# 
# # PA 
# load("./outputs/PA_fmodel.RData")
# # PA_PIA
# load("./outputs/PA_PIA_fmodel.RData")
# 
# str(PA_fit_1)
# 
# # for some reason everytime I upload these again, PA column becomes factor, so change to numeric:
# PA_fit_1$PA <- as.numeric(levels(PA_fit_1$PA))[PA_fit_1$PA]
# PA_fit_2$PA <- as.numeric(levels(PA_fit_2$PA))[PA_fit_2$PA]
# PA_fit_3$PA <- as.numeric(levels(PA_fit_3$PA))[PA_fit_3$PA]
# 
# PA_PIA_fit_1$PA <- as.numeric(levels(PA_PIA_fit_1$PA))[PA_PIA_fit_1$PA]
# PA_PIA_fit_2$PA <- as.numeric(levels(PA_PIA_fit_2$PA))[PA_PIA_fit_2$PA]
# PA_PIA_fit_3$PA <- as.numeric(levels(PA_PIA_fit_3$PA))[PA_PIA_fit_3$PA]
# 
# PA_PIA_eval_1$PA <- as.numeric(levels(PA_PIA_eval_1$PA))[PA_PIA_eval_1$PA]
# PA_PIA_eval_2$PA <- as.numeric(levels(PA_PIA_eval_2$PA))[PA_PIA_eval_2$PA]
# PA_PIA_eval_3$PA <- as.numeric(levels(PA_PIA_eval_3$PA))[PA_PIA_eval_3$PA]
# 
# PA_PIA_eval_1_south$PA <- as.numeric(levels(PA_PIA_eval_1_south$PA))[PA_PIA_eval_1_south$PA]
# PA_PIA_eval_2_south$PA <- as.numeric(levels(PA_PIA_eval_2_south$PA))[PA_PIA_eval_2_south$PA]
# PA_PIA_eval_3_south$PA <- as.numeric(levels(PA_PIA_eval_3_south$PA))[PA_PIA_eval_3_south$PA]
# 
# # full datasets
# PA_fmodel$PA <- as.numeric(levels(PA_fmodel$PA))[PA_fmodel$PA]
# PA_PIA_fmodel$PA <- as.numeric(levels(PA_PIA_fmodel$PA))[PA_PIA_fmodel$PA]



##################################################################################################
######### Load and prepare env var  ##############################################################

### This is done in the 'data_processing' script as well:

# # All variables have been resampled (500 m), projected (epsg:3577) and masked to modelling extent
# vars_stack_all <- raster::stack(list.files(paste0(getwd(),"./spatial_data/variables_model"), pattern = '.asc', full.names=TRUE))
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
######### GLMs full data models ##################################################################

# In the case of GLMs, full data models are fitted first to get the formula of the backwards selected model.

##### PA
GLM_PA <- glm(PA ~ 
                poly(temp_warm,2)+ 
                poly(pp_annual,2)+ 
                poly(pp_driest,2)+ 
                poly(temp_season,2)+ 
                poly(pp_season,2)+ 
                poly(fPAR_mean,2)+ 
                poly(fPAR_var,2)+ 
                poly(GPP_var,2)+ 
                poly(mean_EDI,2)+ 
                poly(tsf,2)+ 
                poly(c_height,2)+ 
                poly(mean_xveg,2),
              family = binomial(link = "logit"), data = PA_fmodel[,c(2,7:18)])

summary(GLM_PA)
# Null deviance: 6945.6  on 5393  degrees of freedom
# Residual deviance: 5791.3  on 5369  degrees of freedom
# AIC: 5841.3


##### PIA
GLM_PIA <- glm(PA ~ 
                 poly(temp_warm,2)+ 
                 poly(pp_annual,2)+ 
                 poly(pp_driest,2)+ 
                 poly(temp_season,2)+ 
                 poly(pp_season,2)+ 
                 poly(fPAR_mean,2)+ 
                 poly(fPAR_var,2)+ 
                 poly(GPP_var,2)+ 
                 poly(mean_EDI,2)+ 
                 poly(tsf,2)+ 
                 poly(c_height,2)+ 
                 poly(mean_xveg,2),
               family = binomial(link = "logit"), data = PA_PIA_fmodel[,c(2,7:18)])

summary(GLM_PIA)
# Null deviance: 17025  on 14283  degrees of freedom
# Residual deviance: 13471  on 14259  degrees of freedom
# AIC: 13521


# scope for backwards selection
glm.scope <- list("temp_warm" = ~1 + temp_warm + poly(temp_warm,2),
                  "pp_annual" = ~1 + pp_annual + poly(pp_annual,2),
                  "pp_driest" = ~1 + pp_driest + poly(pp_driest,2),
                  "temp_season" = ~1 + temp_season + poly(temp_season,2),
                  "pp_season" = ~1 + pp_season + poly(pp_season,2),
                  "fPAR_mean" = ~1 + fPAR_mean + poly(fPAR_mean,2),
                  "fPAR_var" = ~1 + fPAR_var + poly(fPAR_var,2),
                  "GPP_var" = ~1 + GPP_var + poly(GPP_var,2),
                  "mean_EDI" = ~1 + mean_EDI + poly(mean_EDI,2),
                  "tsf" = ~1 + tsf + poly(tsf,2),
                  "c_height" = ~1 + c_height + poly(c_height,2),
                  "mean_xveg" = ~1 + mean_xveg + poly(mean_xveg,2))


##### PA back
GLM_PA_back <- step.Gam(GLM_PA, scope = glm.scope, direction= "backward", trace = T)

summary(GLM_PA_back)
# Call:
# glm(formula = PA ~ 
#       poly(temp_warm, 2) + 
#       poly(pp_annual, 2) + 
#       poly(pp_driest, 2) + 
#       poly(temp_season, 2) + 
#       poly(pp_season, 2) + 
#       poly(fPAR_mean, 2) + 
#       poly(fPAR_var, 2) + 
#       poly(GPP_var, 2) + 
#       mean_EDI + 
#       poly(tsf, 2) + 
#       poly(c_height, 2) + 
#       poly(mean_xveg, 2), family = binomial(link = "logit"), data = PA_fmodel[, c(2, 7:18)], trace = FALSE)

# Null deviance: 6945.6  on 5393  degrees of freedom
# Residual deviance: 5791.5  on 5370  degrees of freedom
# AIC: 5839.5



##### PIA back
GLM_PIA_back <- step.Gam(GLM_PIA, scope = glm.scope, direction= "backward", trace = T)

# save
# save(GLM_PIA_back, file = "./outputs/GLM_PIA_back.RData")

summary(GLM_PIA_back)
# Call:
# glm(formula = PA ~ 
#       poly(temp_warm, 2) + 
#       poly(pp_annual, 2) +  
#       pp_driest + 
#       poly(temp_season, 2) + 
#       poly(pp_season, 2) + 
#       poly(fPAR_mean, 2) + 
#       fPAR_var + 
#       poly(GPP_var, 2) + 
#       poly(mean_EDI, 2) + 
#       poly(tsf, 2) + 
#       poly(c_height, 2) + 
#       mean_xveg, family = binomial(link = "logit"), data = PA_PIA_fmodel[, c(2, 7:18)], trace = FALSE)

# Null deviance: 17025  on 14283  degrees of freedom
# Residual deviance: 13473  on 14262  degrees of freedom
# AIC: 13517


##### Prediction

#### PA
GLM_PA_back_pred <- raster::predict(vars_stack, GLM_PA_back, type="response", filename="./outputs/raster/GLM_PA_back.tif", overwrite=TRUE)

#### PIA
GLM_PIA_back_pred <- raster::predict(vars_stack, GLM_PIA_back, type="response", filename="./outputs/raster/GLM_PIA_back.tif", overwrite=TRUE)




##################################################################################################
######### create blocks for internal evaluation with blockCV #####################################

#### create spatial blocks

### PA
locations_PA <- PA_fmodel[,c(4:5)]
loc_PA_sdf <- SpatialPointsDataFrame(locations_PA[,c("x", "y")], PA_fmodel)
crs(loc_PA_sdf) <- (vars_stack$pp_annual@crs)


sb_PA <- spatialBlock(speciesData = loc_PA_sdf,
                      species = "PA",
                      rasterLayer = vars_stack$pp_annual,
                      theRange = 25000, # size of the blocks
                      k = 10,
                      selection = 'random',
                      maskBySpecies = T,
                      showBlocks = T,
                      biomod2Format = FALSE)

# save(sb_PA, file = "./outputs/sb_PA.RData")

### PA_PIA
locations_PIA <- PA_PIA_fmodel[,c(4:5)]
loc_PIA_sdf <- SpatialPointsDataFrame(locations_PIA[,c("x", "y")], PA_PIA_fmodel)
crs(loc_PIA_sdf) <- (vars_stack$pp_annual@crs)


sb_PIA <- spatialBlock(speciesData = loc_PIA_sdf,
                       species = "PA",
                       rasterLayer = vars_stack$pp_annual,
                       theRange = 25000, # size of the blocks
                       k = 10,
                       selection = 'random',
                       maskBySpecies = T,
                       showBlocks = T,
                       biomod2Format = FALSE)

# save(sb_PIA, file = "./outputs/sb_PIA.RData")




##################################################################################################
######### Internal evaluation ####################################################################

##### PA
folds <- sb_PA$folds

k <- 10
eval_GLM_PA_back_blockcv_test <- vector('list', k)
eval_GLM_PA_back_blockcv_train <- vector('list', k)

for (i in 1:k) {  # at each iteration i, we use fold i for testing
  # extract the training data for iteration i 
  trainSet <- unlist(folds[[i]][1]) # extract the training set indices
  testSet <- unlist(folds[[i]][2]) # extract the testing set indices
  
  # fitting a glm to training fold i
  # formulas chosen by backwards selection.
  GLM_PA_back_blockCV <- glm(PA ~ 
                               poly(temp_warm, 2) +
                               poly(pp_annual, 2) +
                               poly(pp_driest, 2) +
                               poly(temp_season, 2) +
                               poly(pp_season, 2) +
                               poly(fPAR_mean, 2) +
                               poly(fPAR_var, 2) +
                               poly(GPP_var, 2) +
                               mean_EDI +
                               poly(tsf, 2) +
                               poly(c_height, 2) +
                               poly(mean_xveg, 2), family = binomial(link = "logit"), data = PA_fmodel[, c(2, 7:18)][trainSet,])
  
  # extract the testing data for iteration i (presences and absences separately); standardised covariates for GLM:
  testdat <- PA_fmodel[,c(2, 7:18)][testSet,]
  test.pres <- testdat[testdat$PA==1,]
  test.abs  <- testdat[testdat$PA==0,]
  
  traindat <- PA_fmodel[,c(2, 7:18)][trainSet,]
  train.pres <- traindat[traindat$PA==1,]
  train.abs  <- traindat[traindat$PA==0,]
  
  # evaluate model with held-out data
  eval_GLM_PA_back_blockcv_train[[i]] <- evaluate(p=train.pres,a=train.abs,model=GLM_PA_back_blockCV, type="response")
  eval_GLM_PA_back_blockcv_test[[i]] <- evaluate(p=test.pres,a=test.abs,model=GLM_PA_back_blockCV, type="response")  
}



##### PIA
folds <- sb_PIA$folds

k <- 10
eval_GLM_PIA_back_blockcv_test <- vector('list', k)
eval_GLM_PIA_back_blockcv_train <- vector('list', k)

for (i in 1:k) {  # at each iteration i, we use fold i for testing
  # extract the training data for iteration i 
  trainSet <- unlist(folds[[i]][1]) # extract the training set indices
  testSet <- unlist(folds[[i]][2]) # extract the testing set indices
  
  # fitting a glm to training fold i
  # formulas chosen by backwards selection.
  
  GLM_PIA_back_blockCV <- glm(PA ~ 
                                poly(temp_warm, 2) +
                                poly(pp_annual, 2) +
                                pp_driest +
                                poly(temp_season, 2) +
                                poly(pp_season, 2) +
                                poly(fPAR_mean, 2) +
                                fPAR_var +
                                poly(GPP_var, 2) +
                                poly(mean_EDI, 2) +
                                poly(tsf, 2) +
                                poly(c_height, 2) +
                                mean_xveg, family = binomial(link = "logit"), data = PA_PIA_fmodel[, c(2, 7:18)][trainSet,])
  
  # extract the testing data for iteration i (presences and absences separately); standardised covariates for GLM:
  testdat <- PA_PIA_fmodel[,c(2, 7:18)][testSet,]
  test.pres <- testdat[testdat$PA==1,]
  test.abs  <- testdat[testdat$PA==0,]
  
  traindat <- PA_PIA_fmodel[,c(2, 7:18)][trainSet,]
  train.pres <- traindat[traindat$PA==1,]
  train.abs  <- traindat[traindat$PA==0,]
  
  # evaluate model with held-out data
  eval_GLM_PIA_back_blockcv_train[[i]] <- evaluate(p=train.pres,a=train.abs,model=GLM_PIA_back_blockCV, type="response")
  eval_GLM_PIA_back_blockcv_test[[i]] <- evaluate(p=test.pres,a=test.abs,model=GLM_PIA_back_blockCV, type="response")  
}


##### PA
eval_GLM_PA_back_blockcv_test <- sapply(eval_GLM_PA_back_blockcv_test, slot, 'auc')
round(mean(eval_GLM_PA_back_blockcv_test), 3)
round(sd(eval_GLM_PA_back_blockcv_test), 2)

##### PIA
eval_GLM_PIA_back_blockcv_test <- sapply(eval_GLM_PIA_back_blockcv_test, slot, 'auc')
round(mean(eval_GLM_PIA_back_blockcv_test), 3)
round(sd(eval_GLM_PIA_back_blockcv_test), 2)





##################################################################################################
######### External evaluation ####################################################################

# Fit one GLM per each external evaluation fitting set. All of them fitted with the formula obtained from backwards selection

##### PA

# Set_1
GLM_PA_back_1 <- glm(PA ~ 
                       poly(temp_warm, 2) +
                       poly(pp_annual, 2) +
                       poly(pp_driest, 2) +
                       poly(temp_season, 2) +
                       poly(pp_season, 2) +
                       poly(fPAR_mean, 2) +
                       poly(fPAR_var, 2) +
                       poly(GPP_var, 2) +
                       mean_EDI +
                       poly(tsf, 2) +
                       poly(c_height, 2) +
                       poly(mean_xveg, 2), family = binomial(link = "logit"), data = PA_fit_1[, c(2, 7:18)])


# Set_2
GLM_PA_back_2 <- glm(PA ~ 
                       poly(temp_warm, 2) +
                       poly(pp_annual, 2) +
                       poly(pp_driest, 2) +
                       poly(temp_season, 2) +
                       poly(pp_season, 2) +
                       poly(fPAR_mean, 2) +
                       poly(fPAR_var, 2) +
                       poly(GPP_var, 2) +
                       mean_EDI +
                       poly(tsf, 2) +
                       poly(c_height, 2) +
                       poly(mean_xveg, 2), family = binomial(link = "logit"), data = PA_fit_2[, c(2, 7:18)])


# Set_3
GLM_PA_back_3 <- glm(PA ~ 
                       poly(temp_warm, 2) +
                       poly(pp_annual, 2) +
                       poly(pp_driest, 2) +
                       poly(temp_season, 2) +
                       poly(pp_season, 2) +
                       poly(fPAR_mean, 2) +
                       poly(fPAR_var, 2) +
                       poly(GPP_var, 2) +
                       mean_EDI +
                       poly(tsf, 2) +
                       poly(c_height, 2) +
                       poly(mean_xveg, 2), family = binomial(link = "logit"), data = PA_fit_3[, c(2, 7:18)])



##### PIA

# Set_1
GLM_PIA_back_1 <- glm(PA ~ 
                        poly(temp_warm, 2) +
                        poly(pp_annual, 2) +
                        pp_driest +
                        poly(temp_season, 2) +
                        poly(pp_season, 2) +
                        poly(fPAR_mean, 2) +
                        fPAR_var +
                        poly(GPP_var, 2) +
                        poly(mean_EDI, 2) +
                        poly(tsf, 2) +
                        poly(c_height, 2) +
                        mean_xveg, family = binomial(link = "logit"), data = PA_PIA_fit_1[, c(2, 7:18)])

# Set_2
GLM_PIA_back_2 <- glm(PA ~ 
                        poly(temp_warm, 2) +
                        poly(pp_annual, 2) +
                        pp_driest +
                        poly(temp_season, 2) +
                        poly(pp_season, 2) +
                        poly(fPAR_mean, 2) +
                        fPAR_var +
                        poly(GPP_var, 2) +
                        poly(mean_EDI, 2) +
                        poly(tsf, 2) +
                        poly(c_height, 2) +
                        mean_xveg, family = binomial(link = "logit"), data = PA_PIA_fit_2[, c(2, 7:18)])

# Set_3
GLM_PIA_back_3 <- glm(PA ~ 
                        poly(temp_warm, 2) +
                        poly(pp_annual, 2) +
                        pp_driest +
                        poly(temp_season, 2) +
                        poly(pp_season, 2) +
                        poly(fPAR_mean, 2) +
                        fPAR_var +
                        poly(GPP_var, 2) +
                        poly(mean_EDI, 2) +
                        poly(tsf, 2) +
                        poly(c_height, 2) +
                        mean_xveg, family = binomial(link = "logit"), data = PA_PIA_fit_3[, c(2, 7:18)])


##### Evaluation

### PA models

# Set_1
GLM_PA_back_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = GLM_PA_back_1)

# Set_2
GLM_PA_back_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = GLM_PA_back_2)

# Set_3
GLM_PA_back_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = GLM_PA_back_3)


round(GLM_PA_back_eval_1@auc, 3)
round(GLM_PA_back_eval_2@auc, 3)
round(GLM_PA_back_eval_3@auc, 3)

# mean PA
round(mean(c(GLM_PA_back_eval_1@auc, GLM_PA_back_eval_2@auc, GLM_PA_back_eval_3@auc)),3)
round(sd(c(GLM_PA_back_eval_1@auc, GLM_PA_back_eval_2@auc, GLM_PA_back_eval_3@auc)),2)

### PA_PIA models

# Set_1
GLM_PIA_back_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = GLM_PIA_back_1)

# Set_2
GLM_PIA_back_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = GLM_PIA_back_2)

# Set_3
GLM_PIA_back_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = GLM_PIA_back_3)


round(GLM_PIA_back_eval_1@auc, 3)
round(GLM_PIA_back_eval_2@auc, 3)
round(GLM_PIA_back_eval_3@auc, 3)

# mean PIA
round(mean(c(GLM_PIA_back_eval_1@auc, GLM_PIA_back_eval_2@auc, GLM_PIA_back_eval_3@auc)),3)
round(sd(c(GLM_PIA_back_eval_1@auc, GLM_PIA_back_eval_2@auc, GLM_PIA_back_eval_3@auc)),2)



##################################################################################################
######### External evaluation south ##############################################################

##### PA models

# Set_1
GLM_PA_back_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = GLM_PA_back_1)

# Set_2
GLM_PA_back_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = GLM_PA_back_2)

# Set_3
GLM_PA_back_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = GLM_PA_back_3)


round(GLM_PA_back_eval_1_south@auc, 3)
round(GLM_PA_back_eval_2_south@auc, 3)
round(GLM_PA_back_eval_3_south@auc, 3)

# mean PA
round(mean(c(GLM_PA_back_eval_1_south@auc, GLM_PA_back_eval_2_south@auc, GLM_PA_back_eval_3_south@auc)),3)

round(sd(c(GLM_PA_back_eval_1_south@auc, GLM_PA_back_eval_2_south@auc, GLM_PA_back_eval_3_south@auc)),2) 


##### PA_PIA models

# Set_1
GLM_PIA_back_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = GLM_PIA_back_1)

# Set_2
GLM_PIA_back_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = GLM_PIA_back_2)

# Set_3
GLM_PIA_back_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = GLM_PIA_back_3)


round(GLM_PIA_back_eval_1_south@auc, 3)
round(GLM_PIA_back_eval_2_south@auc, 3)
round(GLM_PIA_back_eval_3_south@auc, 3)

# mean PIA
round(mean(c(GLM_PIA_back_eval_1_south@auc, GLM_PIA_back_eval_2_south@auc, GLM_PIA_back_eval_3_south@auc)),3)

round(sd(c(GLM_PIA_back_eval_1_south@auc, GLM_PIA_back_eval_2_south@auc, GLM_PIA_back_eval_3_south@auc)),2)


##################################################################################################
######### Variables importance ###################################################################

##### PA dominance analysis
dom_PA <-dominanceAnalysis(GLM_PA_back)
print(dom_PA)




##### PIA dominance analysisb
dom_PIA <-dominanceAnalysis(GLM_PIA_back)
print(dom_PIA)



##################################################################################################
######### Effects plots ##########################################################################

#### all response curves of PIA

plot(allEffects(GLM_PIA_back), type="response", main="", ylab="", ylim = c(0,1)) # plot with same y axis limits
dev.copy(png, "./outputs/plots/glm_allcurves.png", units = "in", width = 15, height = 9, res = 600) 
dev.off()

#### individual curves for the 4 most important variables according to domimance analysis:
# temp_season
# temp_warm
# mean_xveg   # this is different. In thesis model mean_xveg was 5th
# pp_annual


## temp_season
temp_season <-effect("temp_season",GLM_PIA_back)
png("./outputs/plots/glm_temp_season.png", units = "in", width = 7, height = 5, res = 600)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=1.6)))
plot(temp_season, main="", axes=list(y=list(type="response", lim = c(0,1), lab="", cex=1.5), x=list(cex=1.5)))
dev.off()
rm(temp_season)

## temp_warm
temp_warm <-effect("temp_warm",GLM_PIA_back)
png("./outputs/plots/glm_twarm.png", units = "in", width = 7, height = 5, res = 400)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=1.6)))
plot(temp_warm, main="", axes=list(y=list(type="response", lim = c(0,1), lab="", cex=1.5), x=list(cex=1.5)))
dev.off()
rm(temp_warm)

## mean_xveg
mean_xveg <-effect("mean_xveg",GLM_PIA_back)
png("./outputs/plots/glm_mxveg.png", units = "in", width = 7, height = 5, res = 400)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=1.6)))
plot(mean_xveg, main="", axes=list(y=list(type="response", lim = c(0,1), lab="", cex=1.5), x=list(cex=1.5)))
dev.off()
rm(mean_xveg)

## pp_annual
pp_annual <-effect("pp_annual",GLM_PIA_back)
png("./outputs/plots/glm_ppannual.png", units = "in", width = 7, height = 5, res = 400)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=1.6)))
plot(pp_annual, main="", axes=list(y=list(type="response", lim = c(0,1), lab="", cex=1.5), x=list(cex=1.5)))
dev.off()
rm(pp_annual)