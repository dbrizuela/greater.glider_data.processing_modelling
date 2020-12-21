
### Script to do post-processing of models' outputs. These include:
# Analysis of extrapolation (MESS maps)
# Maps of predictions differences
# Extraction of top 5% predicted values



##################################################################################################
######### Analysis of extrapolation (MESS maps) ##################################################

# Mess requires variables in raster stack and columns of variables in dataframes to be in the same order.

##### Maxent models

mess_rndm <- mess(x=vars_stack, v=SWD_rndm_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_rndm, paste0('./outputs/raster/mess_rndm/', names(mess_rndm), '.tif'),bylayer=TRUE, overwrite=T)

mess_bmodel <- mess(x=vars_stack, v=SWD_bmodel_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_bmodel, paste0('./outputs/raster/mess_bmodel/', names(mess_bmodel), '.tif'),bylayer=TRUE, overwrite=T)

mess_bfile <- mess(x=vars_stack, v=SWD_bfile_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_bfile, paste0('./outputs/raster/mess_bfile/', names(mess_bfile), '.tif'),bylayer=TRUE, overwrite=T)

mess_tgb <- mess(x=vars_stack, v=SWD_tgb_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_tgb, paste0('./outputs/raster/mess_tgb/', names(mess_tgb), '.tif'),bylayer=TRUE, overwrite=T)

##### PA
mess_PA <- mess(x=vars_stack, v=PA_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_PA, paste0('./outputs/raster/mess_PA/', names(mess_PA), '.tif'),bylayer=TRUE, overwrite=T)

##### PIA
mess_PIA <- mess(x=vars_stack, v=PA_PIA_fmodel[c(8:16,7,17,18)], full=T)
writeRaster(mess_PIA, paste0('./outputs/raster/mess_PIA/', names(mess_PIA), '.tif'),bylayer=TRUE, overwrite=T)


### get only celss with extrapolation
# define function
mess_extrapol <- function(layers, suffix="", file="") {
  
  output <- layers
  values(output)[values(output) >= 0] <- NA 
  names(output) <- paste0(names(output), suffix)
  writeRaster(output, paste0(file, names(output), '.tif'),bylayer=TRUE, overwrite=T)
  return(stack(output))
}

mess_rndm_0 <- mess_extrapol(mess_rndm, suffix="_rndm_0", file = './outputs/raster/mess_rndm/')
mess_bmodel_0 <- mess_extrapol(mess_bmodel, suffix="_bmodel_0", file = './outputs/raster/mess_bmodel/')
mess_bfile_0 <- mess_extrapol(mess_bfile, suffix="_bfile_0", file = './outputs/raster/mess_bfile/')
mess_tgb_0 <- mess_extrapol(mess_tgb, suffix="_tgb_0", file = './outputs/raster/mess_tgb/')
mess_PA_0 <- mess_extrapol(mess_tgb, suffix="_PA_0", file = './outputs/raster/mess_PIA')
mess_PIA_0 <- mess_extrapol(mess_tgb, suffix="_PIA_0", file = './outputs/raster/mess_PA/')

### Everytime I close R I have to load raster again!
# mess_rndm <- stack(list.files("./outputs/raster/mess_rndm", '.tif', full.names=TRUE))
# mess_bmodel <- stack(list.files("./outputs/raster/mess_bmodel", '.tif', full.names=TRUE))
# mess_bfile <- stack(list.files("./outputs/raster/mess_bfile", '.tif', full.names=TRUE))
# mess_tgb <- stack(list.files("./outputs/raster/mess_tgb", '.tif', full.names=TRUE))
# mess_PA <- stack(list.files("./outputs/raster/mess_PA", ".tif", full.names=T))
# mess_PIA <- stack(list.files("./outputs/raster/mess_PIA", ".tif", full.names=T))




##################################################################################################
######### Maps of predictions differences ########################################################

# This is only done for selected models

# function to rescale cell values between 0 and 1
rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

# rescale
Maxent_bmodel_01 <- rasterRescale(r = Maxent_bmodel)
Maxent_bfile_01 <- rasterRescale(r = Maxent_bfile)
BRT_PIA_pred_01 <- rasterRescale(r = BRT_PIA_pred)

### diff between bfile and bmodel
diff_bfile_bmodel <- Maxent_bfile_01-Maxent_bmodel_01
writeRaster(diff_bfile_bmodel, "./outputs/raster/diff_bfile_bmodel.tif")

# Diff bmodel to BRT
diff_bmodel_BRT <- Maxent_bmodel_01-BRT_PIA_pred_01
writeRaster(diff_bmodel_BRT, "./outputs/raster/diff_bmodel_BRT.tif")

# Diff bfile to BRT
diff_bfile_BRT <- Maxent_bfile_01-BRT_PIA_pred_01
writeRaster(diff_bfile_BRT, "./outputs/raster/diff_bfile_BRT.tif")


### correlation between pairs
pairs(stack(Maxent_bfile_01,Maxent_bmodel_01)) # diff_bfile_bmodel
pairs(stack(Maxent_bmodel_01,BRT_PIA_pred_01)) # diff_bmodel_BRT
pairs(stack(Maxent_bfile_01,BRT_PIA_pred_01)) # diff_bfile_BRT


### Correlation in south

# stack preds to clip to southern half
preds_south_stack <- stack(Maxent_bfile_01,Maxent_bmodel_01, BRT_PIA_pred_01)
preds_south_stack <- mask(preds_south_stack, vic_nsw)
preds_south_stack <- crop(preds_south_stack, vic_nsw)

names(preds_south_stack) <- c("bfile_allFC_south", "bmodel_allFC_south", "BRT_PIA_south")

pairs(subset(preds_south_stack, c(1,2))) # bfile_bmodel
pairs(subset(preds_south_stack, c(2,3))) # bmodel_BRT
pairs(subset(preds_south_stack, c(1,3))) # bfile_BRT

# write raster stack
writeRaster(preds_south_stack, paste0("./outputs/raster/", names(preds_south_stack), '.tif'),bylayer=TRUE)

# remove stack
rm(preds_south_stack)

#### diff raster won't be necessary in this project anymore, so delete to reduce project weight
# rm(diff_GLM_BRT, diff_bmodel_BRT, diff_bfile_BRT, diff_bfile_bmodel)




##################################################################################################
######### Top 5% predicted values ########################################################

# This is only done for selected models

# define function
top5_prediction <- function(layer, file="") {
  output <- layer
  quant <- as.vector(quantile(layer, probs = 0.95,na.rm=TRUE))
  values.top5 <- c(0, quant, 0,
                   quant, 1, 1)
  rcl.matrix <- matrix(values.top5, 
                       ncol=3, 
                       byrow=TRUE)
  reclass.raster <- raster::reclassify(output, rcl = rcl.matrix, filename=file, overwrite=T, include.lowest=FALSE, right=TRUE)
  return(reclass.raster)
  
}

bmodel_top5 <- top5_prediction(layer=Maxent_bmodel_allFC_pred, file = "./outputs/raster/bmodel_top5.tif")
bfile_top5 <- top5_prediction(layer=Maxent_bfile_allFC_pred, file = "./outputs/raster/bfile_top5.tif")
BRT_top5 <- top5_prediction(layer=BRT_PIA_pred, file = "./outputs/raster/BRT_PIA_top5.tif")


##### Calculate percentage agreement at the top 5%

# define function
percent_agree <- function(layer1, layer2) {
  
  total <- layer1+layer2
  NAs <- cellStats(total ==0, 'sum', na.rm=T)
  ones <- cellStats(total ==1, 'sum', na.rm=T)
  agree <- NAs <- cellStats(total ==2, 'sum', na.rm=T)
  
  # calculate percentage
  percent <- (agree*100)/(ones+agree)
  return(percent)
}

# agreement at top 5% bmodel and BRT_PIA
percent_agree(bmodel.95P, BRT.95P)
percent_agree(bfile.95P, BRT.95P)
percent_agree(bfile.95P, bmodel.95P)

