### Script to fit boosted regression trees models across the entire distribution range of the greater glider

### All required data is loaded in the script for fitting GLMs

#################### Internal evaluation with BlockCV ####################

##### PA

folds <- sb_PA$folds # Get the fold information

k <- 10
eval_BRT_PA_blockCV_train <- eval_BRT_PA_blockCV_test <- rep(NA, k)

for (i in 1:k) {  # at each iteration i, we use fold i for testing
  # extract the training data for iteration i 
  trainSet <- unlist(folds[[i]][1]) # extract the training set indices
  testSet <- unlist(folds[[i]][2])
  #fit models to training data:
  BRT_PA_blockCV <- gbm.step(data=PA_fmodel[trainSet,],
                             gbm.x = c(7:18),
                             gbm.y = 2,
                             family = "bernoulli",
                             tree.complexity = 5, 
                             learning.rate = 0.005,  
                             bag.fraction = 0.75)
  
  
  # extract the testing data for iteration i (presences and absences separately), different data formats for different methods:
  testdat <- PA_fmodel[testSet,]
  pvol.dat.test.pres <- testdat[testdat$PA==1,]
  pvol.dat.test.abs  <- testdat[testdat$PA==0,]
  
  num.train.pres <- table(PA_fmodel[trainSet,"PA"])[2]
  num.train.abs <- table(PA_fmodel[trainSet,"PA"])[1]
  num.test.pres <- nrow(pvol.dat.test.pres)
  num.test.abs <- nrow(pvol.dat.test.abs)
  
  #test data for brt
  Pres <-gbm::predict.gbm(BRT_PA_blockCV,newdata=pvol.dat.test.pres,n.trees=BRT_PA_blockCV$gbm.call$best.trees,type="response")
  Abs <-gbm::predict.gbm(BRT_PA_blockCV,newdata=pvol.dat.test.abs,n.trees=BRT_PA_blockCV$gbm.call$best.trees,type="response")
  mv <- wilcox.test(Pres,Abs)
  eval_BRT_PA_blockCV_test[i] <- as.numeric(mv$statistic) /(num.test.pres*num.test.abs) 
  
  #train data for brt - here use the fitted values, which are on the response scale
  mv <- wilcox.test(BRT_PA_blockCV$fitted[BRT_PA_blockCV$data$y==1],BRT_PA_blockCV$fitted[BRT_PA_blockCV$data$y==0])
  eval_BRT_PA_blockCV_train[i] <- as.numeric(mv$statistic) /(num.train.pres*num.train.abs)
  
}


# Get test AUC and sd
round(mean(eval_BRT_PA_blockCV_test),3)
round(sd(eval_BRT_PA_blockCV_test),2)




##### PIA

folds <- sb_PIA$folds

k <- 10
eval_BRT_PIA_blockCV_train <- eval_BRT_PIA_blockCV_test <- rep(NA, k)

for (i in 1:k) {
  trainSet <- unlist(folds[[i]][1]) 
  testSet <- unlist(folds[[i]][2])
  
  BRT_PIA_blockCV <- gbm.step(data=PA_PIA_fmodel[trainSet,],
                              gbm.x = c(7:18),
                              gbm.y = 2,
                              family = "bernoulli",
                              tree.complexity = 5, 
                              learning.rate = 0.05,  
                              bag.fraction = 0.75)
  
  testdat <- PA_PIA_fmodel[testSet,]
  pvol.dat.test.pres <- testdat[testdat$PA==1,]
  pvol.dat.test.abs  <- testdat[testdat$PA==0,]
  
  num.train.pres <- table(PA_PIA_fmodel[trainSet,"PA"])[2]
  num.train.abs <- table(PA_PIA_fmodel[trainSet,"PA"])[1]
  num.test.pres <- nrow(pvol.dat.test.pres)
  num.test.abs <- nrow(pvol.dat.test.abs)
  
  Pres <-gbm::predict.gbm(BRT_PIA_blockCV,newdata=pvol.dat.test.pres,n.trees=BRT_PIA_blockCV$gbm.call$best.trees,type="response")
  Abs <-gbm::predict.gbm(BRT_PIA_blockCV,newdata=pvol.dat.test.abs,n.trees=BRT_PIA_blockCV$gbm.call$best.trees,type="response")
  mv <- wilcox.test(Pres,Abs)
  eval_BRT_PIA_blockCV_test[i] <- as.numeric(mv$statistic) /(num.test.pres*num.test.abs) 
  
  mv <- wilcox.test(BRT_PIA_blockCV$fitted[BRT_PIA_blockCV$data$y==1],BRT_PIA_blockCV$fitted[BRT_PIA_blockCV$data$y==0])
  eval_BRT_PIA_blockCV_train[i] <- as.numeric(mv$statistic) /(num.train.pres*num.train.abs)
  
}


# Get test AUC and sd
round(mean(eval_BRT_PIA_blockCV_test),3)
round(sd(eval_BRT_PIA_blockCV_test),2)




#################### External evaluation ####################

##### PA

# Set_1
BRT_PA_1 <- gbm.step(data=PA_fit_1,
                   gbm.x = c(7:18),
                   gbm.y = 2,
                   family = "bernoulli",
                   tree.complexity = 5,
                   learning.rate = 0.005,
                   bag.fraction = 0.75)
# Set_2
BRT_PA_2 <- gbm.step(data=PA_fit_2,
                     gbm.x = c(7:18),
                     gbm.y = 2,
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate = 0.005,
                     bag.fraction = 0.75)
# Set_3
BRT_PA_3 <- gbm.step(data=PA_fit_3,
                     gbm.x = c(7:18),
                     gbm.y = 2,
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate = 0.005,
                     bag.fraction = 0.75)




# Set_1
BRT_PA_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = BRT_PA_1, n.trees=BRT_PA_1$gbm.call$best.trees)

# Set_2
BRT_PA_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = BRT_PA_2, n.trees=BRT_PA_2$gbm.call$best.trees)

# Set_3
BRT_PA_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = BRT_PA_3, n.trees=BRT_PA_3$gbm.call$best.trees)


round(BRT_PA_eval_1@auc, 3)
round(BRT_PA_eval_2@auc, 3)
round(BRT_PA_eval_3@auc, 3)

# mean PA
round(mean(c(BRT_PA_eval_1@auc, BRT_PA_eval_2@auc, BRT_PA_eval_3@auc)),3)
round(sd(c(BRT_PA_eval_1@auc, BRT_PA_eval_2@auc, BRT_PA_eval_3@auc)),2)



##### PIA

# Set_1
BRT_PIA_1 <- gbm.step(data=PA_PIA_fit_1,
                     gbm.x = c(7:18),
                     gbm.y = 2,
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate = 0.05,
                     bag.fraction = 0.75)
# Set_2
BRT_PIA_2 <- gbm.step(data=PA_PIA_fit_2,
                     gbm.x = c(7:18),
                     gbm.y = 2,
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate = 0.05,
                     bag.fraction = 0.75)
# Set_3
BRT_PIA_3 <- gbm.step(data=PA_PIA_fit_3,
                     gbm.x = c(7:18),
                     gbm.y = 2,
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate = 0.05,
                     bag.fraction = 0.75)


# Set_1
BRT_PIA_eval_1 <- dismo::evaluate(p = PA_PIA_eval_1[PA_PIA_eval_1$PA ==1, c(7:18)], a = PA_PIA_eval_1[PA_PIA_eval_1$PA ==0, c(7:18)], model = BRT_PIA_1, n.trees=BRT_PIA_1$gbm.call$best.trees)

# Set_2
BRT_PIA_eval_2 <- dismo::evaluate(p = PA_PIA_eval_2[PA_PIA_eval_2$PA ==1, c(7:18)], a = PA_PIA_eval_2[PA_PIA_eval_2$PA ==0, c(7:18)], model = BRT_PIA_2, n.trees=BRT_PIA_2$gbm.call$best.trees)

# Set_3
BRT_PIA_eval_3 <- dismo::evaluate(p = PA_PIA_eval_3[PA_PIA_eval_3$PA ==1, c(7:18)], a = PA_PIA_eval_3[PA_PIA_eval_3$PA ==0, c(7:18)], model = BRT_PIA_3, n.trees=BRT_PIA_3$gbm.call$best.trees)


round(BRT_PIA_eval_1@auc, 3)
round(BRT_PIA_eval_2@auc, 3)
round(BRT_PIA_eval_3@auc, 3)

# mean PIA
round(mean(c(BRT_PIA_eval_1@auc, BRT_PIA_eval_2@auc, BRT_PIA_eval_3@auc)),3)
round(sd(c(BRT_PIA_eval_1@auc, BRT_PIA_eval_2@auc, BRT_PIA_eval_3@auc)),2)



#################### southern half external evaluation ####################

##### PA

# Set_1
BRT_PA_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = BRT_PA_1, n.trees=BRT_PA_1$gbm.call$best.trees)

# Set_2
BRT_PA_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = BRT_PA_2, n.trees=BRT_PA_2$gbm.call$best.trees)

# Set_3
BRT_PA_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = BRT_PA_3, n.trees=BRT_PA_3$gbm.call$best.trees)


round(BRT_PA_eval_1_south@auc, 3)
round(BRT_PA_eval_2_south@auc, 3)
round(BRT_PA_eval_3_south@auc, 3)

# mean PA
round(mean(c(BRT_PA_eval_1_south@auc, BRT_PA_eval_2_south@auc, BRT_PA_eval_3_south@auc)),3)
round(sd(c(BRT_PA_eval_1_south@auc, BRT_PA_eval_2_south@auc, BRT_PA_eval_3_south@auc)),2)


##### PIA

# Set_1
BRT_PIA_eval_1_south <- dismo::evaluate(p = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==1, c(7:18)], a = PA_PIA_eval_1_south[PA_PIA_eval_1_south$PA ==0, c(7:18)], model = BRT_PIA_1, n.trees=BRT_PIA_1$gbm.call$best.trees)

# Set_2
BRT_PIA_eval_2_south <- dismo::evaluate(p = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==1, c(7:18)], a = PA_PIA_eval_2_south[PA_PIA_eval_2_south$PA ==0, c(7:18)], model = BRT_PIA_2, n.trees=BRT_PIA_2$gbm.call$best.trees)

# Set_3
BRT_PIA_eval_3_south <- dismo::evaluate(p = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==1, c(7:18)], a = PA_PIA_eval_3_south[PA_PIA_eval_3_south$PA ==0, c(7:18)], model = BRT_PIA_3, n.trees=BRT_PIA_3$gbm.call$best.trees)


round(BRT_PIA_eval_1_south@auc, 3)
round(BRT_PIA_eval_2_south@auc, 3)
round(BRT_PIA_eval_3_south@auc, 3)

# mean PIA
round(mean(c(BRT_PIA_eval_1_south@auc, BRT_PIA_eval_2_south@auc, BRT_PIA_eval_3_south@auc)),3)
round(sd(c(BRT_PIA_eval_1_south@auc, BRT_PIA_eval_2_south@auc, BRT_PIA_eval_3_south@auc)),2)




##################################################################################################
######### Final (full data) models ###############################################################

BRT_PA <- gbm.step(data=PA_fmodel,
                   gbm.x = c(7:18),
                   gbm.y = 2,
                   family = "bernoulli",
                   tree.complexity = 5,
                   learning.rate = 0.005,
                   bag.fraction = 0.75)



BRT_PIA <- gbm.step(data=PA_PIA_fmodel,
                    gbm.x = c(7:18),
                    gbm.y = 2,
                    family = "bernoulli",
                    tree.complexity = 5,
                    learning.rate = 0.05,
                    bag.fraction = 0.75)



# save 
# save(BRT_PIA, file = "./outputs/models/BRT_PIA.RData")

##### Predictions

#### PA
BRT_PA_pred <- predict(vars_stack, BRT_PA, n.trees=BRT_PA$gbm.call$best.trees, type="response")
writeRaster(BRT_PA_pred, "./outputs/raster/BRT_PA.tif", overwrite = T)

#### PA_PIA
BRT_PIA_pred <- predict(vars_stack, BRT_PIA, n.trees=BRT_PIA$gbm.call$best.trees, type="response")
writeRaster(BRT_PIA_pred, "./outputs/raster/BRT_PIA.tif", overwrite = T)



##################################################################################################
######### Variables' contributions, responses and interactions ###################################


##### Contributions

### PA
BRT_PA$contributions

### PIA
BRT_PIA$contributions



##### Interactions. PIA model only

BRT_PIA_int <- gbm.interactions(BRT_PIA)
BRT_PIA_int$rank.list
