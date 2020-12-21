#### Data processing

### This script contains the steps to process occurrence data.

#### Notes:

### All the obtained records were combined into a single dataframe and columns names standardized by authors. If steps for processing raw exactly as acquired from data sources, please contact Diego Brizuela-Torres at diego.brizuela.t@gmail.com

### sources of data are: 
# Victorian Biodiversity Atlas, Department of Environment, Land, Water and Planning, government of Victoria. Detailed species records, accessed Nov. 19, 2018
# BioNet, Office of Environment and Heritage, government of New South Wales. Species Sightings occurrence data, accessed Oct. 30, 2018.
# WildNet, Department of Environment and Science, government of Queensland. Species profile search, accessed Sep. 18, 2018.
# Atlas of Living Australia. Occurrence download, accessed Aug. 08, 2018.
# Teresa Eyre (Queensland herbarium)
# Brendan Wintle (University of Melbourne)

### Records from the New South Wales dataset 'State Forests Biodata' (SFB) are loaded separately because they are processed differently. 

# Load required libraries 
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(geosphere)
library(blockCV)
library(sf)
library(mapview)
library(spatialEco)

## remove except:  rm(list=(ls()[ls()!= c("occ_all", "bfile_raster")]))

##################################################################################################
#################### General data cleaning #######################################################


##### load data
load("./data/pvol_all.RData") # P. volans all records excepting 'State Forests Biodata'. 31965 records
load("./data/pvol_SFB.RData") # P. volans from 'State Forest Biodata'. 9427 records
load("./data/AM_all.RData") # All AM records excepting 'State Forest Biodata'. 326654 records
load("./data/AM_SFB.RData") # AM records from 'State Forests Biodata'. 22408 records 


##### change name of records listed as P. volans volans to just P. volans. P. volans minor remains. Although we modelled all records of single species (see main text), we thought it was not a bad idea to keep track of the records originally recorded as P.volans minor
pvol_all$species[pvol_all$species %in% "Petauroides volans volans"] <- "Petauroides volans"
pvol_SFB$species[pvol_SFB$species %in% "Petauroides volans volans"] <- "Petauroides volans"

##### merge all records to do a single coordinates projection and a single exclusion outside study area
occ_all <- dplyr::bind_rows(pvol_all, pvol_SFB, AM_all, AM_SFB) # bind_rows merges dataframes despite of having different columns. Need to have different columns because of the additional columns required for processing SFB.

rm(pvol_all, pvol_SFB, AM_all, AM_SFB)


##### Project coordinates to Australian Albers EPSG:3577

# need to round to detect duplicates from Australian states databases that were projected when made available in ALA.
occ_all$latitude.GDA94<-round(occ_all$latitude.GDA94,3) 
occ_all$longitude.GDA94<-round(occ_all$longitude.GDA94,3)

coordinates(occ_all)<- c(8,7)
proj4string(occ_all) <- CRS("+init=epsg:4283")
occ_all <- spTransform(occ_all, CRS("+init=epsg:3577"))


##### remove points outside study area 

# load study area polygon
bckgr_poly <- shapefile("./spatial_data/bckgr3577ok_proj.shp")
crs(bckgr_poly) <- ("+init=epsg:3577")

occ_all_bckgr <- occ_all[bckgr_poly,] # keep only records inside modelling extent 

occ_all <- as.data.frame(occ_all_bckgr) # n = 374741
rm(occ_all_bckgr)

##### remove by year and spatial accuracy

occ_all<- occ_all %>% dplyr::filter(!is.na(occ_all$year)) # remove records with no year. n = 370174
occ_all<- occ_all %>% dplyr::filter(!is.na(occ_all$accuracy)) # remove records with no accuracy. n = 365510 

# remove by year and accuracy pvol
occ_all <- occ_all[occ_all$year >= 1976,] # n = 357210   
occ_all <- occ_all[occ_all$accuracy <= 500,] # n = 321557

# separate presences and absences
pvol_pre <- occ_all %>% filter(species %in% c("Petauroides volans", "Petauroides volans minor") & presence.absence %in% 1) # n = 33841
pvol_abs <- occ_all %>% filter(species %in% c("Petauroides volans", "Petauroides volans minor") & presence.absence %in% 0) #  n = 1366
AM_all <- occ_all %>% filter(!species %in% c("Petauroides volans", "Petauroides volans minor")) # n = 286350


##### Duplicates removal by year and coordinates
pvol_pre <- pvol_pre %>% distinct(year, latitude.GDA94, longitude.GDA94, .keep_all = TRUE) # n = 14745      

pvol_abs <- pvol_abs %>% distinct(year, latitude.GDA94, longitude.GDA94, .keep_all = TRUE) # n = 1109 
# AM uses presence.absence as criterion to keep presences and absences of AM species as separate records
AM_all <- AM_all %>% distinct(presence.absence, year, latitude.GDA94, longitude.GDA94, .keep_all = TRUE) # n =  58073 


##### once duplicates have been removed, create Id column useful for later processes:

# pvol pre
pvol_pre$ID <-  1:nrow(pvol_pre)
pvol_pre$ID <- as.character(pvol_pre$ID)
pvol_pre$ID <- str_c("pvp", pvol_pre$ID, sep = "_", collapse = NULL) 

# pvol abs
pvol_abs$ID <-  1:nrow(pvol_abs)
pvol_abs$ID <- as.character(pvol_abs$ID)
pvol_abs$ID <- str_c("pva", pvol_abs$ID, sep = "_", collapse = NULL)

# arb.sp_all
AM_all$ID <-  1:nrow(AM_all)
AM_all$ID <- as.character(AM_all$ID)
AM_all$ID <- str_c("AM", AM_all$ID, sep = "_", collapse = NULL)


pvol_all <- rbind(pvol_pre, pvol_abs) # merge all pvol records. n = 15854
rm(pvol_pre, pvol_abs) # remove individual objects. 

pvol_all$presence.absence <- as.numeric(levels(pvol_all$presence.absence))[pvol_all$presence.absence] # make PA numeric

AM_all$presence.absence <- as.numeric(levels(AM_all$presence.absence))[AM_all$presence.absence] # make PA numeric

occ_all <- rbind(pvol_all, AM_all)

##################################################################################################
#################### 'State Forests Biodata' (SFB) processing ####################################


# separate SFB pvol records
pvol_all_SFB <- pvol_all %>% dplyr::filter(dataset.name %in% "State Forests Biodata") # n =  4819

# separate SFB AM records (AM_all has 73496)
AM_all_SFB <- AM_all %>% dplyr::filter(dataset.name %in% "State Forests Biodata") # n =  11103 



#################### Step 1: Clean data, only include non-incidental observations and for AM only spotlighting

summary(as.factor(AM_all_SFB$detection.method)) # Many of these are heard, but p.vol this is very few. 
summary(as.factor(pvol_all_SFB$detection.method))

### Get rid of all observations for other arboreal mammals that are not spotlighting.
dat_AM <- AM_all_SFB %>% dplyr::filter(detection.method %in% "spotlighting") # Select only rows det method = spotlighting, gets rid of most (10,000 records). n= 2003

dat_pvol <- pvol_all_SFB  %>% dplyr::filter(detection.method %in% "spotlighting") # n= 4490.   115 of these are absences. 

# Combine, so that can process all data together
dat.all<-rbind(dat_pvol, dat_AM) # n= 6493  


##### get pvol absences and do their matching separately based on LocationKey

pvol_abs_SFB <- dat_pvol %>% dplyr::filter(presence.absence %in% 0) # 115 absences. All absences have the 'SF' pattern in LocationKey
pvol_pre_SFB <- dat_pvol %>% dplyr::filter(presence.absence %in% 1) # 4375  

pvol_abs_SFB$LocationKey_match <- substr(pvol_abs_SFB$LocationKey, start = 1, stop = 7)
pvol_pre_SFB$LocationKey_match <- substr(pvol_pre_SFB$LocationKey, start = 1, stop = 7) # The records with the LP pattern will have an overly reduced LocationKey, but it doesn't matter as I'm matching the absences only to the SF pattern presence records (all absences have the SF pattern). This column can be deleted latter on.

##### also get AM records that match the pvol absences. Same situation: LocationKet_match column will make no sense for LP records, but it doesn't matter because all pvol absences have the SF pattern
AM_abs.match <- dat_AM # create additional object to do matching of explicit absences

AM_abs.match$LocationKey_match <- substr(AM_abs.match$LocationKey, start = 1, stop = 7)

# see how many matches are obtained with LocationKey
pre_match <- pvol_pre_SFB %>% dplyr::filter(LocationKey_match %in% pvol_abs_SFB$LocationKey_match)  # 18 matches

# see how many Location Keys matched
unique(pre_match$LocationKey_match)
# LocationKeys that matched:
# SF12620: 1 absence (2003-03-25), 8 presences. Out of these presences, 5 match well with the date of the absence. 3 of the presences do not match well, they were recorded about one month after. (2003-04-24). -3
# SF17871: 1 absence (1993-04-20), 6 presences match with the LocationKey but dates are not reasonable (1996-12-16, 1996-12-19, 1997-02-28). -6
# SF17870: 2 absences (1993-08-06), 2 presences match with LocationKey. 1 presences matches well with date, the other one doesn't (1993-02-22) -1
# SF73617: 1 absence (1999-11-25), 3 presences match with the LocationKey and exact date.

# Discard the matched presences with unreasonable dates.
pre_match$date <- as.character(pre_match$date) # apparently need to be transformed to character to be able to filter by date
pre_match <- pre_match %>% dplyr::filter(!date %in% c("2003-04-24", "1996-12-16", "1996-12-19", "1997-02-28", "1993-02-22")) # 8 presences

AM_abs.match <- AM_abs.match %>% dplyr::filter(LocationKey_match %in% pvol_abs_SFB$LocationKey_match) # 268 records

# remove dates from dat_AM_abs.match that do not match dates in pvol_abs_SFB. Just two dates
AM_abs.match$date <- as.character(AM_abs.match$date)
AM_abs.match <- AM_abs.match %>% dplyr::filter(!date %in% c("1993-02-12", "2003-04-24")) # 268-3 = 265 records

# All but 9 are AM explicit absences. These will be removed from evaluation dataset.
AM_abs.match <- AM_abs.match %>% dplyr::filter(presence.absence %in% 1) # 9 presences matching

# remove LocationKey_match column
pvol_abs_SFB <- pvol_abs_SFB[, -17]
pre_match <- pre_match[, -17]
AM_abs.match <- AM_abs.match[, -17]

# Merge pvol presences and AM records that matched with pvol abs. ## In these matches I'm incorporating presences and absences that matched reasonably as well as the 'inferred' absences that matched these records
pvol_abs_pre.AM.match <- rbind(pvol_abs_SFB, pre_match, AM_abs.match) # n=132 

pvol_abs_pre.AM.match$presence.absence[!pvol_abs_pre.AM.match$species %in% "Petauroides volans"] <- 0

pvol_abs_pre.AM.match$species[!pvol_abs_pre.AM.match$species %in% "Petauroides volans"] <- "inferred_absence"



#################### Step 2: Identify location keys that have > 2 records

# for the rest of the analysis, remove records that were already matched to the pvol_abs_SFB. should be 6493-132
dat.all <- dat.all %>% dplyr::filter(!ID %in% pvol_abs_pre.AM.match$ID) # 6361 records ok

summary(as.factor(dat.all$presence.absence)) # Mostly presences. Absences at this point are only from AM species, as I matched explicit pvol absences above

# remove 347 AM absences beforehand.  Should be  6361-347
dat.all <- dat.all %>% dplyr::filter(!presence.absence %in% 0) # 6014 

# Create new variable with first X characters of LocationKey
sur.type<-as.factor(substr(dat.all$LocationKey, start = 1, stop = 2)) # Figure out which rows are SF and which LP
LP.rows<-which(sur.type=="LP") # Get a list of rows that are LP surveys

dat.all$SurveyID<-substr(dat.all$LocationKey, start = 1, stop = 7) # Make a new variable SurveyID with the first 7 characters 
dat.all$SurveyID[LP.rows]<-substr(dat.all$LocationKey[LP.rows],start=1,stop=10) # Fix up LP rows so have 10 characters

# Add new variable for aggregating
dat.all$Count<-1 

survey.n<-aggregate(Count~SurveyID, data=dat.all, FUN=sum) # Count how many obs for each survey ID. 1871 surveys identified (after removing the surveys that had pvol absences and removing AM absences)

# Now get rid of surveys with 1 record only 
survey.keep<-survey.n[survey.n$Count>1,] # Gets rid of ~ half, left with DB:  1012         

# Reduce initial dataset so only have these surveys
dat.survey<-dat.all[dat.all$SurveyID %in% survey.keep$SurveyID,] #  5155



#################### Step 3: Check that have at least 2 records that are > 500m of each other (recorded within X days) 

# First make sure have proper dates column & add useful columns
dat.survey$date<-lubridate::ymd(dat.survey$date)

# Add new columns to store calculated data
dat.survey$min.timediff<-NULL # Create a new variable where put the minimum time difference to another obs in survey that is > 500m away
dat.survey$max.timediff<-NULL # Create a new variable where put the max time difference to another obs in survey that is > 500m away
dat.survey$min.distdiff<-NULL # Create a new variable where put the min distance to another obs in survey
dat.survey$max.distdiff<-NULL # Create a new variable where put the max distance to another obs in survey

# Convert to spatial for calculating distances. 
xy <- dat.survey[,c(14,15)]
dat.survey.spdf <- SpatialPointsDataFrame(coords = xy, data = dat.survey,proj4string=CRS("+init=epsg:3577"))

# Loop through each record and find the times and distances to all other records with the same surveyID. 
# Record the minimum time difference to a record > 500m away.
for(i in 1:length(dat.survey.spdf[,1])){
  p<-dat.survey.spdf[i,] # Get just the row interested in
  s<-p$SurveyID # extract survey name
  ID<-p$ID # extract record id
  dat<-dat.survey.spdf[dat.survey.spdf$SurveyID==s & dat.survey.spdf$ID != ID,] # Get the remaining records from that survey
  distances <- gDistance(p,dat,byid=TRUE) #Calculate the distances from the focal  obs to all other records in the survey
  times<-as.numeric(abs(difftime(p$date, dat$date, unit="days"))) # Get the absolute time difference to all other records in the survey 
  
  # # Now want to work out if more than 2 records > 500m apart within 1 week
  # Find the minimum time to a point > 500m away
  dist500<-which(distances>500) # which of the listed distances are >500
  dat.survey.spdf$min.timediff[i]<-min(times[dist500]) # Select the minimum time difference, using only the points that are > 500m away 
  dat.survey.spdf$max.timediff[i]<-max(times[dist500])
  dat.survey.spdf$min.distdiff[i]<-min(distances) # Save the min and max distances to other records in the survay
  dat.survey.spdf$max.distdiff[i]<-max(distances)
  dat.survey.spdf$n.survey[i]<-length(dat)+1 # Save the number of records in the survey
}

# Note that in the loop above, if no obs > 500m from point get 'Inf' for the recorded time diff.


# Add a column to make P/A based on species and type of obs 
dat.survey.spdf$PA<-0 # Make all absent
dat.survey.spdf$PA[dat.survey.spdf$species =="Petauroides volans" & dat.survey.spdf$presence.absence == 1]<-1 # Convert to present if p.vol and recorded as present in dataset

# Get rid of records that do not have another obs within 7 days that is >500m away
dat.survey.7days<-dat.survey.spdf[dat.survey.spdf$min.timediff <=7,] # Could check how impacts results if allow longer time period. 3856      

spplot(dat.survey.7days, "PA", alpha=0.05)

# check surveys
dat.survey.7days <- as.data.frame(dat.survey.7days) # 3856 records

## sustitue 'presence.absence' column for 'PA' and use 'PA' in all subsequent sets
dat.survey.7days$presence.absence <- dat.survey.7days$PA

#### merge dat.survey.7days with the previously subset pvol absencs and their matches.
SFB_processed_all <- dplyr::bind_rows(dat.survey.7days, pvol_abs_pre.AM.match) # n = 3988   
SFB_processed_all$species[!SFB_processed_all$species %in% "Petauroides volans"] <- "inferred_absence"

# remove columns that won't be necessary anymore
SFB_processed_all <- dplyr:: select(SFB_processed_all,  -c("SightingKey", "LocationKey", "Description", "SightingNotes", "SurveyID", "Count", "min.timediff", "max.timediff", "min.distdiff", "max.distdiff", "n.survey", "PA", "longitude.GDA94.1", "latitude.GDA94.1"))

# remove unnecesary objects
rm(AM_abs.match, pvol_abs_SFB, pvol_pre_SFB, pre_match, dat.all, dat.survey, dat.survey.7days, dat.survey.spdf, distances, p, pvol_abs_pre.AM.match, survey.keep, survey.n, xy)



##################################################################################################
#################### Obatin the rest of the PA survey data #######################################

##### PA data

# select surveys where all data is contained in a single dataset name
pvol_PAsurveys <- pvol_all %>% dplyr::filter(dataset.name %in% c("Arboreal mammal project 2015", 
                                                                 "Estimating the density of Greater Gliders in the Hume, Port Phillip and Gippsland Regions of Victoria", 
                                                                 "B.Wintle_Eden_single.visit", 
                                                                 "B.Wintle_Eden_repeated.visits", 
                                                                 "SWB", 
                                                                 "GEN", 
                                                                 "NQ", 
                                                                 "BMRG GG", 
                                                                 "DAWSON", 
                                                                 "SBB", 
                                                                 "EPASBB", 
                                                                 "CRA",
                                                                 "SPG",
                                                                 "TE")) # n = 1399

# select and merge surveys that are divided among multiple dataset names
pvol_PAsurveys <- rbind(pvol_PAsurveys,(dplyr::filter(pvol_all, grepl(c('MIL|BUN0|BAR1|BAR0|MH|SER0|STM|FEEDTR'),dataset.name)))) # n = 1467. This is all the data with PA explicit absences (plus the few explicit absences from SFB).

# remove non-spotlighting records
pvol_PAsurveys <- pvol_PAsurveys[pvol_PAsurveys$detection.method %in% "spotlighting", ] # n = 1466

# get names from PA_surveys
pvol_PAsurveys_names <- unique(pvol_PAsurveys$dataset.name) # 81



##################################################################################################
#################### Create inferred absences and match with presences ###########################

#get data set names from pvol_pre
pvol_dset.name <- unique(pvol_all$dataset.name) # 364 datasets name

# remove from datasets names the names of the PA surveys. 
pvol_dset.name <- pvol_dset.name[!pvol_dset.name %in% pvol_PAsurveys_names] # n = 283

pvol_dset.name <- pvol_dset.name[!pvol_dset.name %in% "State Forests Biodata"] # also remove from SFB.  n = 282

# AM records that come from the same surveys as pvol presences
AM_dset_pvol <- AM_all[AM_all$dataset.name %in% pvol_dset.name,] # n = 3568

# keep only presences (to consider only these as inferred absences)
AM_abs.dsetname <- AM_all[AM_all$presence.absence == 0,] # 515 AM absences
AM_abs.dsetname <- unique(AM_abs.dsetname$dataset.name) # n = 14 dataset names that registered AM absences


# remove AM records from datasets that also registered absences
# Here, it could be the case that some of the 'datasets' that include AM absences are also a compilation of many surveys. Only 14 datasets that include AM absences. Feasible to examine them individually and see which of these datasets seem to be a single survey (and exclude them, as they were registering absences), and which seem to be a compilation of different surveys. In this case, I'll keep the AM records to get the inferred absences and just remove the AM absences

# AM_abs.dsetname: 

# [1] "Arboreal mammal project 2015"              "Data Priorities Fauna Survey"             
# [3] "LNE Summer Survey 1997/98"                 "Lower North East CRA Survey"              
# [5] "OEH Data from Scientific Licences dataset" "OEH Default Sightings"                    
# [7] NA                                          "Koorawatha Area Survey"                   
# [9] "Misc AVW surveys from the 2000s"           "Lane Cove Fauna Survey"                   
# [11] "Sydney Region Fauna Surveys"               "63"                                       
# [13] "Pilliga NR Community Survey"               "State Forests Biodata" 

# [1] "Arboreal mammal project 2015": already discarded as includes pvol absences **drop**

# [2] "Data Priorities Fauna Survey": Records from many years (2003 to 2011), so is compilatory, however there's ony one AM absence. so apparently in most of the compiled survyes absences are not reported, so I'll keep it and remove AM absence **keep**

# [3] "LNE Summer Survey 1997/98": Seems compilatory as well. surveys in the two years (1997 and 1998) are close in dates to each other, so seems like this was a single 'surveying project'. There's one AM absence (the only absence from that dataset name), so seems like they were not really recording absences constantly. I'll keep it and remove AM absence   **keep**

# [4] "Lower North East CRA Survey". Exact same situation as above, excepting the years (1996 to 1998)  **keep** 

# [5] "OEH Data from Scientific Licences dataset": Compilatory from many years (1977-2018). just two AM absences, so seems like the majority didn't record/upload absences. Keep and just remove AM absences   **keep**

# [6] "OEH Default Sightings" Same situation as above. 3 AM absences, remove them.   **keep**

# [7] NA: NAs come from VBA. Won't use them, although they wouldn't match because there are no pvol datasets with NA   **drop**

# [8] "Koorawatha Area Survey" Single survey and apparently recording absences. No pvol match anyway. **drop**

# [9] "Misc AVW surveys from the 2000s": compiltaory (2000 to 2009). Only 3 AM absences. Drop these and keep.    **keep** 

# [10] "Lane Cove Fauna Survey" Surveys span from Feb to May 2004. one AM absence recorded, but seems reasonable that is a single surveying project. Only 2 pvol presences to match anyway.     **drop**

# [11] "Sydney Region Fauna Surveys". Compilatory (2002 to 2017). Only 1 AM absence. Drop this and keep.    **keep** 

# [12] "63": from VBA. Apparently single 'surveying project' only recorded Trichosurus cunningham and there's one absence. No matches in pvol anyway.  **drop**

# [13] "Pilliga NR Community Survey" Apparently single survey with AM absence recorded. No pvol matches anyway.     **drop

# keep in list only those datasets that recorded AM absences to then exclude them from AM_dset_pvol i.e. to don't use them to infer absences
AM_abs.dsetname <- AM_abs.dsetname[c(1,7,8,10,12,13)]


# remove from AM_dset_pvol the AM records from datasets that recorded AM absences.
AM_dset_pvol <- AM_dset_pvol[!AM_dset_pvol$dataset.name %in% AM_abs.dsetname,] # n = 3567 # just one dropped

# remove individual AM absences (12 absences)
AM_dset_pvol <- AM_dset_pvol[AM_dset_pvol$presence.absence %in% 1,] # n = 3567-12 = 3555

AM_dset_pvol_dset.name <- unique(AM_dset_pvol$dataset.name)

# pvol presence from which I'll inferr absences (those with the same dataset names as AM_dset_pvol)

pvol_inf <- pvol_all[pvol_all$dataset.name %in% AM_dset_pvol_dset.name ,] # n = 7629


# as a last attempt to remove incidental records, remove all records that have the word 'incidental' in the dataset name, from both AM and pvol before matching (310 records)

# column has to be factor to use the grepl function
pvol_inf <- dplyr::filter(pvol_inf, !(grepl('incidental|Incidental|unknown', pvol_inf$dataset.name))) # n = 7059
AM_dset_pvol <- dplyr::filter(AM_dset_pvol, !(grepl('incidental|Incidental|unknown', AM_dset_pvol$dataset.name))) # n = 3084


##### Spatial process of excluding AM presences less than 500 m from pvol pre

# need to use as spatial points
coordinates(pvol_inf) <- c(14,15)
crs(pvol_inf) <- ("+init=epsg:3577")

coordinates(AM_dset_pvol) <- c(14,15)
crs(AM_dset_pvol) <- ("+init=epsg:3577")

# create buffer around pvol presences
pvol_pre_buff <- buffer(pvol_inf, width=500, dissolve=TRUE)

# save shp just to make sure it's 500 radius (500m from the point outwards)
# shapefile(pvol_pre_buff, "./outputs/pvol_pre_buff.shp") # i checked in ArcMap and it is indeed 500 m radius, so OK!

# get AM records that that overlapp buffers polygons
AM_dset_pvol_buff <- over(AM_dset_pvol, pvol_pre_buff, returnList = FALSE)
AM_dset_pvol_buff <- as.data.frame(AM_dset_pvol_buff) 

# transform AM records contained in pvol datasets back to data frame
AM_dset_pvol <- as.data.frame(AM_dset_pvol)

# merge AM_dset_pvol  with column of id. records with NA in ts_dset_pvol_buff mean no overlapping
AM_buff_id <- cbind(AM_dset_pvol, AM_dset_pvol_buff)
AM_buff_id <- as.data.frame(AM_buff_id)

# get records with no overlapp
inf_abs <- AM_buff_id[is.na(AM_buff_id$AM_dset_pvol_buff),] # with removing by detection method: n =  2373.  

# remove column that had the buffer Id info
inf_abs <- inf_abs[,-17]

# change name to inf abs and make 0s for PA
inf_abs$species <- "inferred_absence"
inf_abs$presence.absence <- 0

pvol_inf <- as.data.frame(pvol_inf)

########## PA and PIA datasets. These datasets are yet to be combined to then be subset spatially into training and evaluation sets. Then, fitting sets will have records from evaluation blocks removed.

##### PIA dataset:
PIA <- rbind(inf_abs, pvol_inf) # n = 9432
# remove SFB columns from PIA
PIA <- PIA %>% dplyr::select(-c("SightingKey", "LocationKey", "Description", "SightingNotes"))

##### PA dataset. SFB_processed_all includes the inferred absences from AM records from SFB, but is effectively treated as PA:
# first remove SFB columns from pvol_PAsurveys
pvol_PAsurveys <- pvol_PAsurveys %>% dplyr::select(-c("SightingKey", "LocationKey", "Description", "SightingNotes"))

PA <- rbind(pvol_PAsurveys, SFB_processed_all) # n = 5454

##### PA + PIA dataset. This is the one that will be subset spatially into train and eval
PA_PIA <- rbind(PA, PIA) # n = 14886

summary(as.factor(PA_PIA$detection.method))

# remove records from remains identification and scat identification
PA_PIA <- dplyr::filter(PA_PIA, !detection.method %in% c("remains identification", "scat identification")) # n = 14665




##################################################################################################
#################### Create PB sets #############################################################

##### fitting presences. One records per 500m grid cell

pvol_all <- dplyr::select(pvol_all, -c("count", "date", "accuracy", "detection.method", "SightingKey", "LocationKey", "Description", "SightingNotes")) # only keepn necessary columns at this point

# subset presences
pvol_pre <- pvol_all[pvol_all$presence.absence == 1, ] # n = 14745

# load raster mask
mask500m <- raster("./spatial_data/variables_model/anuc01_ok.asc")
crs(mask500m) <- ("+init=epsg:3577")

# keep one presence per cell
pvol.cell500 <- cellFromXY(mask500m, pvol_pre[6:7])  
pvol.cell500.dups <- which(duplicated(pvol.cell500)) 
pvol_pre500 <- pvol_pre[-c(pvol.cell500.dups),]  # 9469 records

rm(pvol.cell500, pvol.cell500.dups)



##### random background points
points_rndm <- as.data.frame(dismo::randomPoints(mask500m, 180000))


##### Create TGB

occ_all <- dplyr::select(occ_all, -c("count", "date", "accuracy", "detection.method", "SightingKey", "LocationKey", "Description", "SightingNotes"))

# keep one occurrence per cell
occ_all.cell500 <- cellFromXY(mask500m, occ_all[6:7])  
occ_all.cell500.dups <- which(duplicated(occ_all.cell500)) 
tgb <- occ_all[-c(occ_all.cell500.dups),] # 42171 records

rm(occ_all.cell500, occ_all.cell500.dups)



##### Create bmodel background points

# load bmodel (see main text) raster. For details on this raster creation please contact Diego Brizuela-Torres at diego.brizuela.t@gmail.com 
bmodel_raster <- raster("./spatial_data/pred_roads_lcover.tif")
crs(bmodel_raster) <- ("+init=epsg:3577")

# create 180000 points. Create more than 100000 (final amount of backgroun points to use) because I'll loose more points along the rest of the data processing steps. Once all these steps are done, I'll sandomnly subset the final 100000 points.
points_bmodel <- as.data.frame(xyFromCell(bmodel_raster, sample(which(!is.na(values(bmodel_raster))), size=180000, prob=values(bmodel_raster)[!is.na(values(bmodel_raster))],replace=FALSE)))



##### Create bfile points

# create bfile raster

library(MASS)
#points to create bfile based on all unique locations of occurrences
points <- occ_all[,c(6,7)] # start from set of unique occurences I used for tgb

#rasterize points
points.ras <- rasterize(points, mask500m, 1)
rm(points)

# two-dimensional kernel density estimate. in thia case, it considers absences too. So it's a general representation of all sampling effort, not only of presences density.
points_1 <- which(values(points.ras) == 1)
points.locs <- coordinates(points.ras)[points_1, ]

dens <- kde2d(points.locs[,1], points.locs[,2], n = c(nrow(points.ras), ncol(points.ras)))
dens.ras <- raster(dens)

# standardize from 0 to 1
dens.ras_std <- raster.transformation(dens.ras, trans = "norm", smin = 0, smax = 1)

# crop raster to polygon
bfile_raster <- mask(dens.ras_std, bckgr_poly)

rm(dens.ras_std, dens, dens.ras, points_1, points.locs) # won't need these anymore

plot(bfile_raster, main="Bias file") # final bfile raster

writeRaster(bfile_raster, "./spatial_data/bfile_raster.tif") # advisable to save raster, as it takes a long time to create

##### create 180000 points
points_bfile <- as.data.frame(xyFromCell(bfile_raster, sample(which(!is.na(values(bfile_raster))), size=180000, prob=values(bfile_raster)[!is.na(values(bfile_raster))],replace=FALSE)))




##################################################################################################
#################### Assign TSF to background sets with no year (rndm, bmodel, bfile) ############

# steps:
# 1. load blocks: study area has been divided in 10 rows of equal size to get years of occurrence data and assign these to background points in each block. These blocks were created in ArcMap
# 2. separate occurrences and background points into blocks. For occurrences I'm using occ_all, which is oll occurrences (after removing by year and accuracy) excluding duplicates. it aims to represent total surveying effort.
# 3. sample years from occurrences and assign them to background points in the same block.
# 4. assign TSF

##### 1. load spatial blocks

fire_blocks.list <- list.files(paste0(getwd(),"./spatial_data/fire_blocks"), pattern = "\\.shp$", full.names=TRUE)

fire_blocks <- lapply(fire_blocks.list, raster::shapefile)

for(i in 1:length(fire_blocks)){
  raster::crs(fire_blocks[[i]]) <- ("+init=epsg:3577")
}



##### 2. separate occurrences by blocks

# create list
occ_block <- list()
length(occ_block) <- 10

# occ_all into spatial points
coordinates(occ_all) <- c(6,7)
crs(occ_all) <- ("+init=epsg:3577")

# separate by blocks
for(i in 1:length(fire_blocks)){
  occ_block[[i]] <- occ_all[fire_blocks[[i]],]
}


##### 3. sample years from occurrences and assign them to background points in the same block.

### random points
points_rndm$year <- NA # create year column
points_rndm$bckgr <- "rndm" # id column to identify the sets of backgr points after I merge themm all for following steps

coordinates(points_rndm) <- c(1,2)
crs(points_rndm) <- ("+init=epsg:3577")

rndm_block <- list()
length(rndm_block) <- 10

for(i in 1:length(fire_blocks)){
  rndm_block[[i]] <- points_rndm[fire_blocks[[i]],]
  rndm_block[[i]] <- as.data.frame(rndm_block[[i]])
  
  rndm_block[[i]]$year <- sample(occ_block[[i]]$year, size = nrow(rndm_block[[i]]), replace = TRUE) # sample points from occurrences and assign them to background points in the same block
}


### bmodel points
# create year column
points_bmodel$year <- NA
points_bmodel$bckgr <- "bmodel" 

coordinates(points_bmodel) <- c(1,2)
crs(points_bmodel) <- ("+init=epsg:3577")

bmodel_block <- list()
length(bmodel_block) <- 10

for(i in 1:length(fire_blocks)){
  bmodel_block[[i]] <- points_bmodel[fire_blocks[[i]],]
  bmodel_block[[i]] <- as.data.frame(bmodel_block[[i]])
  
  bmodel_block[[i]]$year <- sample(occ_block[[i]]$year, size = nrow(bmodel_block[[i]]), replace = TRUE) 
}


### bfile points
# create year column
points_bfile$year <- NA
points_bfile$bckgr <- "bfile" 

coordinates(points_bfile) <- c(1,2)
crs(points_bfile) <- ("+init=epsg:3577")

bfile_block <- list()
length(bfile_block) <- 10

for(i in 1:length(fire_blocks)){
  bfile_block[[i]] <- points_bfile[fire_blocks[[i]],]
  bfile_block[[i]] <- as.data.frame(bfile_block[[i]])
  
  bfile_block[[i]]$year <- sample(occ_block[[i]]$year, size = nrow(bfile_block[[i]]), replace = TRUE) 
}




#### 4. Assign TSF

# combine rndm, bmodel and bfiles to assign tsf in a single run
bckgr_all_tsf <- do.call("rbind", c(rndm_block, bmodel_block, bfile_block))

#Add empty column for tsf values
bckgr_all_tsf$tsf<-NA

# load fire polygons
fire_poly <- shapefile("./spatial_data/WF_allyears_bckgr.shp")
crs(fire_poly) <- ("+init=epsg:3577")

valNF <- 117    # set value for oldest fires in the study area.

# spatial points
coordinates(bckgr_all_tsf) <- c(1,2)
crs(bckgr_all_tsf) <- ("+init=epsg:3577")


## Extract fires
fires.ext_bckgr <- over(bckgr_all_tsf, fire_poly, returnList = TRUE) # extract all fires that overlap with points
sites_fires_bckgr <- sapply(fires.ext_bckgr, nrow) # Check how many fires overlap with each point, use these to set tsf to valNF if is 0


for(i in 1:length(sites_fires_bckgr)){ 
  if(sites_fires_bckgr[[i]]==0){ # If no recorded fires at location, then set to valNF
    bckgr_all_tsf$tsf[i]<-valNF
  }else{
    tsf<-bckgr_all_tsf$year[i] - fires.ext_bckgr[[i]]$FireYr # Otherwise calculate the interval from year of observation to year of each fire
    tsf[tsf<0]<-NA # if fire occured after observation (i.e. value is negative) then set to NA so not selected
    tsf<-c(tsf,valNF) # Add the valNF here to avoid errors when the only recorded fires occur after the observation
    bckgr_all_tsf$tsf[i]<-min(tsf,na.rm=T) # Select the minimum value and save as tsf
  }
}


# separate info rndm and bias background
bckgr_all_tsf <- as.data.frame(bckgr_all_tsf) # 

# remove records with 0 values of TSF. 
bckgr_all_tsf <- filter(bckgr_all_tsf, tsf >= 1) # n =  535923




#### Now assign tsf to tgb

tgb_tsf <- tgb

tgb_tsf$tsf<-NA

# valNF <- 117    # oldest fires for the study area.   Already set

# spatial points
coordinates(tgb_tsf) <- c(6,7)
crs(tgb_tsf) <- ("+init=epsg:3577")


## Extract fires
fires.ext_tgb <- over(tgb_tsf, fire_poly, returnList = TRUE) 
sites_fires_tgb <- sapply(fires.ext_tgb, nrow) 


for(i in 1:length(sites_fires_tgb)){ 
  if(sites_fires_tgb[[i]]==0){ 
    tgb_tsf$tsf[i]<-valNF
  }else{
    tsf<-tgb_tsf$year[i] - fires.ext_tgb[[i]]$FireYr 
    tsf[tsf<0]<-NA 
    tsf<-c(tsf,valNF)
    tgb_tsf$tsf[i]<-min(tsf,na.rm=T)
  }
}


# separate info rndm and bias background
tgb_tsf <- as.data.frame(tgb_tsf) 

# remove records with 0 values of TSF. 
tgb_tsf <- filter(tgb_tsf, tsf >= 1) # n =  41650




##################################################################################################
#################### Extract TSF to fitting datasets #############################################

# first create id column
PA$id_set <- "PA"
PA_PIA$id_set <- "PA_PIA"
pvol_pre500$id_set <- "PB"

# remove unnecesary columns
PA_PIA <- dplyr::select(PA_PIA, -c("count", "date", "accuracy", "detection.method"))
PA <- dplyr::select(PA, -c("count", "date", "accuracy", "detection.method"))

fit.set_all <- rbind(PA, PA_PIA, pvol_pre500) # bind to process altogether

fit.set_all$tsf<-NA

# valNF <- 117    # oldest fires for the study area.   Already set

# spatial points
coordinates(fit.set_all) <- c(6,7)
crs(fit.set_all) <- ("+init=epsg:3577")

## Extract fires
fires.ext_fit.set_all <- over(fit.set_all, fire_poly, returnList = TRUE) 
sites_fires_fit.set_all <- sapply(fires.ext_fit.set_all, nrow) 

for(i in 1:length(sites_fires_fit.set_all)){ 
  if(sites_fires_fit.set_all[[i]]==0){ 
    fit.set_all$tsf[i]<-valNF
  }else{
    tsf<-fit.set_all$year[i] - fires.ext_fit.set_all[[i]]$FireYr 
    tsf[tsf<0]<-NA 
    tsf<-c(tsf,valNF)
    fit.set_all$tsf[i]<-min(tsf,na.rm=T)
  }
}


# separate info rndm and bias background
fit.set_all <- as.data.frame(fit.set_all) 

# remove records with 0 values of TSF. 
fit.set_all <- filter(fit.set_all, tsf >= 1) 




##################################################################################################
#################### Extract variables to all datasets ###########################################


##### Prepare variables

# All variables have been resampled (500 m), projected (epsg:3577) and masked to modelling extent
vars_stack_all <- raster::stack(list.files(paste0(getwd(),"./spatial_data/variables_model"), pattern = '.asc', full.names=TRUE))
crs(vars_stack_all) <- ("+init=epsg:3577")

names(vars_stack_all)

# Change names to variables in stack
vars.names <- c("temp_mean", 
                "temp_warm", 
                "temp_cold",
                "pp_annual", 
                "pp_driest", 
                "temp_season", 
                "pp_season", 
                "mean_drought", 
                "fPAR_mean", 
                "fPAR_var", 
                "GPP_mean", 
                "GPP_var", 
                "mean_EDI", 
                "tsf", 
                "c_height",
                "mean_xveg")

names(vars_stack_all) <- vars.names
names(vars_stack_all)

# pairs(vars_stack_all, maxpixels = 10000) # check correlation

### variables dropped

# Mean temperature: Highly correlated with temperature and precipitation variables that better describe extreme conditions (mean_temp_warmest(0.97), pp_driest(0.8), pp_seasonality(0.8), mean_drought_run(0.8)). It was preffered to use these extreme variables rather than mean conditions along the year.

# Mean Temperature of Coldest Quarter: highly correlated with temp_warm (0.88) and pp_season(0.85) High temperatures (temp_warm) are more likely a limiting factor, so preference is given to it

# Mean drough run: highly correalted with other precipitation and temperature variables that were created at a finner resolution (ANUCLIM variables created by J. Elith at ~271 m).

# GPP mean: Highly correlated with fPAR mean. It was decided to use fPAR over GPP as fPAR, which summarizes vegetation greenness is a more direct surrogate of nitrogen content in foliage.


# drop variables
names(vars_stack_all)
vars_stack <- dropLayer(vars_stack_all, c(1,3,8,11))
names(vars_stack)
rm(vars_stack_all)

# pairs(vars_stack, maxpixels = 10000) # verify



##### Extract variables to all datasets 

#### Do extraction to all fitting datasets and background datasets together.
## datasets so far:
# bckgr_all_tsf: rndm, bmodel and bfile background points
# tgb_tsf: tgb with tsf
# fit.set_all: all fitting datasets (evaluation dataset has not yet been set aside)

### Make all columns the same for all sets

str(bckgr_all_tsf)
str(tgb_tsf)
str(fit.set_all)

bckgr_all_tsf$species <- "bckgr_point"
bckgr_all_tsf$PA <- "bckgr_point"
colnames(bckgr_all_tsf)[4] <- "id_set" # same name for id of dataset as the other datasets

colnames(tgb_tsf)[c(2,6,7)] <- c("PA", "x","y")
tgb_tsf$species <- "tgb_point"
tgb_tsf$PA <- "bckgr_point"
tgb_tsf$id_set <- "tgb"
tgb_tsf <- dplyr::select(tgb_tsf, -c("source", "dataset.name")) # these won't be necessary anymore

colnames(fit.set_all)[c(2,6,7)] <- c("PA", "x","y")
fit.set_all <- dplyr::select(fit.set_all, -c("source", "dataset.name", "ID")) # these won't be necessary anymore

# same order of columns
bckgr_all_tsf <- dplyr::select(bckgr_all_tsf, c("species", "PA", "year", "x", "y", "id_set", "tsf")) 
tgb_tsf <- dplyr::select(tgb_tsf, c("species", "PA", "year", "x", "y", "id_set", "tsf")) 

# check
str(bckgr_all_tsf)
str(tgb_tsf)
str(fit.set_all)

# merge all together to do extraction
extract_all <- rbind(bckgr_all_tsf, tgb_tsf, fit.set_all)

# extract excluding tsf raster(to keep values extracted previously)
coordinates(extract_all) <- c(4,5)
crs(extract_all) <- ("+init=epsg:3577")

extract_all <- as.data.frame(raster::extract(subset(vars_stack, c(1:9,11,12)), extract_all, sp=TRUE))

extract_all <- na.omit(extract_all)

# save(extract_all, file = "./outputs/extract_all.RData") # can save 

#### separate into sets

unique(extract_all$id_set)
# rndm   
# bmodel 
# bfile  
# tgb    
# PA     
# PA_PIA 
# PB     

rndm_extract <- extract_all[extract_all$id_set %in% "rndm", ]
bmodel_extract <- extract_all[extract_all$id_set %in% "bmodel", ]
bfile_extract <- extract_all[extract_all$id_set %in% "bfile", ]
tgb_extract <- extract_all[extract_all$id_set %in% "tgb", ]
PA_extract <- extract_all[extract_all$id_set %in% "PA", ]
PA_PIA_extract <- extract_all[extract_all$id_set %in% "PA_PIA", ]
PB_extract <- extract_all[extract_all$id_set %in% "PB", ]

# Save these as the datasets to fit final models with ALL data available. 
# Background datasets are subsampled to 100,000 points here. These 100,000 points are sampled irrespective of external evaluation blocks as these objects will be used to fit final models only.

rndm_bckgr_fmodel <- dplyr::sample_n(rndm_extract, size = 100000, replace = FALSE)
bmodel_bckgr_fmodel <- dplyr::sample_n(bmodel_extract, size = 100000, replace = FALSE)
bfile_bckgr_fmodel <- dplyr::sample_n(bfile_extract, size = 100000, replace = FALSE)

tgb_bckgr_fmodel <- as.data.frame(tgb_extract)
PA_fmodel <- as.data.frame(PA_extract)
PA_PIA_fmodel <- as.data.frame(PA_PIA_extract)
PB_fmodel <- as.data.frame(PB_extract)

# # save objects
# save(rndm_bckgr_fmodel, file="./outputs/rndm_bckgr_fmodel.RData")
# save(bmodel_bckgr_fmodel, file="./outputs/bmodel_bckgr_fmodel.RData")
# save(bfile_bckgr_fmodel, file="./outputs/bfile_bckgr_fmodel.RData")
# save(tgb_bckgr_fmodel, file="./outputs/tgb_bckgr_fmodel.RData")
# save(PA_fmodel, file="./outputs/PA_fmodel.RData")
# save(PA_PIA_fmodel, file="./outputs/PA_PIA_fmodel.RData")
# save(PB_fmodel, file="./outputs/PB_fmodel.RData")




##################################################################################################
######### Subset spatially into training and and evaluation datasets for external evaluation #####

# prepare data to use blockCV
PA_PIA_blocks <- PA_PIA_extract
PA_PIA_blocks <- dplyr::select(PA_PIA_blocks, c("PA", "x", "y"))

PA_PIA_blocks_pre <- PA_PIA_blocks %>% dplyr::filter(PA %in% 1)
PA_PIA_blocks_abs <- PA_PIA_blocks %>% dplyr::filter(PA %in% 0)
id_eval <- c(rep(1, nrow(PA_PIA_blocks_pre)), rep(0, nrow(PA_PIA_blocks_abs)))

locations_eval <- as.data.frame(rbind(PA_PIA_blocks_pre[,2:3], PA_PIA_blocks_abs[,2:3]))
loc4block_eval <- data.frame(id_eval, locations_eval)

loc.df_eval  <- SpatialPointsDataFrame(loc4block_eval[,c("x", "y")], loc4block_eval)
crs(loc.df_eval) <- ("+init=epsg:3577")

# Note: Allocation of presence and absence records to folds will change from run to run, since 'maskBySpecies' was set to FALSE. This allows for a small degree of stochasticity in blocks/folds allocation. Hence, we've provided the set of blocks we used for the models presented in the paper. Also, generation of blocks can take a long time.

# blocks with systematic asignation
# sb_eval_syst <- spatialBlock(speciesData = loc.df_eval,
#                              species = "id_eval",
#                              rasterLayer = mask500m,
#                              theRange = 25000, # 25 km size
#                              k = 3,
#                              selection = 'systematic',
#                              maskBySpecies = F,
#                              showBlocks = T,
#                              biomod2Format = FALSE)
#
#
# 
# 
# # mapview(sb_eval_syst$blocks)
# 
# #### separate the 3 sets of external evaluation blocks
# blocks_syst <- sb_eval_syst$blocks
# 
## extract the three sets of blocks
# 
# # 1
# block1_syst <- blocks_syst[blocks_syst$folds == 1,]
# 
# # 2
# block2_syst <- blocks_syst[blocks_syst$folds == 2,]
# 
# # 3
# block3_syst <- blocks_syst[blocks_syst$folds == 3,]
# 
# rm(blocks_syst) # too heavy, remove

# load blocks
block1_syst <- shapefile("./spatial_data/block1_syst.shp")
block2_syst <- shapefile("./spatial_data/block2_syst.shp")
block3_syst <- shapefile("./spatial_data/block3_syst.shp")

#### get location of blocks that have evaluation data

coordinates(PA_PIA_extract) <- c(4,5)
crs(PA_PIA_extract) <- ("+init=epsg:3577")

### set 1
blocks_eval.data_syst1 <- sp::over(x =PA_PIA_extract, y = block1_syst)
# remove NAs 
blocks_eval.data_syst1 <- na.omit(blocks_eval.data_syst1) # These are the centroids of the eval blocks. Use these to get final evaluation blocks

coordinates(blocks_eval.data_syst1) <- c(3,4)
crs(blocks_eval.data_syst1) <- ("+init=epsg:3577")

#### get blocks that overlap PA_PIA data
block1_syst_eval <- block1_syst[blocks_eval.data_syst1,]
plot(block1_syst_eval)



### set 2
blocks_eval.data_syst2 <- sp::over(x =PA_PIA_extract, y = block2_syst)
# remove NAs 
blocks_eval.data_syst2 <- na.omit(blocks_eval.data_syst2) 

coordinates(blocks_eval.data_syst2) <- c(3,4)
crs(blocks_eval.data_syst2) <- ("+init=epsg:3577")

block2_syst_eval <- block2_syst[blocks_eval.data_syst2,]
plot(block2_syst_eval)



### set 3
blocks_eval.data_syst3 <- sp::over(x =PA_PIA_extract, y = block3_syst)
# remove NAs 
blocks_eval.data_syst3 <- na.omit(blocks_eval.data_syst3) 

coordinates(blocks_eval.data_syst3) <- c(3,4)
crs(blocks_eval.data_syst3) <- ("+init=epsg:3577")

block3_syst_eval <- block3_syst[blocks_eval.data_syst3,]
plot(block3_syst_eval)

# keep only eval blocks with PA_PIA data. remove sets of blocks and centroids
rm(block1_syst, block2_syst, block3_syst, blocks_eval.data_syst1, blocks_eval.data_syst2, blocks_eval.data_syst3)

#### Get PA_PIA data from evaluation blocks = Evaluation datasets

## Set 1
PA_PIA_eval_1 <- PA_PIA_extract[block1_syst_eval,]

## Set 2
PA_PIA_eval_2 <- PA_PIA_extract[block2_syst_eval,]

## Set 3
PA_PIA_eval_3 <- PA_PIA_extract[block3_syst_eval,]





##################################################################################################
######### Remove external evaluation records from fitting sets ###################################

#### Remove fitting records from eval blocks. fitting datasets: 
# PA_extract
# PA_PIA_extract
# PB_extract
# rndm_extract 
# bmodel_extract 
# bfile_extract
# tgb_extract

#### Do for each of the 3 sets of external evaluation data and then review spatially (manually) that there is no overlap.

# At the end there will be 3 sets of each model fitting dataset. E.g. PA_fit1, PA_fit2, PA_fit3


######## Set 1 ######## 

##### PA_extract
coordinates(PA_extract) <- c(4,5)
crs(PA_extract) <- "+init=epsg:3577"

dat_in <-which(!is.na(over(PA_extract, block1_syst_eval)[,1])) 
PA_fit_1 <- PA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially.
PA_veri <-c(0,1)
test<-which(over(PA_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri) #  Verify for the corresponding evaluation data set. E.g. PA_PIA_eval_1, PA_PIA_eval_2, PA_PIA_eval_3
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(PA_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PA_PIA_extract
dat_in <-which(!is.na(over(PA_PIA_extract, block1_syst_eval)[,1])) 
PA_PIA_fit_1 <- PA_PIA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PA_PIA_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(PA_PIA_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PB_extract
coordinates(PB_extract) <- c(4,5)
crs(PB_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(PB_extract, block1_syst_eval)[,1])) 
PB_fit_1 <- PB_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PB_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(PB_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### rndm_extract
coordinates(rndm_extract) <- c(4,5)
crs(rndm_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(rndm_extract, block1_syst_eval)[,1])) 
rndm_fit_1 <- rndm_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(rndm_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(rndm_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bmodel_extract
coordinates(bmodel_extract) <- c(4,5)
crs(bmodel_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bmodel_extract, block1_syst_eval)[,1])) 
bmodel_fit_1 <- bmodel_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bmodel_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(bmodel_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bfile_extract
coordinates(bfile_extract) <- c(4,5)
crs(bfile_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bfile_extract, block1_syst_eval)[,1])) 
bfile_fit_1 <- bfile_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bfile_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(bfile_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### tgb_extract
coordinates(tgb_extract) <- c(4,5)
crs(tgb_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(tgb_extract, block1_syst_eval)[,1])) 
tgb_fit_1 <- tgb_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(tgb_fit_1,PA_PIA_eval_1)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_1, add=TRUE)
plot(tgb_fit_1[test,], add=TRUE, col="red") # no dups
rm(PA_veri)

#### review plots
dev.off() # remove plots



######## Set 2 ######## 

##### PA_extract
dat_in <-which(!is.na(over(PA_extract, block2_syst_eval)[,1])) 
PA_fit_2 <- PA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially.
PA_veri <-c(0,1)
test<-which(over(PA_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri) #  Verify for the corresponding evaluation data set. E.g. PA_PIA_eval_1, PA_PIA_eval_2, PA_PIA_eval_3
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(PA_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PA_PIA_extract
dat_in <-which(!is.na(over(PA_PIA_extract, block2_syst_eval)[,1])) 
PA_PIA_fit_2 <- PA_PIA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PA_PIA_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(PA_PIA_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PB_extract
coordinates(PB_extract) <- c(4,5)
crs(PB_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(PB_extract, block2_syst_eval)[,1])) 
PB_fit_2 <- PB_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PB_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(PB_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### rndm_extract
coordinates(rndm_extract) <- c(4,5)
crs(rndm_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(rndm_extract, block2_syst_eval)[,1])) 
rndm_fit_2 <- rndm_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(rndm_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(rndm_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bmodel_extract
coordinates(bmodel_extract) <- c(4,5)
crs(bmodel_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bmodel_extract, block2_syst_eval)[,1])) 
bmodel_fit_2 <- bmodel_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bmodel_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(bmodel_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bfile_extract
coordinates(bfile_extract) <- c(4,5)
crs(bfile_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bfile_extract, block2_syst_eval)[,1])) 
bfile_fit_2 <- bfile_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bfile_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(bfile_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### tgb_extract
coordinates(tgb_extract) <- c(4,5)
crs(tgb_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(tgb_extract, block2_syst_eval)[,1])) 
tgb_fit_2 <- tgb_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(tgb_fit_2,PA_PIA_eval_2)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_2, add=TRUE)
plot(tgb_fit_2[test,], add=TRUE, col="red") # no dups
rm(PA_veri)

#### review plots
dev.off() # remove plots



######## Set 3 ########

##### PA_extract
dat_in <-which(!is.na(over(PA_extract, block3_syst_eval)[,1])) 
PA_fit_3 <- PA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PA_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(PA_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PA_PIA_extract
dat_in <-which(!is.na(over(PA_PIA_extract, block3_syst_eval)[,1])) 
PA_PIA_fit_3 <- PA_PIA_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PA_PIA_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(PA_PIA_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### PB_extract
coordinates(PB_extract) <- c(4,5)
crs(PB_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(PB_extract, block3_syst_eval)[,1])) 
PB_fit_3 <- PB_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(PB_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(PB_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### rndm_extract
coordinates(rndm_extract) <- c(4,5)
crs(rndm_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(rndm_extract, block3_syst_eval)[,1])) 
rndm_fit_3 <- rndm_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(rndm_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(rndm_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bmodel_extract
coordinates(bmodel_extract) <- c(4,5)
crs(bmodel_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bmodel_extract, block3_syst_eval)[,1])) 
bmodel_fit_3 <- bmodel_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bmodel_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(bmodel_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### bfile_extract
coordinates(bfile_extract) <- c(4,5)
crs(bfile_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(bfile_extract, block3_syst_eval)[,1])) 
bfile_fit_3 <- bfile_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(bfile_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(bfile_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


##### tgb_extract
coordinates(tgb_extract) <- c(4,5)
crs(tgb_extract) <- ("+init=epsg:3577")
dat_in <-which(!is.na(over(tgb_extract, block3_syst_eval)[,1])) 
tgb_fit_3 <- tgb_extract[-c(dat_in),] 
rm(dat_in)

# Verify spatially
PA_veri <-c(0,1)
test<-which(over(tgb_fit_3,PA_PIA_eval_3)[,2] %in% PA_veri)
plot(bckgr_poly)
plot(PA_PIA_eval_3, add=TRUE)
plot(tgb_fit_3[test,], add=TRUE, col="red") # no dups
rm(PA_veri)


#### review plots
dev.off() # remove plots



### get the final sample of 100000 for background points

### Set 1
rndm_fit_1 <- as.data.frame(rndm_fit_1)
bmodel_fit_1 <- as.data.frame(bmodel_fit_1)
bfile_fit_1 <- as.data.frame(bfile_fit_1)

rndm_bckgr_1 <- dplyr::sample_n(rndm_fit_1, size = 100000, replace = FALSE)
bmodel_bckgr_1 <- dplyr::sample_n(bmodel_fit_1, size = 100000, replace = FALSE)
bfile_bckgr_1 <- dplyr::sample_n(bfile_fit_1, size = 100000, replace = FALSE)

### Set 2
rndm_fit_2 <- as.data.frame(rndm_fit_2)
bmodel_fit_2 <- as.data.frame(bmodel_fit_2)
bfile_fit_2 <- as.data.frame(bfile_fit_2)

rndm_bckgr_2 <- dplyr::sample_n(rndm_fit_2, size = 100000, replace = FALSE)
bmodel_bckgr_2 <- dplyr::sample_n(bmodel_fit_2, size = 100000, replace = FALSE)
bfile_bckgr_2 <- dplyr::sample_n(bfile_fit_2, size = 100000, replace = FALSE)

### Set 3
rndm_fit_3 <- as.data.frame(rndm_fit_3)
bmodel_fit_3 <- as.data.frame(bmodel_fit_3)
bfile_fit_3 <- as.data.frame(bfile_fit_3)

rndm_bckgr_3 <- dplyr::sample_n(rndm_fit_3, size = 100000, replace = FALSE)
bmodel_bckgr_3 <- dplyr::sample_n(bmodel_fit_3, size = 100000, replace = FALSE)
bfile_bckgr_3 <- dplyr::sample_n(bfile_fit_3, size = 100000, replace = FALSE)


# ### save R objects
# 
# ### Set 1 
# save(PA_fit_1, file = "./outputs/PA_fit_1.RData")
# save(PA_PIA_fit_1, file = "./outputs/PA_PIA_fit_1.RData")
# save(PB_fit_1, file = "./outputs/PB_fit_1.RData")
# save(rndm_bckgr_1, file = "./outputs/rndm_bckgr_1.RData")
# save(bmodel_bckgr_1, file = "./outputs/bmodel_bckgr_1.RData")
# save(bfile_bckgr_1, file = "./outputs/bfile_bckgr_1.RData")
# save(tgb_bckgr_1, file = "./outputs/tgb_bckgr_1.RData")
# 
# ### Set 2 
# save(PA_fit_2, file = "./outputs/PA_fit_2.RData")
# save(PA_PIA_fit_2, file = "./outputs/PA_PIA_fit_2.RData")
# save(PB_fit_2, file = "./outputs/PB_fit_2.RData")
# save(rndm_bckgr_2, file = "./outputs/rndm_bckgr_2.RData")
# save(bmodel_bckgr_2, file = "./outputs/bmodel_bckgr_2.RData")
# save(bfile_bckgr_2, file = "./outputs/bfile_bckgr_2.RData")
# save(tgb_bckgr_2, file = "./outputs/tgb_bckgr_2.RData")
# 
# ### Set 3 
# save(PA_fit_3, file = "./outputs/PA_fit_3.RData")
# save(PA_PIA_fit_3, file = "./outputs/PA_PIA_fit_3.RData")
# save(PB_fit_3, file = "./outputs/PB_fit_3.RData")
# save(rndm_bckgr_3, file = "./outputs/rndm_bckgr_3.RData")
# save(bmodel_bckgr_3, file = "./outputs/bmodel_bckgr_3.RData")
# save(bfile_bckgr_3, file = "./outputs/bfile_bckgr_3.RData")
# save(tgb_bckgr_3, file = "./outputs/tgb_bckgr_3.RData")
# 
# # evaluation datasets
# save(PA_PIA_eval_1, file = "./outputs/PA_PIA_eval_1.RData")
# save(PA_PIA_eval_2, file = "./outputs/PA_PIA_eval_2.RData")
# save(PA_PIA_eval_3, file = "./outputs/PA_PIA_eval_3.RData")



##################################################################################################
######### Evaluation datasets VIC and NSW ########################################################

# obtain evaluation records from southern half of study area

# load VIC and NSW polygon
vic_poly <- shapefile("./spatial_data/vic3577.shp")
crs(vic_poly) <- ("+init=epsg:3577")

nsw_poly <- shapefile("./spatial_data/vic3577.shp")
crs(nsw_poly) <- ("+init=epsg:3577")

vic_nsw <- bind(vic_poly, nsw_poly) # combine int single polygon

# Set_1
coordinates(PA_PIA_eval_1) <- c(4,5)
crs(PA_PIA_eval_1) <- ("+init=epsg:3577")
PA_PIA_eval_1_south <- as.data.frame(PA_PIA_eval_1[vic_nsw,])

# Set_2
coordinates(PA_PIA_eval_2) <- c(4,5)
crs(PA_PIA_eval_2) <- ("+init=epsg:3577")
PA_PIA_eval_2_south <- as.data.frame(PA_PIA_eval_2[vic_nsw,])

# Set_3
coordinates(PA_PIA_eval_3) <- c(4,5)
crs(PA_PIA_eval_3) <- ("+init=epsg:3577")
PA_PIA_eval_3_south <- as.data.frame(PA_PIA_eval_3[vic_nsw,])

# # save objects
# save(PA_PIA_eval_1_south, file ="./outputs/PA_PIA_eval_1_south.RData")
# save(PA_PIA_eval_2_south, file ="./outputs/PA_PIA_eval_2_south.RData")
# save(PA_PIA_eval_3_south, file ="./outputs/PA_PIA_eval_3_south.RData")
