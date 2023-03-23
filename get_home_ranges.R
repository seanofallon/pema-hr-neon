# This script uses NSF-NEON small mammal box trapping capture data to calculate home ranges
# for Peromyscus maniculatus (PEMA) individuals, and assign traits and conditions to those individuals
# for the analysis performed in 'Uncovering multiple influences on space use by deer mice using
# NEON data' by S. O'Fallon, K. Mabry, and N. Pinter-Wollman. For details on NEON data collection
# and data products, refer to their website: https://data.neonscience.org/

# INPUTS:
# pemajn_caps.csv - NEON-recorded data for each PEMA capture w/ tagID in their 'small mammal box
#                   trapping' protocol. Used to calculate home ranges and assign traits/conditions.
# mnka-bp.csv - Minimum number (of PEMA) known alive, representing population density, at each 
#               small mammal sampling plot at each sampling event.

# OUTPUTS:
# ao_ranges - dataframe containing home range area calculations and trait/condition assignments
#             for all eligible PEMA individuals analyzed in this study. Used as input for statistical
#             analysis and data visualization in 'analyze_home_ranges.R'


##### Start up: clear env, load packages, read in data
# clear env
rm(list=ls())

# load packages
library(dplyr)
library(neonUtilities)
library(neonOS)
library(raster)
library(adehabitatHR)
library(stringr)

# read in pemajn_caps and mnka-bp (if these are in your working directory, run as is. Otherwise,
# change file name to the full path name to where you've saved pemajn_caps.csv and mnka-bp.csv)
pemajn_caps<-read.csv("pemajn_caps.csv")
mnka.bp<-read.csv("mnka-bp.csv")


### Reduce pemajn_caps df to PEMA captured at least 5 times
# find tagIDs w/ at least 5 caps
recap_ids <- pemajn_caps %>% count(tagID) %>% filter(n >= 5)
# filter capture df to only tags captured >4x
pema_recaps <- filter(pemajn_caps, tagID %in% recap_ids$tagID)

### Set outlier mice to be removed - found in later steps to have data entry errors
# remove '148' weight that was entered incorrectly in data entry - replace w/ NA
pema_recaps["weight"][pema_recaps["weight"] == 148] <- NA
# remove 'NEON.MAM.D01.R1761', who had errors in pregnancyStatus
pema_recaps["weight"][pema_recaps["tagID"] == 'NEON.MAM.D01.R1761'] <- NA

### Add col w/ lifestage based on weight to restrict dataset to adult mice
for (i in 1:nrow(pema_recaps)) {
  if (is.na(pema_recaps[i,'weight']) == T) {
    pema_recaps[i,'lifeStageByWeight']<-'unknown'
  }else if (pema_recaps[i,'weight']<14) {
    pema_recaps[i,'lifeStageByWeight']<-'juvenile'
  }else if (pema_recaps[i,'weight']>=14 & pema_recaps[i,'weight']<=15.5) {
    pema_recaps[i,'lifeStageByWeight']<-'subadult'
  }else if (pema_recaps[i,'weight']>15.5) {
    pema_recaps[i,'lifeStageByWeight']<-'adult'
  } else
    pema_recaps[i,'lifeStageByWeight']<-'unknown'
}

# loop goes through all indivs in recap_ids, replacing all lifeStageByWeight
# after first 'adult' w/ 'adult' (mice stay adults after growing to adulthood)
for (i in 1:nrow(recap_ids)) {
  imouse<-filter(pema_recaps,tagID==recap_ids$tagID[i]) # take all caps for ith mouse
  if ('adult' %in% imouse$lifeStageByWeight) { # if we ever assign 'adult' to imouse
    # assign date of first 'adult' assignment to first.adult
    first.adult<-min(imouse[imouse$lifeStageByWeight=='adult','collectDate.x'])
    # put 'adult' for lifeStageByWeight in imouse for all caps after 1st 'adult' assignment
    imouse[imouse$collectDate.x >= first.adult,'lifeStageByWeight'] <- 'adult'
    # apply to pema_recaps
    pema_recaps[pema_recaps$tagID==recap_ids$tagID[i],'lifeStageByWeight']<-
      imouse$lifeStageByWeight
  }
}

### Filter pema_recaps to mice captured at least 5x as an adult - our dataset
# remove non-adult captures
ao_recaps<-filter(pema_recaps,lifeStageByWeight=='adult')
# find tagIDs w/ at least 5 adult caps
ao_ids <- ao_recaps %>% count(tagID) %>% filter(n >= 5)
# filter capture df to only tags captured >4x
ao_recaps <- filter(ao_recaps, tagID %in% ao_ids$tagID)


##### Calculate home ranges
# set up df to hold home ranges
ao_ranges <- data.frame(tagID=character(),
                        MCP=numeric(),
                        UD=numeric())
# Function that converts letters to numbers, to deal w/ NEON coords:
LETTER2num <- function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}
# make objects for HR loop
data=ao_recaps # data to be fed to loop
unq_id=unique(data$tagID) # all IDs - to dictate # of iterations

# get HR of all individuals:
for(i in 1:length(unq_id)){
  ix=unq_id[i]
  ao_ranges[i,'tagID']<-ix
  traps=data$trapCoordinate[data$tagID==ix] # take the column with the trap coordinates of a single individual from the NEON data
  splt= data.table::tstrsplit(traps, split="")
  x_coor=10*unlist(lapply(splt[[1]],FUN=LETTER2num)) # turn A,B,C etc to numbers and multiply by 10 to get units in meters
  if(max(x_coor)>200){
    x_coor[x_coor==240]=NA  
  }
  if(length(splt)==2){ # if there are no traps at row 10
    y_coor=10*as.numeric(unlist(splt[[2]])) # multiplied by 10 to get units in meters
  }
  if(length(splt)==3){# if there are traps at row 10
    y_coor=10*as.numeric(unlist(splt[[2]])) # multiplied by 10 to get units in meters
    y_coor[!is.na(splt[[3]])]=100
  }
  
  coor=cbind(x_coor,y_coor) # bind xy into a single variable
  coor=coor[!is.na(coor[,1]),] #remove NAs - if X was entered instead of a letter, this will remove those entries
  coor=coor[!is.na(coor[,2]),] #remove NAs - if X was entered instead of a number, this will remove those entries
  coor = as.data.frame(coor)
  for_samp = dim(coor)[1]
  # get observed HR
  if(dim(unique(coor))[1]>2 & nrow(coor)>4){
    coor = sp::SpatialPoints(coor) # coordinates need to be in the format that adehabitatHR likes:
    # get MCP area:
    ao_ranges[i,'MCP']=as.numeric(mcp.area(coor, percent = 95,unin = "m", unout = "m2", plotit = FALSE)) # get the area - in m^2 - or the MCP. you can add as.numeric() around this and save in a variable that has all the areas for all the animals
    # get area based on kernel
    ud_est= kernelUD(coor) # make the UD of the animal
    ao_ranges[i,'UD']=as.numeric(kernel.area(ud_est, percent = 50, unin = "m", unout = "m2")) # get the area of the 50% UD of the animal. you can add as.numeric() around this and save in a variable that has all the areas for all the animals
  }else{
    ao_ranges[i,'MCP']=NA
    ao_ranges[i,'UD'] = NA
  }
}

# remove NAs from ao_ranges
ao_ranges <- na.omit(ao_ranges)


### Attach vars to ao_ranges for modeling
for (i in 1:nrow(ao_ranges)) { # for each range calculated:
  # take adult caps of ith mouse
  imouse <- filter(ao_recaps,tagID==ao_ranges$tagID[i])
  ### SEX
  # if mouse was ever pregnant, assign it F in ao_ranges
  if ('pregnant' %in% imouse$pregnancyStatus) {
    ao_ranges[i,'sex']<-'F'}
  # if mouse assigned M > F, assign it M in ao_ranges
  else if ('M' > 'F' %in% imouse$sex) {
    ao_ranges[i,'sex']<-'M'}
  # if mouse assigned F > M, assign it F in ao_ranges
  else if ('M' < 'F' %in% imouse$sex) {
    ao_ranges[i,'sex']<-'F'}
  # that should cover all cases - but if it doesn't, assign U for unknown
  else {
    ao_ranges[i,'sex']<-'U'}
  ### SITE
  ao_ranges[i,'site'] <- imouse[1,"siteID"]
  ### PLOT ID (needed for MNKA assignment...)
  ao_ranges[i,'plotID'] <- imouse[1,"plotID"]
  ### VEG TYPE
  # if imouse caps occured in forests/woodlands, assign to vegType 'forest'
  if (imouse$nlcdClass[1]=='deciduousForest' | imouse$nlcdClass[1]=='mixedForest'
      | imouse$nlcdClass[1]=='evergreenForest' | imouse$nlcdClass[1]=='woodyWetlands') {
    ao_ranges[i,'vegType']<-'forest'}
  # if imouse caps occured in grassy habitat, assign to vegType 'grassland'
  else if (imouse$nlcdClass[1]=='cultivatedCrops' | imouse$nlcdClass[1]=='grasslandHerbaceous') {
    ao_ranges[i,'vegType']<-'grassland'}
  # if imouse caps occured in shrubby habitat, assign to vegType 'shrubland'
  else if (imouse$nlcdClass[1]=='shrubScrub') {
    ao_ranges[i,'vegType']<-'shrubland'}
  # that should cover everything, but in case of errors, assign U for unknown
  else  {
    ao_ranges[i,'vegType']<-'U'
  }
  ### HINDFOOT LENGTH
  ao_ranges[i,'meanHindfootLength']<-mean(imouse[imouse$pregnancyStatus!='pregnant' & imouse$pregnancyStatus!='unknown',]$hindfootLength,na.rm=T)  
  ### WEIGHT
  ao_ranges[i,'meanWeight']<-mean(imouse[imouse$pregnancyStatus!='pregnant' & imouse$pregnancyStatus!='unknown',]$weight,na.rm=T)
  ### RANGE OF DATES COLLECTED
  ao_ranges[i,'firstCollected']<-min(imouse$collectDate.x)
  ao_ranges[i,'lastCollected']<-max(imouse$collectDate.x)
  ### YEAR
  ao_ranges[i,'year']<-str_trunc(ao_ranges[i,'lastCollected'],4,side='right',ellipsis = '')
  ### LATITUDE
  ao_ranges[i,'latitude'] <- imouse[1,'decimalLatitude']
  ### MIN NUMBER KNOWN ALIVE (POP DENSITY)
  mouse_alive<-pemajn_caps[pemajn_caps$collectDate.x >= ao_ranges$firstCollected[i] &
                             pemajn_caps$collectDate.x <= ao_ranges$lastCollected[i] &
                             pemajn_caps$site == ao_ranges$site[i] &
                             pemajn_caps$plotID == ao_ranges$plotID[i],]

  all_mnka4mouse<-mnka.bp[mnka.bp$eventID %in% unique(mouse_alive$eventID) &
                            mnka.bp$plotID == unique(mouse_alive$plotID),'n']
  ao_ranges[i,'meanMNKA']<-mean(all_mnka4mouse,na.rm=T)
  ao_ranges[i,'maxMNKA']<-max(all_mnka4mouse,na.rm=T)
}


### Calculate bodyCondition based on residuals of relationship btwn meanHindfootLength and meanWeight
# first, eliminate NAs from ao_ranges as they cannot be included in lm
ao_ranges<-na.omit(ao_ranges)
# calculate relationship btwn vars
ao.mw.by.mhl<-lm(meanWeight ~ meanHindfootLength, ao_ranges)
# match residuals to ranges
for (i in 1:nrow(ao_ranges)) {
  ao_ranges[i,'bodyCondition']<-ao.mw.by.mhl$residuals[i]
}
