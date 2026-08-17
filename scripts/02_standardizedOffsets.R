# Calculates standardized offsets from spatio-temporal offsets, as well as donor and recipient importance 
# based on  methodology described in Lachmuth et al. (2024).
# Code originally written by Susanne Lachmuth, UMD, 2020-2023
# Code last edited by Nora Heaphy, UVM, 17 Aug 2026.

# load packages
library(adehabitatLT)
library(extendedForest) 
library(fields)
library(gradientForest)
library(reshape2)
library(parallel)
library(doParallel)
library(data.table)
library(tidyverse)
library(sf)
library(terra)
library(tmap)

# set parameters
nCores <- 50
nBreak <- 100

# read in spatiotemporal offset matrix
load("results/allOffset_df_2071_SSP245.RData", verbose = TRUE)

# load in ecdf function
load("results/ecdfSpatOffset.RData", verbose = TRUE)

# SCALE SPATIO-TEMPORAL OFFSETS -------------------------------------
# Some of these steps are quite computationally intensive and use up a lot of memory. Consider
# deleting intermediate objects whenever possible.

## re-express raw offsets as probability based on spatial offset ecdf  -------------

# prep for loop
breakIt_prob <- split(1:nrow(allOffset_df), cut(1:nrow(allOffset_df), nBreak, labels = FALSE))

# run in parallel
cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_prob_loop <- foreach(i = 1:length(breakIt_prob), .packages = c("stats")) %dopar%{
  allOffset_prob_sub <- allOffset_df[breakIt_prob[[i]],]
  for(j in 1:ncol(allOffset_prob_sub)) {
    allOffset_prob_sub[ , j] <- ecdfSpatOffset(allOffset_prob_sub[ , j]) 
  } 
  return(allOffset_prob_sub)
}

stopCluster(cl)

# reformat data
allOffset_prob_list <- do.call(rbind, allOffset_prob_loop)
allOffset_prob_df <- as.data.frame(allOffset_prob_list)
rm(allOffset_prob_list, allOffset_prob_loop)

# sanity check: should be numbers between 0 and 1
allOffset_prob_df[1:5,1:5]

## re-express raw offsets as quantiles (sigmas) of chi distribution with 1 degree of freedom (half-normal distribution)  -------------

breakIt_qchi <- split(1:nrow(allOffset_prob_df), cut(1:nrow(allOffset_prob_df), nBreak, labels = FALSE))

cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_sigma_loop <- foreach(i = 1:length(breakIt_qchi), .packages=c("adehabitatLT")) %dopar%{
  allOffset_sigma_sub <- allOffset_prob_df[breakIt_qchi[[i]],]
  for(j in 1:ncol(allOffset_sigma_sub)) {       
    allOffset_sigma_sub[ , j] <- adehabitatLT::qchi(allOffset_sigma_sub[ , j],1) 
  } 
  return(allOffset_sigma_sub)
}

stopCluster(cl)

# reformat data
allOffset_sigma_list <- do.call(rbind, allOffset_sigma_loop)
allOffset_sigma_df <- as.data.frame(allOffset_sigma_list)
rm(allOffset_prob_df, allOffset_sigma_loop, allOffset_sigma_list)

# sanity check: should be numbers between 0 and approx 3
allOffset_sigma_df[1:5,1:5]

save(allOffset_sigma_df, file = "results/allOffset_standardized_2071_SSP245.RData")

### DONOR & RECIPIENT IMPORTANCE --------------------------------------------

## apply sigma threshold ---------------------------------------------------
# set "not to exceed" offset threshold (sigma = 1, 1.5, or 2) and calculate transferability matrix
sigmaTH <- 1.5

breakIt_sigma <- split(1:nrow(allOffset_sigma_df), cut(1:nrow(allOffset_sigma_df), nBreak, labels = FALSE))

cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_sigmaTH_loop <- foreach(i = 1:length(breakIt_sigma)) %dopar%{
  allOffset_sigmaTH_sub<-allOffset_sigma_df[breakIt_sigma[[i]],]
  for  (j in 1:ncol(allOffset_sigmaTH_sub)){
    allOffset_sigmaTH_sub[,j]<-ifelse(allOffset_sigmaTH_sub[,j]<=sigmaTH,1,0)
  }
  return(allOffset_sigmaTH_sub)
}

stopCluster(cl)

# reformat data
allOffset_sigmaTH_list <- do.call(rbind, allOffset_sigmaTH_loop)
allOffset_sigmaTH_df <-as.data.frame(allOffset_sigmaTH_list)
row.names(allOffset_sigmaTH_df) <- row.names(allOffset_sigma_df)
rm(allOffset_sigmaTH_list, allOffset_sigmaTH_loop, allOffset_sigma_df)

# sanity check: should be all 1s and 0s
allOffset_sigmaTH_df[1:5,1:5]

save(allOffset_sigmaTH_df, file = "results/allOffset_sigmaTH15_2071_SSP245.RData")

## calculate donor importance ----------------------------------------------

# read in current spruce presences (candidates)
clim_sprucePres_xy <- fread("data/clim_sprucePres_reduced.csv", stringsAsFactors = T, header = T)

# make donor importance dataframe
donImp_xy <- data.frame(clim_sprucePres_xy[,c("x","y","cellindex")])

# calculate donor importance
donImp_xy$donImp <- colSums(allOffset_sigmaTH_df, na.rm = T) 
donImp_xy$percDonImp <- donImp_xy$donImp/nrow(allOffset_sigmaTH_df)*100

write.csv(donImp_xy, "results/donImp_sigmaTH15_2071_SSP245.csv", row.names = F)

## calculate recipient importance ------------------------------------------

# read in future grid cells (recipients)
futClim_studyArea_xy <- fread("data/future_reduced.csv", stringsAsFactors = T, header = T)

# make recipient importance dataframe
recImp_xy <- data.frame(futClim_studyArea_xy[, c("x", "y", "cellindex")]) 

# calculate recipient importance
recImp_xy$recImp<-rowSums(allOffset_sigmaTH_df,na.rm = T)
recImp_xy$percRecImp<-recImp_xy$recImp/ncol(allOffset_sigmaTH_df)*100

write.csv(recImp_xy, "results/recImp_sigmaTH15_2071_SSP245.csv", row.names = F)