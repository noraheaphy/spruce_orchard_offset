# Calculates recipient importance for top candidate donors as well as actual scion collected for orchard.
# Code last edited by Nora Heaphy, UVM, 16 Nov 2026.

library(tidyterra)
library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(parallel)
library(doParallel)
library(adehabitatLT)
library(gradientForest)
library(extendedForest)
library(viridis)
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(ggsci)
library(ggridges)
library(ggpp)
library(ggpubr)
library(scales)

# REC IMP FOR TOP G.O. CANDIDATES -----------------------------------------------------

# read in top candidates selected by genomic offset model
donors <- read.csv("results/top_cands.csv")

# read in spatiotemporal offset matrix
load("results/allOffset_sigmaTH15_2071_SSP245.RData", verbose = TRUE)

# read in all donors
clim_sprucePres_xy <- read.csv("data/clim_sprucePres_reduced.csv")

# use candidate xy coords to find index numbers of those coords in clim_sprucePres_xy
candidate_indices <- which(clim_sprucePres_xy$x %in% donors$x & clim_sprucePres_xy$y %in% donors$y)

# add a capital V in front of every index in the candidate_indices vector to match allOffset_df columns
candidate_indices <- paste0("V", candidate_indices, sep = "")

# use indices to subset allOffset_df to only include columns with those indices
candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_indices)]

# make a dataframe for recipient importance
recImp_xy <- read.csv("results/recImp_sigmaTH15_2071_SSP245.csv") %>%
    select(x, y, cellindex)

# calculate cumulative offsets for recipient points with the 233 best candidate donors
recImp_xy$recImp <- rowSums(candidate_offsets, na.rm = T) 
recImp_xy <- recImp_xy %>% mutate(percRecImp = (recImp/ncol(candidate_offsets))*100)

write.csv(recImp_xy, "results/recImp_top_cands.csv")

# SCION OFFSETS AND REC IMP -----------------------------------------------------

scions_values <- read.csv("data/scion_collections_env_extract_SSP245_2071_2100.csv")

# gradient forest model (explained further in 01_rawGenomicOffsets.R)
climGF <- read.csv(file = "data/clim_sprucePops.csv")
snpScores <- read.csv(file = "data/snpScores_sprucePops.csv")
maxLevel <- log2(0.368*nrow(climGF)/2)
gfMod <- gradientForest(cbind(climGF, snpScores), predictor.vars = colnames(climGF),
                           response.vars = colnames(snpScores), ntree = 500,
                           maxLevel = maxLevel, trace = T, corr.threshold = 0.5)

# transform climate data
scions_clim <- scions_values[,c(9:19)]
transClim <- predict(gfMod, scions_clim)
transClim_xy <- transClim
transClim_xy$x <- scions_values$Long
transClim_xy$y <- scions_values$Lat
transClim_xy$Tree.ID <- scions_values$Tree.ID

# transform future climate data for recipient points
futClim_xy <- read.csv("data/future_reduced.csv")
futClim <- futClim_xy[,c(2:12)]
futClim_trans_xy <- data.frame(futClim_xy[,c("x","y","cellindex")], predict(gfMod, futClim))
futClim_trans_xy <- futClim_trans_xy[order(futClim_trans_xy$cellindex),]
currClim_trans <- transClim_xy[,-(12:14)]
futClim_trans <- futClim_trans_xy[,c(4:14)]

# calculate raw offsets between all scion (current) and recipient (future) cells 
nCores <- 50
nBreak <- 100
breakItAll <- split(1:nrow(futClim_trans), cut(1:nrow(futClim_trans), nBreak, labels = FALSE))
cl <- makeCluster(nCores)
registerDoParallel(cl)
allOffset <- foreach(i = 1:length(breakItAll), .packages = c("fields")) %dopar%{
  mat_allOffset <- rdist(futClim_trans[breakItAll[[i]],], currClim_trans)
  return(mat_allOffset)
}
stopCluster(cl)
allOffset_list <- do.call(rbind, allOffset)
rm(allOffset) 
allOffset_df <- as.data.frame(allOffset_list)
rm(allOffset_list)

# standardize offsets
load("results/ecdfSpatOffset.RData", verbose = TRUE)
breakIt_prob <- split(1:nrow(allOffset_df), cut(1:nrow(allOffset_df), nBreak, 
                                                labels = FALSE))
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
rm(allOffset_df, breakIt_prob)
allOffset_prob_list <- do.call(rbind, allOffset_prob_loop)
rm(allOffset_prob_loop)
allOffset_prob_df <- as.data.frame(allOffset_prob_list)
rm(allOffset_prob_list)
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
rm(allOffset_prob_df, breakIt_qchi)
allOffset_sigma_list <- do.call(rbind, allOffset_sigma_loop)
rm(allOffset_sigma_loop)
allOffset_sigma_df <- as.data.frame(allOffset_sigma_list)
rm(allOffset_sigma_list)
allOffset_sigma_df[1:5,1:5]
save(allOffset_sigma_df, file = "results/scions_allOffset_sigma_df_2071_SSP245.RData")

# donor and recipient importance
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
rm(breakIt_sigma)
allOffset_sigmaTH_list <- do.call(rbind, allOffset_sigmaTH_loop)
rm(allOffset_sigmaTH_loop)
allOffset_sigmaTH_df <-as.data.frame(allOffset_sigmaTH_list)
rm(allOffset_sigmaTH_list)
row.names(allOffset_sigmaTH_df) <- row.names(allOffset_sigma_df)
save(allOffset_sigmaTH_df, file = "results/scions_allOffset_sigmaTH15_df_2071_SSP245.RData")

donImp_xy <- data.frame(scions_values[,c("Long","Lat","Tree.ID")])
donImp_xy$donImp <- colSums(allOffset_sigmaTH_df, na.rm = T) 
donImp_xy$percDonImp <- donImp_xy$donImp/nrow(allOffset_sigmaTH_df)*100

recImp_xy <- data.frame(futClim_xy[, c("x", "y", "cellindex")]) 
recImp_xy$recImp <- rowSums(allOffset_sigmaTH_df,na.rm = T)
recImp_xy$percRecImp <- recImp_xy$recImp/ncol(allOffset_sigmaTH_df)*100

write.csv(donImp_xy, "results/scions_donor_importance_TH15_2071_SSP245.csv", row.names = F)
write.csv(recImp_xy, "results/scions_recipient_importance_TH15_2071_SSP245.csv", row.names = F)