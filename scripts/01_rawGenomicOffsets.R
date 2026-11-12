# Script to model genomic variation with gradient forest to calculate raw genomic offsets as described in 
# Lachmuth et al. (2024).
# Code originally written by Susanne Lachmuth, UMD, 2020-2023
# Code last edited by Nora Heaphy, UVM, 11 Nov 2026.

# load packages
# note that extendedForest and gradientForest need to be manually repaired in the package directory after
# installation in order to work
library(adehabitatLT)
library(extendedForest)
library(fields)
library(gradientForest)
library(reshape2)
library(parallel)
library(doParallel)
library(data.table)
library(sf)
library(tidyverse)
library(tmap)
library(ggplot2)
library(viridis)
library(ggmap)


# GRADIENT FOREST MODELING -----------------------------------------------------
# all this only needs to be done once
# for different future time periods and climate projections, start from spatio-temporal offsets section

# read in data for 11 climate variables for 64 sampled red spruce populations
climGF <- read.csv(file = "data/clim_sprucePops.csv")

# read in allele frequencies of 240 candidate SNPs for 64 sampled red spruce populations
# these frequencies are ancestral vs. derived not major vs. minor allele
# pops are in same order as in climate data file
snpScores <- read.csv(file = "data/snpScores_sprucePops.csv")

# model of allele frequency turnover along climatic gradients
maxLevel <- log2(0.368*nrow(climGF)/2) # accounts for correlations, see ?gradientForest

# fit gf models for candidate SNPs
gfMod <- gradientForest(cbind(climGF, snpScores), predictor.vars = colnames(climGF),
                           response.vars = colnames(snpScores), ntree = 500,
                           maxLevel = maxLevel, trace = T, corr.threshold = 0.5)

#### transform current climate data into gene-environment space ####

# read data for 11 climate variables for all geographic grid cells (2.5 arc-min) with red spruce presence
# records (8821 lat/lons)
# this file already includes the orchard site with cellindex = 999999
clim_sprucePres_xy <- fread("data/currClim_sprucePres_1960_1990.csv", stringsAsFactors = T, header = T)

# subset to climate data only
clim_sprucePres <- clim_sprucePres_xy[,-c(1:3)]

# transform climate data using gradient forest model predictions
transClim_xy <- predict(gfMod, clim_sprucePres)

# put back in xy coordinates
transClim_xy$x <- clim_sprucePres_xy$x
transClim_xy$y <- clim_sprucePres_xy$y
transClim_xy$cellindex <- clim_sprucePres_xy$cellindex

# subset to just climate data
transClim <- transClim_xy %>% dplyr::select(CMD:PET)


# CONTEMPORARY SPATIAL OFFSETS --------------------------------------------
# all this only needs to be done once here to generate ECDF for standardizing spatial offsets
# for different future time periods and climate projections, start from spatio-temporal offsets section

# settings
nCores <- 50
nBreak <- 100

# run in parallel
breakIt <- split(1:nrow(transClim), cut(1:nrow(transClim), nCores, labels = FALSE))

cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)

# calculate pairwise offset matrix between all presences under current climates
spatOffset <- foreach(i = 1:length(breakIt), .packages = c("fields")) %dopar%{
  mat <- rdist(transClim[breakIt[[i]],], transClim)
  return(mat)
}

stopCluster(cl)

# reformat
spatOffset_list <- do.call(rbind, spatOffset)
spatOffset_df <- as.data.frame(spatOffset_list)
colnames(spatOffset_df) <- rownames(transClim)
spatOffset_mat <- as.matrix(spatOffset_df)
spatOffset_lon <- reshape2::melt(spatOffset_mat)[reshape2::melt(lower.tri(spatOffset_mat))$value,]

# get empirical cumulative density function
ecdfSpatOffset <- ecdf(spatOffset_lon$value)

# save for standardized offsets script
save(spatOffset_df, file = "results/spatOffset_df.RData")
save(ecdfSpatOffset, file = "results/ecdfSpatOffset.RData")


# CHANGE EXTENT --------------------------------------------
# crop extent to just study area: PA, NY, NJ, ME, NH, VT, MA, CT, RI

# read in US states shapefiles (NAD83)
states <- st_read("data/US_states.shp")

# subset to only the states we want
study_states <- states %>% filter(NAME %in% c("Maine", "New Hampshire", "Vermont", "Massachusetts",
                                              "Connecticut", "New York", "Pennsylvania", "New Jersey",
                                              "Rhode Island"))

# convert transClim into spatial points
transClim_pts <- st_as_sf(transClim_xy, coords = c("x", "y"), crs = 4326)

# re-project points to NAD83 to match US states CRS
transClim_pts_proj <- st_transform(transClim_pts, "+proj=longlat +datum=NAD83")

# use state boundaries to clip points
transClim_reduced <- st_intersection(transClim_pts_proj, st_union(study_states))

# retrieve xy coords
transClim_xy <- transClim_reduced %>% st_drop_geometry() %>%
    mutate(x = st_coordinates(transClim_reduced)[,1], y = st_coordinates(transClim_reduced)[,2])

# subset to just climate data
transClim <- transClim_xy %>% dplyr::select(CMD:PET)

# clip untransformed points to study area for later
clim_sprucePres_pts <- st_as_sf(clim_sprucePres_xy, coords = c("x", "y"), crs = 4326)
clim_sprucePres_pts_proj <- st_transform(clim_sprucePres_pts, "+proj=longlat +datum=NAD83")
clim_sprucePres_reduced <- st_intersection(clim_sprucePres_pts_proj, st_union(study_states))
clim_sprucePres_xy <- clim_sprucePres_reduced %>% st_drop_geometry() %>%
    mutate(x = st_coordinates(clim_sprucePres_reduced)[,1], y = st_coordinates(clim_sprucePres_reduced)[,2])

# save for standardized offsets script
write.csv(clim_sprucePres_xy, "data/clim_sprucePres_reduced.csv", row.names = F)

# CONTEMPORARY SPATIAL OFFSETS --------------------------------------------
# re-calculate contemporary spatial offset just for points in northeastern USA

# settings
nCores <- 50
nBreak <- 100

# run in parallel
breakIt <- split(1:nrow(transClim), cut(1:nrow(transClim), nCores, labels = FALSE))

cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)

# calculate pairwise offset matrix between all presences under current climates
spatOffset <- foreach(i = 1:length(breakIt), .packages = c("fields")) %dopar%{
  mat <- rdist(transClim[breakIt[[i]],], transClim)
  return(mat)
}

stopCluster(cl)

# reformat
spatOffset_list <- do.call(rbind, spatOffset)
spatOffset_df <- as.data.frame(spatOffset_list)
colnames(spatOffset_df) <- rownames(transClim)
spatOffset_mat <- as.matrix(spatOffset_df)
spatOffset_lon <- reshape2::melt(spatOffset_mat)[reshape2::melt(lower.tri(spatOffset_mat))$value,]

# save for Bristol offsets
save(spatOffset_df, file = "results/spatOffset_studyArea_df.RData")


# SPATIO-TEMPORAL OFFSETS -------------------------------------------------
# Spatio-temporal offsets describe the G - E disruption for a (donor) population being
# transferred from one grid cell under current climate to any other grid cell (recipient) under future climate.
# repeat for each future time period and climate projection SSP
# SSP2-45 (moderate emissions) and 2071-2100 shown here as an example

# read in future climate (2071-2100, SSP245)
futClim_studyArea_xy <- fread("data/futClim_2071_SSP245.csv", stringsAsFactors = T, header = T)

# transform into gene-environment space
futClim_studyArea <- futClim_studyArea_xy[,CMD:PET]
futClim_trans_xy <- data.frame(futClim_studyArea_xy[,c("x","y","cellindex")], predict(gfMod,futClim_studyArea))
futClim_trans_xy <- futClim_trans_xy[order(futClim_trans_xy$cellindex),]

#### change extent ####

# convert to spatial points
future_trans_pts <- st_as_sf(futClim_trans_xy, coords = c("x", "y"), crs = 4326)

# re-project to NAD83 to match U.S. states
future_trans_pts_proj <- st_transform(future_trans_pts, "+proj=longlat +datum=NAD83")

# use state boundaries to clip future_points
future_trans_reduced <- st_intersection(future_trans_pts_proj, st_union(study_states))

# retrieve xy coords
future_trans_reduced_xy <- future_trans_reduced %>% st_drop_geometry() %>%
    mutate(x = st_coordinates(future_trans_reduced)[,1], y = st_coordinates(future_trans_reduced)[,2])

## clip untransformed future points for later
future_pts <- st_as_sf(futClim_studyArea_xy, coords = c("x", "y"), crs = 4326)
future_pts_proj <- st_transform(future_pts, "+proj=longlat +datum=NAD83")
future_reduced <- st_intersection(future_pts_proj, st_union(study_states))
future_reduced_xy <- future_reduced %>% st_drop_geometry() %>%
    mutate(x = st_coordinates(future_reduced)[,1], y = st_coordinates(future_reduced)[,2])

# save for standardized offsets script
write.csv(future_reduced_xy, "data/future_reduced.csv", row.names = F)

# rename for downstream analyses
currClim_trans_xy <- transClim_xy
futClim_trans_xy <- future_trans_reduced_xy

# remove xy coordinates and cell index
currClim_trans <- currClim_trans_xy[,-(12:14)]
futClim_trans <- futClim_trans_xy[,c(2:12)]

### calculate raw offsets between all donor (current) and recipient (future) cells ------------------

# settings
nCores <- 50
nBreak <- 100

breakItAll <- split(1:nrow(futClim_trans), cut(1:nrow(futClim_trans), nBreak, labels = FALSE))

cl <- makeCluster(nCores)
registerDoParallel(cl)

# calculate pairwise offset matrix between all current donors and future recipients
allOffset <- foreach(i = 1:length(breakItAll), .packages = c("fields")) %dopar%{
  mat_allOffset <- rdist(futClim_trans[breakItAll[[i]],], currClim_trans)
  return(mat_allOffset)
}

stopCluster(cl)

# reformat data
allOffset_list <- do.call(rbind, allOffset)
rm(allOffset) # clear more memory, the next one needs a lot
allOffset_df <- as.data.frame(allOffset_list)
rm(allOffset_list) # clear more memory

# save for standardized offsets script
save(allOffset_df, file = "results/allOffset_df_2071_SSP245.RData")