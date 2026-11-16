# Compares recipient importance based on genomic offset-based candidate donor selection vs. 
# alternative methods of seed sourcing (random, local to orchard site, within same seed zone, 
# one seed zone north, from within seed zones containing top genomic offset candidates)
# Code last edited by Nora Heaphy, UVM, 16 Nov 2026.

library(tidyterra)
library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(parallel)
library(doParallel)
library(adehabitatLT)
library(viridis)
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(ggsci)
library(ggridges)
library(ggpp)
library(ggpubr)
library(scales)


# REC IMP FOR TOP G.O. DONORS FOR COMPARISON -------------------------------------
# this was also calculated in 03_recImp_donorsScions.R

# read in background shapefile
states <- st_read("data/US_states.shp")
study_states <- states %>% filter(NAME %in% c("Maine", "New Hampshire", "Vermont", "Massachusetts",
                                              "Connecticut", "New York", "Pennsylvania", "New Jersey",
                                              "Rhode Island"))

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
recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4326)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_pts, aes(color = percRecImp), size = 0.5, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    limits = c(0, 100), name = "Recipient importance (%)") +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Donors selected by genomic offset model",
  subtitle = "N = 233") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())

# NULL MODEL: DONORS LOCAL TO ORCHARD SITE -------------------------------------

clim_sprucePres_pts <- st_as_sf(clim_sprucePres_xy, coords = c("x", "y"), crs = 4326) 

# project to a meters coordinate system (NAD83 / UTM zone 19N)
clim_sprucePres_proj <- st_transform(clim_sprucePres_pts, crs = 26919)
bristol <- clim_sprucePres_proj[4710, ]
buffer <- st_buffer(bristol, dist = 50000)  # 50 km buffer
bristol_local_points <- clim_sprucePres_proj[st_within(clim_sprucePres_proj, buffer, sparse = FALSE), ]
bristol_local_pts <- st_transform(bristol_local_points, crs = 4269)
bristol_local <- clim_sprucePres_xy %>% 
  filter(clim_sprucePres_xy$cellindex %in% bristol_local_points$cellindex)

candidate_indices <- which(clim_sprucePres_xy$x %in% bristol_local$x & clim_sprucePres_xy$y %in% bristol_local$y)
candidate_indices <- paste0("V", candidate_indices, sep = "")
candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_indices)]
recImp_xy <- read.csv("results/recImp_sigmaTH15_2071_SSP245.csv") %>%
    select(x, y, cellindex)
recImp_xy$recImp <- rowSums(candidate_offsets, na.rm = T) 
recImp_xy$percRecImp <- recImp_xy$recImp/ncol(candidate_offsets)*100
recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4269)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_pts, aes(color = percRecImp), size = 0.6, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    name = "Recipient importance (%)", limits = c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Donors local to orchard site",
  subtitle = "N = 168") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())

# NULL MODEL: RANDOM SELECTION -------------------------------------

x <- 1:ncol(allOffset_sigmaTH_df)
candidate_indices <- sample(x, size = 200, replace = FALSE)
candidate_indices <- sort(candidate_indices)
candidate_indices <- paste0("V", candidate_indices, sep = "")
candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_indices)]
recImp_xy <- read.csv("results/recImp_sigmaTH15_2071_SSP245.csv") %>%
    select(x, y, cellindex)
recImp_xy$recImp <- rowSums(candidate_offsets, na.rm = T) 
recImp_xy$percRecImp <- recImp_xy$recImp/ncol(candidate_offsets)*100
recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4269)

# loop
recImp_aggr <- read.csv("results/recImp_sigmaTH15_2071_SSP245.csv") %>%
    select(x, y, cellindex)
x <- 1:ncol(allOffset_sigmaTH_df)
for (i in 1:200) {
  candidate_indices <- sample(x, size = 200, replace = FALSE)
  candidate_indices <- sort(candidate_indices)
  candidate_indices <- paste0("V", candidate_indices, sep = "")
  candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_indices)]
  new_col1 <- paste0("recImp", i)
  recImp_aggr[[new_col1]] <- rowSums(candidate_offsets, na.rm = T) 
  new_col2 <- paste0("percRecImp", i)
  recImp_aggr[[new_col2]] <- recImp_aggr[[new_col1]]/ncol(candidate_offsets)*100
  recImp_aggr[[new_col1]] <- NULL
}
recImp_aggr$avg_recImp <- rowMeans(recImp_aggr[, 4:203])
recImp_pts <- st_as_sf(recImp_aggr, coords = c("x", "y"), crs = 4269)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_pts, aes(color = avg_recImp), size = 0.6, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    name = "Recipient importance (%)", limits = c(0,100)) +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Randomly selected donors",
  subtitle = "N = 200") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())

# NULL MODEL: SAME SEED ZONE -------------------------------------

seed_zones <- st_read("data/ESZ_v2_2.shp")
clim_sprucePres_pts <- st_as_sf(clim_sprucePres_xy, coords = c("x", "y"), crs = 4326)
clim_sprucePres_proj <- st_transform(clim_sprucePres_pts, crs = st_crs(seed_zones))
clim_sprucePres_zones <- st_join(clim_sprucePres_proj, seed_zones)
currClim_zones_sp <- st_transform(clim_sprucePres_zones, crs = 4326)
currClim_zones <- st_drop_geometry(currClim_zones_sp) %>%
      mutate(x = st_coordinates(currClim_zones_sp)[,1],
             y = st_coordinates(currClim_zones_sp)[,2])

futClim_xy <- read.csv("data/future_reduced.csv")
futClim_pts <- st_as_sf(futClim_xy, coords = c("x", "y"), crs = 4326)
futClim_proj <- st_transform(futClim_pts, crs = st_crs(seed_zones))
futClim_zones <- st_join(futClim_proj, seed_zones)
futClim_zones_sp <- st_transform(futClim_zones, crs = 4326)
futClim_zones <- st_drop_geometry(futClim_zones_sp) %>%
      mutate(x = st_coordinates(futClim_zones_sp)[,1],
             y = st_coordinates(futClim_zones_sp)[,2])

currClim_zones$SeedZone <- as.factor(currClim_zones$SeedZone)
futClim_zones$SeedZone <- as.factor(futClim_zones$SeedZone)
# seed zone 37 has no candidate donors, so pull from nearest adjacent zone 16
futClim_zones[which(futClim_zones$SeedZone == "37"),]$SeedZone = "16"
futClim_zones$SeedZone <- droplevels(futClim_zones$SeedZone)
# manual fixes for zones that failed to join
futClim_zones[which(futClim_zones$cellindex %in% c(2830860,2885535)),]$SeedZone <- "1"
futClim_zones[which(futClim_zones$cellindex %in% c(2927998,2937115)),]$SeedZone <- "25"
futClim_zones[which(futClim_zones$cellindex %in% c(3165050,3162016,3095207,3095209,
                                 3092172,3086097,3073946,3070908,
                                 3067870,3067869,3064831,3061792,
                                 3052716,3055761,3046647,3025416)),]$SeedZone <- "16"
futClim_zones[which(futClim_zones$cellindex %in% c(3055811,3058847,3086135)),]$SeedZone <- "17"
futClim_zones[which(is.na(futClim_zones$SeedZone)),]$SeedZone <- "4"
currClim_zones[which(is.na(currClim_zones$SeedZone)),]$SeedZone <- "4"
futClim_zones$SeedZone <- droplevels(futClim_zones$SeedZone)
currClim_zones$SeedZone <- droplevels(currClim_zones$SeedZone)
table(currClim_zones$SeedZone)

### loop
recImp_xy <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(recImp_xy) <- c("x", "y", "cellindex", "recImp", "percRecImp",
                         "seedZone")
donors <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(donors) <- c("x", "y", "cellindex")
seed_zone_nums <- c("1","14","15","16","17","18","2","25","3","4","5",
                    "92","93","94","95","97","98")
for (i in 1:length(seed_zone_nums)){
  print(paste0("Starting seed zone ", seed_zone_nums[i]))
  candidate_indices <- which(currClim_zones$SeedZone == seed_zone_nums[i])
  if (length(candidate_indices) > 200) {
    candidate_indices <- sample(candidate_indices, 200, replace = FALSE)
    print(paste0("> 200 donors in seed zone ", seed_zone_nums[i]))
  } else {
    print(paste0("< 200 donors in seed zone ", seed_zone_nums[i]))
  }
  candidate_names <- paste0("V", candidate_indices, sep = "")
  candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_names)]
  candidate_offsets <- as.data.frame(candidate_offsets)
  candidate_coords <- clim_sprucePres_xy[candidate_indices,] %>%
        dplyr::select(x, y, cellindex)
  donors <- rbind(donors, candidate_coords)
  recip_indices <- which(futClim_zones$SeedZone == seed_zone_nums[i])
  recip_offsets <- candidate_offsets[which(rownames(candidate_offsets) %in% recip_indices),]
  recip_offsets <- as.data.frame(recip_offsets)
  futClim_sub <- futClim_zones %>% filter(SeedZone == seed_zone_nums[i])
  recImp_temp <- data.frame(futClim_sub[,c("x","y","cellindex")])
  recImp_temp$recImp <- rowSums(recip_offsets, na.rm = T) 
  recImp_temp$percRecImp <- (recImp_temp$recImp/200)*100
  recImp_temp$seedZone <- seed_zone_nums[i]
  recImp_xy <- rbind(recImp_xy, recImp_temp)
  print(paste0("Finished seed zone ", seed_zone_nums[i]))
}

recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4326) 
recImp_proj <- st_transform(recImp_pts, crs = 4269)
recImp_xy$seedZone <- as.factor(recImp_xy$seedZone)
summ <- recImp_xy %>% dplyr::select(percRecImp, seedZone) %>% 
  group_by(seedZone) %>%
  summarize(avgRecImp = mean(percRecImp), medianRecImp = median(percRecImp)) %>%
  arrange(medianRecImp)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_proj, aes(color = percRecImp), size = 0.6, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    name = "Recipient importance (%)", limits = c(0,100)) +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Donors from within the same seed zone",
  subtitle = "N = up to 200 per seed zone") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())

# NULL MODEL: ONE SEED ZONE NORTH -------------------------------------

# seed zone donors are drawn from the zone bordering the recipient zone that has the farthest southern extent
recImp_xy <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(recImp_xy) <- c("x", "y", "cellindex", "recImp", "percRecImp", "seedZone")
seed_zone_recips <- c("1","2","3","4","5","14","15","16","17","18","25","92","93","94","95","97","98")
seed_zone_donors <- c("3","94","15","14","98","97","16","16","98","16","14","14","94","14","14","98","16")
donors <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(donors) <- c("x", "y", "cellindex")

for (i in 1:length(seed_zone_recips)){
  print(paste0("Starting seed zone ", seed_zone_recips[i]))
  candidate_indices <- which(currClim_zones$SeedZone == seed_zone_donors[i])
  if (length(candidate_indices) > 200) {
    candidate_indices <- sample(candidate_indices, 200, replace = FALSE)
    print(paste0("> 200 donors in seed zone ", seed_zone_donors[i]))
  } else {
    print(paste0("< 200 donors in seed zone ", seed_zone_donors[i]))
  }
  candidate_names <- paste0("V", candidate_indices, sep = "")
  candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_names)]
  candidate_offsets <- as.data.frame(candidate_offsets)
  candidate_coords <- clim_sprucePres_xy[candidate_indices,] %>%
        dplyr::select(x, y, cellindex)
  donors <- rbind(donors, candidate_coords)
  recip_indices <- which(futClim_zones$SeedZone == seed_zone_recips[i])
  recip_offsets <- candidate_offsets[which(rownames(candidate_offsets) %in% recip_indices),]
  recip_offsets <- as.data.frame(recip_offsets)
  futClim_sub <- futClim_zones %>% filter(SeedZone == seed_zone_recips[i])
  recImp_temp <- data.frame(futClim_sub[,c("x","y","cellindex")])
  recImp_temp$recImp <- rowSums(recip_offsets, na.rm = T) 
  recImp_temp$percRecImp <- (recImp_temp$recImp/200)*100
  recImp_temp$seedZone <- seed_zone_recips[i]
  recImp_xy <- rbind(recImp_xy, recImp_temp)
  print(paste0("Finished seed zone ", seed_zone_recips[i]))
}

recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4326) 
recImp_proj <- st_transform(recImp_pts, crs = 4269)
recImp_xy$seedZone <- as.factor(recImp_xy$seedZone)
summ <- recImp_xy %>% dplyr::select(percRecImp, seedZone) %>% 
  group_by(seedZone) %>%
  summarize(avgRecImp = mean(percRecImp), medianRecImp = median(percRecImp)) %>%
  arrange(medianRecImp)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_proj, aes(color = percRecImp), size = 0.6, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    name = "Recipient importance (%)", limits = c(0,100)) +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Donors from one seed zone south",
  subtitle = "N = up to 200 per seed zone") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())

# NULL MODEL: SEED ZONES CONTAINING TOP GO DONORS -------------------------------------

seed_zones <- st_read("data/ESZ_v2_2.shp")
cands <- read.csv("results/top_cands.csv")
cands_pts <- st_as_sf(cands, coords = c("x", "y"), crs = 4326)
cands_proj <- st_transform(cands_pts, crs = st_crs(seed_zones))
cands_zones <- st_join(cands_proj, seed_zones)
cands_zones <- st_transform(cands_zones, crs = 4326)
summ <- cands_zones %>% group_by(SeedZone) %>%
  summarize(N = n())
# these are the seed zones with >5 G.O. top donors
seed_zones_sub <- seed_zones %>% filter(SeedZone %in% c(4,14,15,95))
sub_zones_proj <- st_transform(seed_zones_sub, crs = st_crs(study_states))
cands_proj <- st_transform(cands_pts, crs = st_crs(study_states))

candidate_indices <- which(currClim_zones$SeedZone %in% c(4,14,15,95))
candidate_indices <- sample(candidate_indices, 200, replace = FALSE)
candidate_coords <- clim_sprucePres_xy[candidate_indices,] %>%
        dplyr::select(x, y, cellindex)
candidate_indices <- paste0("V", candidate_indices, sep = "")
candidate_offsets <- allOffset_sigmaTH_df[, which(names(allOffset_sigmaTH_df) %in% candidate_indices)]
recImp_xy <- data.frame(futClim_xy[,c("x","y","cellindex")])
recImp_xy$recImp <- rowSums(candidate_offsets, na.rm = T)
recImp_xy$percRecImp <- (recImp_xy$recImp/200)*100
recImp_pts <- st_as_sf(recImp_xy, coords = c("x", "y"), crs = 4326) 
recImp_proj <- st_transform(recImp_pts, crs = 4269)

ggplot() +
  geom_sf(data = study_states, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = recImp_proj, aes(color = percRecImp), size = 0.6, alpha = 0.9) +
  scale_color_viridis(option = "viridis", direction = -1,
    name = "Recipient importance (%)", limits = c(0,100)) +
  scale_x_continuous(labels = function(x) x) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(x = "Longitude", y = "Latitude", title = "Donors from seed zones containing top G.O. candidates",
  subtitle = "N = 200") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(family = "sans", color = "black"),
    axis.title = element_text(family = "sans", color = "black"),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank(),
    legend.title = element_text(margin = ggplot2::margin(b = 15)),
    legend.key = element_blank())