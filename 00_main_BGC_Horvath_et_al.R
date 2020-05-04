#### Improving the representation of high-latitude vegetation in Dynamic Global Vegetation Models ####
# Dynamical Vegetation model version CLM-DV4.5
# updated: 30.4.2020
# author: Peter Horvath
# based on data created also by Hui Tang
# part of scripts reused from: http://www.neonscience.org/field-data-polygons-centroids

#### 00 ####
#### LOAD LIBRARIES #######################################

library(raster)
library(maptools)
library(sf)
library(rgdal)
library(rgeos) 
library(sp)
library(ggplot2)
#library(data.table) 
#library(dplyr)
library(vegan)

#ploting
library(ggplot2)
library(scales)
library(gridExtra)
#### Define paths ####
# all AR18x18 dissolved data

proj_folder <- "E:/Project_2"
in_NOR_data <- "E:/Project_2/HUI_CLM_layers/Norway_borders_N50"
in_veg_data <- "E:/Project_2/HUI_CLM_layers/AR18x18_plots"
in_pts_data <- "E:/Project_2/Single point plots CLM4.5"
in_RS_data <- "E:/Project_2/LAND COVER DATA/Satveg_deling_nd_partnere_09_12_2009/tiff/"
in_binary_raster <- "E:/Project_1_RUNS/new_MODEL_RUN/07_proportion_raster/Binary_raster/"
in_probab_raster <- "E:/Project_1_RUNS/new_MODEL_RUN/04_predict_raster/"
in_raster <- "E:/Project_1_FINAL_layers/MASKED_TIFF/"
in_dm_data <- "E:/Project_2/OUTPUT/"
in_csv_files <- "C:/Users/peterhor/Documents/GitHub/Project_2/"
out_folder <- "E:/Project_2_RUNS/OUTPUT"
out_ggplots <- "E:/Project_2_RUNS/OUTPUT/results/ggplots"
dir.create(file.path(out_folder),showWarnings=F)

#### CRS:  WGS84 UTM zone 33N ####
project_crs <- crs("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
sr <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#### Read in data ####
PFT <- read.csv(paste0(in_csv_files,"02_Translation_schemes/PFT.csv"), sep = ";", dec = ".")
# NOR_border <- readOGR(in_NOR_data, "Landmf_N50_dissolved")
# summary(NOR_border)
# NOR_border_df <- fortify(NOR_border) #for use with ggplot2
# plot(NOR_border, col="lightgreen")

# Remote sensing rasters from SatVeg (https://kartkatalog.miljodirektoratet.no/MapService/Details/satveg)
RS_north <- raster(paste0(in_RS_data,"SatVeg_NN_30m.tif"))
#crs(RS_north) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
RS_south <- raster(paste0(in_RS_data,"SatVeg_SN_30m.tif"))
RS_Norway <-merge(RS_south, RS_north, filename=paste0(in_RS_data,"/RS_Norway"), format="GTiff", overwrite=TRUE)
crs(RS_Norway)
# RS_Norway <- raster(paste0(in_RS_data,"RS_Norway.tif"))
RS_Norway33 <- projectRaster(RS_Norway, crs=sr, filename=paste0(in_RS_data,"/RS_Norway_utm33"), format="GTiff", overwrite=TRUE)
#RS_Norway33 <- raster(paste0(in_RS_data,"RS_Norway_utm33.tif"))
#crs(RS_south) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# plot(RS_north)
# read in DM *preliminary results from 3 DM attempts
DM_maxval <- raster(paste0(out_folder,"/DM_maxVal/DM_maxval.tif"))
DM_preval <- raster(paste0(out_folder,"/DM_max_preval/DM_preval_fill.tif"))
DM_AUC <- raster(paste0(out_folder,"/DM_max_AUC/DM_AUC_fill.tif"))
DM_maxval_PFT <- raster(paste0(out_folder,"/DM_maxVal/DM_maxval_PFT.tif"))
DM_preval_PFT <- raster(paste0(out_folder,"/DM_max_preval/DM_preval_PFT.tif"))
DM_AUC_PFT <- raster(paste0(out_folder,"/DM_max_AUC/DM_AUC_PFT.tif"))
crs(DM_maxval)
# read in DM *preliminary results from BinaryTreshold attempt


#### Create 1km boxes 20x####
veg_plots <- readOGR(in_veg_data, "Hui_dissolved_plots")
veg_plots@proj4string
#need for transformation

# veg_plots <- spTransform(veg_plots_trans, CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# veg_plots@proj4string
#subset
# plot ID numbers to be subset
test_20_plots <- c(	405,	513,	622,	801,	922,	1131,
                    1304,	1322,	1623,	2015,	2108,	2238,	2332,	2425,
                    2948,	2962,	4268,	5369,	6380, 6473)

#only 20 plots chosen
# extract subset from veg_plots original data
veg_plot_20_test <- subset(veg_plots, veg_plots@data$FLATE_NR %in% test_20_plots)

writeOGR(obj=veg_plot_20_test, layer = 'AR_subset_20plots', dsn=paste0(out_folder,"AR18x18"),
         overwrite_layer = TRUE, driver = 'ESRI Shapefile')

# Procedure for creating 1km square shapefiles around subset of AR18x18 centroids
center_20_test <- gCentroid(veg_plot_20_test,byid=TRUE, id = veg_plot_20_test@data$FLATE_NR)
ID=veg_plot_20_test@data$FLATE_NR
center_20_test_df <- SpatialPointsDataFrame(center_20_test,data.frame(id=ID, row.names=ID))
writeOGR(obj=center_20_test_df, layer = 'AR_centroid_20plots', dsn=paste0(out_folder,"AR18x18"),
         overwrite_layer = TRUE, driver = 'ESRI Shapefile')
# plot(center_20_test, add=TRUE)
# plot(veg_plots, col="red", add=TRUE)

# set the radius for the plots
radius <- 500 # radius in meters

# define the plot edges based upon the plot radius. 
yPlus <- center_20_test@coords[,2] + radius
xPlus <- center_20_test@coords[,1]+radius
yMinus <- center_20_test@coords[,2]-radius
xMinus <- center_20_test@coords[,1]-radius

# calculate polygon coordinates for each plot centroid. 
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus,yMinus,  # SE corner
             xMinus,yMinus, # SW corner
             xMinus,yPlus)  # NW corner again - close ploygon

# Extract the plot ID information
ID=veg_plot_20_test@data$FLATE_NR


# First, initialize a list that will later be populated
# a, as a placeholder, since this is temporary
a <- vector('list', length(2))

# loop through each centroid value and create a polygon
# this is where we match the ID to the new plot coordinates
for (i in 1:length(center_20_test@coords[,1])) {  # for each for in object centroids
  a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), ID[i]) 
  # make it an Polygon object with the Plot_ID from object ID
}

# convert a to SpatialPolygon and assign CRS
single_cell_1km<-SpatialPolygons(a,proj4string=CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Create SpatialPolygonDataFrame -- this step is required to output multiple polygons.
single_cell_1km_df <- SpatialPolygonsDataFrame(single_cell_1km, data.frame(id=ID, row.names=ID))
#single_cell_1km_df <- readOGR(layer = 'single_cell_1km_df', out_folder)
crs(single_cell_1km_df)
#crs(single_cell_1km_df) <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# write 1x1 into the shapefiles 
writeOGR(single_cell_1km_df, layer = 'single_cell_1km_df', out_folder,
         overwrite_layer = TRUE, driver = 'ESRI Shapefile')

#### Create 1km boxes 1081x####
#here, since DGVM is not available, we can create boxes corresponding to AR18x18 (600*1500)
# or in fact we can just crop with "Hui_dissolved_plots"
veg_plots <- readOGR(in_veg_data, "Hui_dissolved_plots")
center_1081 <- gCentroid(veg_plots,byid=TRUE, id = veg_plots@data$FLATE_NR)
radius <- 500 # radius in meters
# define the plot edges based upon the plot radius. 
yPlus <- center_1081@coords[,2]+radius
xPlus <- center_1081@coords[,1]+radius
yMinus <- center_1081@coords[,2]-radius
xMinus <- center_1081@coords[,1]-radius
# calculate polygon coordinates for each plot centroid. 
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus,yMinus,  # SE corner
             xMinus,yMinus, # SW corner
             xMinus,yPlus)  # NW corner again - close ploygon
# Extract the plot ID information
ID=veg_plots@data$FLATE_NR
# First, initialize a list that will later be populated
# a, as a placeholder, since this is temporary
a <- vector('list', length(2))
# loop through each centroid value and create a polygon
# this is where we match the ID to the new plot coordinates
for (i in 1:length(center_1081@coords[,1])) {  # for each for in object centroids
  a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), ID[i]) 
  # make it an Polygon object with the Plot_ID from object ID
}
# convert a to SpatialPolygon and assign CRS
single_cell_1km_1081<-SpatialPolygons(a,proj4string=CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Create SpatialPolygonDataFrame -- this step is required to output multiple polygons.
single_cell_1km_1081_df <- SpatialPolygonsDataFrame(single_cell_1km_1081, data.frame(id=ID, row.names=ID))
#single_cell_1km_df <- readOGR(layer = 'single_cell_1km_df', out_folder)
crs(single_cell_1km_1081_df)
#crs(single_cell_1km_df) <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# write 1x1 into the shapefiles 
writeOGR(single_cell_1km_1081_df, layer = 'single_cell_1km_1081_df', paste0(out_folder, "/single_cells_1081_compare"),
         overwrite_layer = TRUE, driver = 'ESRI Shapefile')

## frequency TEMP & PRECIP distribution diagrams ####
#20 vs 1081 plots along temp and precip gradient 
# load in rasters for temp and precip
bioclim_1 <- raster(paste0(in_raster,"bioclim_1.tif"))
bioclim_12 <- raster(paste0(in_raster,"bioclim_12.tif"))
plot(bioclim_1)
plot(center_1081, add=TRUE)
# extract values
c20_bioclim_1 <- extract(bioclim_1, center_20_test)
c1081_bioclim_1 <- extract(bioclim_1, center_1081)
c20_bioclim_12 <- extract(bioclim_12, center_20_test)
c1081_bioclim_12 <- extract(bioclim_12, center_1081)
# plot histograms 
hist(c20_bioclim_1)
hist(c1081_bioclim_1)
hist(c20_bioclim_12)
hist(c1081_bioclim_12)
# various tests
var.test(c20_bioclim_1, c1081_bioclim_1)
t.test(c20_bioclim_1, c1081_bioclim_1)
ks.test(c20_bioclim_1, c1081_bioclim_1)

var.test(c20_bioclim_12, c1081_bioclim_12)
t.test(c20_bioclim_12, c1081_bioclim_12)
ks.test(c20_bioclim_12, c1081_bioclim_12)
chisq.test(c20_bioclim_12, c1081_bioclim_12)

# publication plots in ggplot
temperature_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_bioclim_1,c1081_bioclim_1))
precipitation_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_bioclim_12,c1081_bioclim_12))
library(ggplot2)
ggplot(temperature_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
ggplot(precipitation_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
library(plyr)
t_mean <- ddply(temperature_data, "type", summarise, rating.mean=mean(data, na.rm=T))
t_mean
p_mean <- ddply(precipitation_data, "type", summarise, rating.mean=mean(data, na.rm=T))
p_mean

t_plot <- ggplot(temperature_data, aes(x=data, fill=type)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_density(alpha=.3) + 
  geom_vline(data = t_mean, aes(xintercept = rating.mean, colour = type),
          linetype = "longdash", size=1) +
  labs(title="frequency distribution", 
       #caption="Source: NCAR, SatVeg",
       x="Annual Mean Temperature (°C)",
       y="Density")

p_plot <- ggplot(precipitation_data, aes(x=data, fill=type)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_density(alpha=.3) + 
  geom_vline(data = p_mean, aes(xintercept = rating.mean, colour = type),
             linetype = "longdash", size=1) +
  labs(title="frequency distribution", 
       #caption="Source: NCAR, SatVeg",
       x="Annual Precipitation (mm)",
       y="Density")

#export
ggsave(plot = t_plot, filename = paste0("Temp_distribution_20vs1081",".png"), device = "png", path = paste0(out_ggplots),dpi = 300 )
ggsave(plot = p_plot, filename = paste0("Precip_distribution_20vs1081",".png"), device = "png", path = paste0(out_ggplots),dpi = 300 )

#### climate data test####
# cordex vs seNorge
cordex_check <-read.csv2("C:/Users/peterhor/Google Drive/2015 - University of Oslo/Project 2/20_plots_lat_long_alt_temp_precip_CORDEX.csv", sep = ";", dec = ".")
plot(cordex_check$SeNorge_precip~cordex_check$CORDEX_precip)
abline(a=1, b=1)
plot(cordex_check$SeNorge_temp~cordex_check$CORDEX_temp)
abline(a=1, b=1)

library(ggplot2)

climate_precip <- ggplot(cordex_check, aes(x=SeNorge_precip, y=CORDEX_precip)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_point(data = cordex_check, size=2) +
  geom_abline(intercept =  , slope = 1, color="red", 
              linetype="dashed", size=0.5) +
  labs(title="CORDEX vs SeNorge", 
       #caption="Source: NCAR, SatVeg",
       x="SeNorge - Annual Precipitation (mm)",
       y="CORDEX - Annual Precipitation (mm)")
climate_precip
ggsave(plot = climate_precip, filename = paste0("CORDEX_vs_seNorge_precip.png"), device = "png", width = 10, height = 7, path = paste0(path),dpi = 300 )

climate_temp <- ggplot(cordex_check, aes(x=SeNorge_temp, y=CORDEX_temp)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_point(data = cordex_check, size=2) +
  geom_abline(intercept =  , slope = 1, color="red", 
              linetype="dashed", size=0.5) +
labs(title="CORDEX vs SeNorge", 
     #caption="Source: NCAR, SatVeg",
     x="SeNorge - Annual Mean Temperature (°C)",
     y="CORDEX - Annual Mean Temperature (°C)")
climate_temp

ggsave(plot = climate_temp, filename = paste0("CORDEX_vs_seNorge_temp.png"), device = "png", width = 10, height = 7, path = paste0(path),dpi = 300 )

### ******VT*****#####
#####....###############################
#### ** AR dataset ####
#AR18x18 Validation dataset
ar18x18_raw <- readOGR(in_veg_data, "n18x18")
# transform to the CRS of our choice
#ar18x18_raw <- spTransform(ar18x18_raw, crs("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ar18x18_raw@proj4string

#create summary with xtabs, saving as matrix
unique(ar18x18_raw@data$VEG1) # 61 unique VEG
# Clean dataset - get rid of irrelevant VEG
#ar18x18_raw <- ar18x18_raw[!ar18x18_raw@data$VEG1 %in% c("FI", "Fin", "Hav", "RU", "Sv", "Swe", "SWE"),]
AR_VT_1081_polygons <- subset(ar18x18_raw, !(VEG1 %in% c("FI", "Fin", "Hav", "RU", "Sv", "Swe", "SWE")))
AR_VT_1081 <- as.data.frame.matrix(xtabs(AR_VT_1081_polygons$AREA~AR_VT_1081_polygons$VEG1+AR_VT_1081_polygons$FLATE_NR))
row.names(AR_VT_1081) # wrong order # noquote(row.names(tab_21_raw))
#adjust order of rows so that not 10a but 1a comes first
AR_VT_1081 <- AR_VT_1081[match(transl_rule_AR[,1], rownames(AR_VT_1081)),]
row.names(AR_VT_1081) # changed order
# adjust to percentages per plot
ARtemp_df <- list()
ARtemp <- sapply(names(AR_VT_1081), function(x) {
  ARtemp_df[paste0(x, "_pct")] <<- (AR_VT_1081[x] / sum(AR_VT_1081[x]))*100
  })
AR_VT_1081_perc<-as.data.frame(ARtemp)
rownames(AR_VT_1081_perc) <- rownames(AR_VT_1081)
AR_VT_1081_perc <- round(AR_VT_1081_perc, 2)
names(AR_VT_1081_perc) <- paste0("X", names(AR_VT_1081))
#AR_VT_1081_perc[1:4]
write.csv2(AR_VT_1081_perc, paste0(out_folder,"/AR18x18/AR_VT_1081.csv"))

#filter out relevant plots for this study
AR_VT_20_polygons <- subset(AR_VT_1081_polygons, AR_VT_1081_polygons@data$FLATE_NR %in% test_20_plots)
#create summary with xtabs, saving as matrix
AR_VT_20 <- as.data.frame.matrix(xtabs(AR_VT_20_polygons$AREA~AR_VT_20_polygons$VEG1+AR_VT_20_polygons$FLATE_NR))
#adjust order of rows so that not 10a but 1a comes first
AR_VT_20 <- AR_VT_20[match(transl_rule_AR[,1], rownames(AR_VT_20)),]
row.names(AR_VT_20) # changed order
ARtemp_df <- list()
ARtemp <- sapply(names(AR_VT_20), function(x) {
  ARtemp_df[x] <<- (AR_VT_20[x] / sum(AR_VT_20[x]))*100
})
AR_VT_20_perc<-as.data.frame(ARtemp)
rownames(AR_VT_20_perc) <- rownames(AR_VT_20)
AR_VT_20_perc <- round(AR_VT_20_perc, 2)
names(AR_VT_20_perc) <- paste0("X", colnames(AR_VT_20))
#AR_VT_20_perc[1:4]
write.csv2(AR_VT_20_perc, paste0(out_folder,"/AR18x18/AR_VT_20.csv"))

#### test representativness of 20 plots ####
#sum areas
AR_VT_1081_freq <- rowSums(AR_VT_1081/sum(AR_VT_1081)*100) #
write.csv2(AR_VT_1081_freq, paste0(out_folder,"/AR18x18/AR_VT_1081_freq.csv"))
AR_VT_20_freq <- rowSums(AR_VT_20/sum(AR_VT_20)*100)
write.csv2(AR_VT_20_freq, paste0(out_folder,"/AR18x18/AR_VT_20_freq.csv"))

#AR_stats_20/AR_stats_1081
#AR_perc_20/AR_perc_1081
AR_VT_repre <- cbind(as.data.frame(AR_VT_1081_freq), as.data.frame(AR_VT_20_freq))
var.test(AR_VT_1081_freq, AR_VT_20_freq)
t.test(AR_VT_1081_freq, AR_VT_20_freq)
ks.test(AR_VT_1081_freq, AR_VT_20_freq)
cor.test(AR_VT_1081_freq, AR_VT_20_freq)
chisq.test(AR_VT_1081_freq, AR_VT_20_freq)

write.csv2(AR_VT_repre, paste0(out_folder,"/AR18x18/AR_VT_representative.csv"))
plot(AR_VT_1081_freq~AR_VT_20_freq)


# problem with self-intersection solved by following https://gis.stackexchange.com/questions/119474/clipping-two-spatialpolygons-error-in-rgeosbintopofunc
ar18x18_test <-gBuffer(ar18x18_test, byid=TRUE, width=0)
single_cell_1km_df <- gBuffer(single_cell_1km_df, byid=TRUE, width=0)
ar_clip <- crop(ar18x18_test, single_cell_1km, byid=TRUE)
#update areas of polygons
ar_clip$AREA_2 <- area(ar_clip)
    # head(clip_test@data[,c("AREA_2","AREA")], n=100)
    # head(sort(clip_test@data[,"AREA_2"]), n=10)
    # clip_test@data$AREA[1:40]
    #clip_test@data$AREA_2[1:40]
plot(subset(ar_clip, ar_clip@data$FLATE_NR==222))

ar_clip_data <- data.table(ar_clip@data[,c("VEG1","FLATE_NR","AREA_2")])# save only relevant data as dataframe
# examine coverage for each FLATE_NR
# area stats for each plot
#create summary with xtabs, saving as matrix
ar_clip_vt_flate <- as.data.frame.matrix(xtabs(ar_clip_data$AREA~ar_clip_data$VEG1+ar_clip_data$FLATE_NR))
row.names(ar_clip_vt_flate) # noquote(row.names(tab_21_raw))
#adjust order of rows so that not 10a but 1a comes first
ar_clip_vt_flate <- ar_clip_vt_flate[match(transl_rule_AR[,1], rownames(ar_clip_vt_flate)),]
row.names(ar_clip_vt_flate) # changed order

#### .... #####################
#### ** RS dataset ####
#### reclassify RS map of VT's into PFT map
#load in translation rules
transl_rule_RS <- read.csv(paste0(in_csv_files, "02_Translation_schemes/RS_to_PFT_final.csv"), sep = ";", dec = ".")
transl_RS <- matrix(nrow= 27, ncol = 2)
# make into reclassify matrix # can also to like this cbind(transl_rule_DM$VT_ID, transl_rule_DM$PFT_ID)
for(i in 1:NROW(transl_rule_RS)){
  transl_RS[i,] <- c(transl_rule_RS$RS.code[i], transl_rule_RS$PFT_ID[i])
  # print(transl)
}
#transl <- cbind(transl_rule_RS$RS_code, transl_rule_RS$PFT_ID)

RS_Norway_PFT <- reclassify(RS_Norway33, transl_RS, include.lowest=TRUE, overwrite=TRUE,
                            filename = paste0(out_folder,"/RS/RS_PFT.tif")) 

RS_freq <- freq(RS_Norway33)
write.csv2(RS_freq, paste0(out_folder,"/RS/RS_freq.csv"))
RS_PFT_freq <- freq(RS_Norway_PFT)
write.csv2(RS_PFT_freq, file = paste0(out_folder,"RS/RS_PFT_freq.csv"))
#writeRaster(RS_Norway_PFT, filename = paste0(out_folder,"/RS/RS_Norway_PFT.tif"), format="GTiff", overwrite=TRUE)

# then mask out only 1x1km from raster (this takes 15min)
RS_PFT_single <- mask(RS_Norway_PFT, single_cell_1km_df, filename=paste0(out_folder,"/RS_PFT_single_cells_R"), format="GTiff", overwrite=TRUE)
#RS_S_freq <- freq(RS_south)

# extract values for all 1km squares
single_cells_extr <- extract(RS_PFT_single, single_cell_1km_df, df=TRUE)
  # single_cells_extr <- extract(RS_single_cells, single_cell_1km_df)
  # names(single_cells_extr) <- test_20_plots

sort(unique(single_cells_extr$RS_PFT_cells_R)) # missing representation of RS_type 15 and 24
RS_table <- as.data.frame.matrix(t(table(single_cells_extr, useNA = "ifany"))) # transpose columns and rows and save summarized table, 
colnames(RS_table) <- test_20_plots # assign column names with plot ID
colSums(RS_table)

# table(subset(single_cells_extr,single_cells_extr$ID==21)) # output only one locality *(ID=1)
write.csv2(RS_table, paste0(out_folder,"/RS_table.csv"))


#### .... #####################     
#### ** DM dataset ####
# load all rasters into a stack
# read all raster files inside HOME folder and add them to a list
r.list <- list.files(in_cropped_data_DM, pattern="tif$", full.names=TRUE)
r.list <- mixedsort(r.list)
r.stack <- stack(r.list)
# r.list <- list.files(in_probab_raster, pattern="tif$", full.names=FALSE)
# r.list.path <- list.files(in_probab_raster, pattern="tif$", full.names=TRUE)
names(r.stack)
basename(r.list)
crs(r.stack)

# all norway
time_pred_start <- Sys.time()
max_val <- which.max(r.stack)
writeRaster(max_val, filename = paste0(out_folder_maxval, "DM_maxval.tif"), format="GTiff", overwrite=TRUE)
# max_val <- raster(paste0(out_folder_maxval, "DM_maxval.tif"))
# max_val <- raster(paste0(in_dm_data,"DMtest_whole_WGS_UTM32.tif"))
time_pred_end <- Sys.time()
pred_time <- (time_pred_end-time_pred_start)
print( paste( "calculation finished in", pred_time, sep=" "))
plot(max_val)
# r.freq <- table(getValues(max_val))
raster_freq <- freq(max_val)
write.csv2(raster_freq, paste0(out_folder_maxval,"DM_maxval_freq.csv"))
max_val <- as.factor(max_val)
raster_freq_names <- cbind(raster_freq, vt_names)
plot(max_val)
table(getValues(max_val))

### plot into file#
library(RColorBrewer)
#hi-res PDF into default folder
pdf(file = paste0(out_folder_maxval,"DM_maxVal.pdf"))
plot(max_val) #c("grey",rev(rainbow(28)))); ,col= brewer.pal(28, "Spectral")
title(main = paste("Maximum Value - unique category"))
dev.off()

#low res TIFF
tiff(filename = paste0(out_folder_maxval,"/DM_maxVal_lowres.tiff",
                       height = 27, width = 17, units = 'cm', res=300))
plot(max_val,col= c("grey",rev(rainbow(9))))
title(main = paste("Maximum Value - unique category"))
dev.off()

#### .... #####################     
## function RASTER to PFT table ####
rast_to_table <- function(rast, plots){
  # mask out only 1x1km from raster (this takes 15min)
  rast_mask <- mask(rast, plots)
  # extract values for all 1km squares
  print("finished masking")
  rast_mask_df <- extract(rast_mask, plots, df=TRUE)
  # make table
  print("finished extracting")
  rast_mask_table <- as.data.frame.matrix(t(table(rast_mask_df))) # transpose columns and rows and save summarized table, 
  rast_mask_table_excl <- rast_mask_table[7,]
  rast_mask_table <- rast_mask_table[1:6,]
  # divide by sum of each column (number of raster cells) * 100
  temp_df <- list()
  ARtemp <- sapply(colnames(rast_mask_table), function(x) {
    temp_df[paste0(x)] <<- (rast_mask_table[x] / sum(rast_mask_table[x]))*100
  })
  rast_mask_table<-as.data.frame(ARtemp)
  # round to 2 decimals
  # rast_mask_table <- round(rast_mask_table, 2)
  rownames(rast_mask_table) <- PFT[1:6,2]
  colnames(rast_mask_table) <- paste0("X", plots@data$id) # assign column names with plot ID
  # to decide between 20 and 1081 dataset (for ease of saving reasons)
  if( nrow(plots@data)<30){
    print("saving 20 plots")
    subfolder=20} else {
        print("saving 1081 plots")
        subfolder=1081}
  #check size of plots to distinguish between 1km and 0.9km
  if( area(plots)[1]==1e+06){
    print("saving 1km squares")
    subtype="single_cells_"
    filetype="_sc"
    colnames(rast_mask_table) <- paste0("X", plots@data$id)
    } else {
      print("saving AR rectangles")
      subtype="single_AR_plots_"
      filetype="_ar"
      colnames(rast_mask_table) <- paste0("X", plots@data$FLATE_NR)}
  # table(subset(single_cells_extr,single_cells_extr$ID==21)) # output only one locality *(ID=1)
    write.csv2(rast_mask_table, paste0(out_folder,subtype, subfolder,"_compare/",names(rast), filetype, subfolder,".csv"))
    write.csv2(rast_mask_table_excl, paste0(out_folder, subtype, subfolder,"_compare/",names(rast), filetype, subfolder,"_excl.csv"))
      # Save output
    saveRDS(rast_mask_table, paste0(out_folder,subtype, subfolder,"_compare/",names(rast), filetype, subfolder))
    writeRaster(rast_mask,filename=paste0(out_folder,subtype, subfolder,"_compare/", names(rast), filetype, subfolder), format="GTiff", overwrite=TRUE)
  
  #return(rast_mask_table)
}
# executing rast_to_table####
#load rasters for use in function

RS_PFT <- raster(paste0(out_folder,"/RS/RS_PFT.tif"))
DM_AUC_PFT <- raster(paste0(out_folder_AUC,"DM_AUC_PFT.tif"))
DM_maxval_PFT <- raster(paste0(out_folder_maxval,"DM_maxval_PFT.tif"))
DM_preval_PFT <- raster(paste0(out_folder_preval,"DM_preval_PFT.tif"))
# execute function 20squares
rast_to_table(RS_PFT, single_cell_1km_df)
rast_to_table(DM_AUC_PFT, single_cell_1km_df)
rast_to_table(DM_maxval_PFT, single_cell_1km_df)
rast_to_table(DM_preval_PFT, single_cell_1km_df)
# execute function 1081 squares
rast_to_table(RS_PFT, single_cell_1km_1081_df)
rast_to_table(DM_AUC_PFT, single_cell_1km_1081_df)
rast_to_table(DM_maxval_PFT, single_cell_1km_1081_df)
rast_to_table(DM_preval_PFT, single_cell_1km_1081_df)
# execute function for 20 plots
rast_to_table(RS_PFT, veg_plot_20_test)
rast_to_table(DM_AUC_PFT, veg_plot_20_test)
rast_to_table(DM_maxval_PFT, veg_plot_20_test)
rast_to_table(DM_preval_PFT, veg_plot_20_test)
# execute function for 1081 plots
rast_to_table(RS_PFT, veg_plots)
rast_to_table(DM_AUC_PFT, veg_plots)
rast_to_table(DM_maxval_PFT, veg_plots)
rast_to_table(DM_preval_PFT, veg_plots)


# execute funtion for nonPFT rasters
#rast_to_table(RS_Norway33, single_cell_1km_df)
# rast_to_table(DM_AUC, single_cell_1km_df)
# rast_to_table(DM_maxval, single_cell_1km_df)
# rast_to_table(DM_preval, single_cell_1km_df)



#### .... ##################### 
#### TRANSLATION SCHEME ####
#### ******PFT*****#####
# prepare or load data
# DM had to be prepared in EXCEL because output classification from models didn't match with 31 VTs *(only 19)
DM <- read.csv(paste0(in_csv_files,"01_Comparison_schemes/DM_maxval_PFT.csv"), sep = ";", dec = ".")
      colnames(DM)[1] <- c("VT_code") #, test_20_plots
DGVM <- read.csv(paste0(in_csv_files,"01_Comparison_schemes/DGVM_PFT.csv"), sep = ";", dec = ".")
AR <- read.csv(paste0(in_csv_files,"01_Comparison_schemes/AR18x18_PFT.csv"), sep = ";", dec = ".")
RS <- read.csv(paste0(in_csv_files,"01_Comparison_schemes/RS_PFT.csv"), sep = ";", dec = ".")
      #colnames(RS) <- c("RS_code", test_20_plots)
      
# THERE IS A PROBLEM WITH MISSING VALUES IN RS DATASET *(15 AND 24 ARE NOT REPRESENTED) 
#NEED TO ADD MANUALLY IN CSV
#reuse RS from previous section
names(RS)[1] <- "RS_code"
RS_compare <- RS
 
library(dplyr)
  #load in translation rules
transl_rule_DM <- read.csv(paste0(in_csv_files, "02_Translation_schemes/VEG_to_PFT_final.csv"), sep = ";", dec = ".")
transl_rule_RS <- read.csv(paste0(in_csv_files, "02_Translation_schemes/RS_to_PFT_final.csv"), sep = ";", dec = ".")
transl_rule_AR <- read.csv(paste0(in_csv_files, "02_Translation_schemes/AR18x18_to_PFT_final.csv"), sep = ";", dec = ".")

# * DM ####
DM_compare <- data.frame(append(DM, list(transl_rule_DM$PFT_code), after=match("VT_code", names(DM))))
DM_compare <- data.frame(append(DM_compare, list(transl_rule_DM$PFT_name), after=match("VT_code", names(DM))))
names(DM_compare)[3] <- "PFT_name"
names(DM_compare)[4] <- "PFT_code"
DM_pft <- aggregate(DM_compare[5:length(DM_compare)], by=list(PFT_code=DM_compare$PFT_code), FUN=sum)
#after being aggregated, need to harmonize sums to fractions of 100 (%)
DM_pft_100 <- DM_pft
str(DM_pft)
# divide by sum of each column (number of raster cells) * 100
DM_pft_100[2:NCOL(DM_pft_100)] <- DM_pft_100[2:NCOL(DM_pft_100)]/colSums(DM_pft[2:NCOL(DM_pft)])*100
# round to 1 decimal
DM_pft_100[2:NCOL(DM_pft_100)] <- round(DM_pft_100[2:NCOL(DM_pft_100)], 1)

# * RS ####
RS_compare <- data.frame(append(RS, list(rownames(RS)), after=0))
names(RS_compare)[1] <- "RS_code"
RS_compare <- data.frame(append(RS_compare, list(transl_rule_RS$PFT_code), after=match("RS_code", names(RS_compare))))
RS_compare <- data.frame(append(RS_compare, list(transl_rule_RS$PFT_name), after=match("RS_code", names(RS))))
names(RS_compare)[2] <- "PFT_name"
names(RS_compare)[3] <- "PFT_code"
RS_pft <- aggregate(RS_compare[4:length(RS_compare)], by=list(PFT_code=RS_compare$PFT_code), FUN=sum)
lapply(RS_pft[2:22], sum)
colSums(RS_pft[2:22])

#after being aggregated, need to harmonize sums to fractions of 100 (%)
RS_pft_100 <- RS_pft
  # divide by sum of each column (number of raster cells) * 100
RS_pft_100[2:NCOL(RS_pft_100)] <- RS_pft[2:NCOL(RS_pft)]/colSums(RS_pft[2:NCOL(RS_pft)])*100
  # round to 1 decimal
RS_pft_100[2:NCOL(RS_pft_100)] <- round(RS_pft_100[2:NCOL(RS_pft_100)], 1)

# * AR ####
# first reclassify the shapefile into PFT for later use in QGIS

#save data into temp datatable
#tmp_ar <- data.table(AR_VT_1081_polygons@data)
#merge based on VT
unique(ar18x18_raw$VEG1) 
unique(transl_rule_AR$VT_code)
# always pass the merge command your Spatial*DataFrame object.
AR_PFT_1081_polygons <- merge(ar18x18_raw , transl_rule_AR, by.x = "VEG1",by.y = "VT_code")
#save the file
writeOGR(AR_PFT_1081_polygons, layer =  "AR_PFT_1081", paste0(out_folder, "AR18x18"), overwrite_layer = TRUE, driver="ESRI Shapefile")
# cut out the 20 plots
AR_PFT_20_polygons <- subset(AR_PFT_1081_polygons, AR_PFT_1081_polygons@data$FLATE_NR %in% test_20_plots)
writeOGR(AR_PFT_20_polygons, layer =  "AR_PFT_20", paste0(out_folder, "AR18x18"), overwrite_layer = TRUE, driver="ESRI Shapefile")

# add column with row names equal to first column
AR_VT_20rownames <- data.frame(append(AR_VT_20, list(rownames(AR_VT_20)), after=0))
names(AR_VT_20rownames)[1] <- "VT_code"
# append PFT code and name for TRANSLATION 
AR_VT_20rownames <- data.frame(append(AR_VT_20rownames, list(transl_rule_AR$PFT_code), after=match("VT_code", names(AR_VT_20rownames))))
AR_VT_20rownames <- data.frame(append(AR_VT_20rownames, list(transl_rule_AR$PFT_name), after=match("VT_code", names(AR_VT_20rownames))))
names(AR_VT_20rownames)[2] <- "PFT_name"
names(AR_VT_20rownames)[3] <- "PFT_code"
#colnames(AR_VT_20rownames)[4:ncol(AR_VT_20rownames)] <- colnames(AR_VT_20)
write.csv2(AR_VT_20rownames, paste0(out_folder,"/AR18x18/AR_VT_20rownames.csv"))
# aggregate according to translation scheme
AR_PFT_20 <- aggregate(AR_VT_20rownames[4:length(AR_VT_20rownames)], by=list(PFT_code=AR_VT_20rownames$PFT_code), FUN=sum)
# changed order
rownames(AR_PFT_20) <- AR_PFT_20[,1]
AR_PFT_20 <- AR_PFT_20[2:ncol(AR_PFT_20)]
AR_PFT_20 <- AR_PFT_20[match(PFT[,2], rownames(AR_PFT_20)),]
write.csv2(AR_PFT_20, paste0(out_folder,"/AR18x18/AR_PFT_20.csv"))
#colnames(AR_PFT_20) <- colnames(AR_VT_20)
#after being aggregated, need to harmonize sums to fractions of 100 (%)
# divide by sum of each column (number of raster cells) * 100
ARtemp_df <- list()
ARtemp <- sapply(names(AR_PFT_20), function(x) {
  ARtemp_df[paste0(x)] <<- (AR_PFT_20[x] / sum(AR_PFT_20[x]))*100
})
AR_PFT_20_perc<-as.data.frame(ARtemp)
rownames(AR_PFT_20_perc) <- rownames(AR_PFT_20)
AR_PFT_20_perc <- round(AR_PFT_20_perc, 2)
names(AR_PFT_20_perc) <- paste0("X", names(AR_VT_20))
#AR_VT_1081_perc[1:4]
# we are excluding excluded types
write.csv2(AR_PFT_20_perc[7,], paste0(out_folder,"/single_cells_20_compare/AR_PFT_20_excl.csv"))
# recalculate the percentages again without the excluded 7th line

AR_PFT_20_perc <- AR_PFT_20_perc[-7,]
ARtemp_df <- list()
ARtemp <- sapply(names(AR_PFT_20_perc), function(x) {
  ARtemp_df[paste0(x)] <<- (AR_PFT_20_perc[x] / sum(AR_PFT_20_perc[x]))*100
})
AR_PFT_20_perc<-as.data.frame(ARtemp)
rownames(AR_PFT_20_perc) <- rownames(AR_PFT_20)[1:6]
#AR_PFT_20_perc <- round(AR_PFT_20_perc, 2)
names(AR_PFT_20_perc) <- paste0("X", names(AR_VT_20))
write.csv2(AR_PFT_20_perc, paste0(out_folder,"/single_cells_20_compare/AR_PFT_20_perc.csv"))


AR_PFT_20_freq <- rowSums(AR_PFT_20/sum(AR_PFT_20)*100)
write.csv2(AR_PFT_20_freq, paste0(out_folder,"/AR18x18/AR_PFT_20_freq.csv"))


#### test representativenes of PFTs 
# repeat the steps above with the whole 1081 dataset polygons
#ar18x18_vt_flate <- ar18x18_raw_str
AR_VT_1081rownames <- data.frame(append(AR_VT_1081, list(rownames(AR_VT_1081)), after=0))
names(AR_VT_1081rownames)[1] <- "VT_code"
# append PFT code and name for TRANSLATION 
AR_VT_1081rownames <- data.frame(append(AR_VT_1081rownames, list(transl_rule_AR$PFT_code), after=match("VT_code", names(AR_VT_1081rownames))))
AR_VT_1081rownames <- data.frame(append(AR_VT_1081rownames, list(transl_rule_AR$PFT_name), after=match("VT_code", names(AR_VT_1081rownames))))
names(AR_VT_1081rownames)[2] <- "PFT_name"
names(AR_VT_1081rownames)[3] <- "PFT_code"
write.csv2(AR_VT_1081rownames, paste0(out_folder,"/AR18x18/AR_VT_1081rownames.csv"))
# aggregate according to translation scheme
AR_PFT_1081 <- aggregate(AR_VT_1081rownames[4:length(AR_VT_1081rownames)], by=list(PFT_code=AR_VT_1081rownames$PFT_code), FUN=sum)
#after being aggregated, need to harmonize sums to fractions of 100 (%)
rownames(AR_PFT_1081) <- AR_PFT_1081[,1]
AR_PFT_1081 <- AR_PFT_1081[,2:ncol(AR_PFT_1081)]
AR_PFT_1081 <- AR_PFT_1081[match(PFT[,2], rownames(AR_PFT_1081)),]
write.csv2(AR_PFT_1081, paste0(out_folder,"/AR18x18/AR_PFT_1081.csv"))

ARtemp_df <- list()
ARtemp <- sapply(names(AR_PFT_1081), function(x) {
  ARtemp_df[paste0(x)] <<- (AR_PFT_1081[x] / sum(AR_PFT_1081[x]))*100
})
AR_PFT_1081_perc<-as.data.frame(ARtemp)
rownames(AR_PFT_1081_perc) <- rownames(AR_PFT_1081)
AR_PFT_1081_perc <- round(AR_PFT_1081_perc, 2)
names(AR_PFT_1081_perc) <- names(AR_PFT_1081)
#AR_PFT_1081_perc[1:4]
write.csv2(AR_PFT_1081_perc, paste0(out_folder,"/AR18x18/AR_PFT_1081_perc.csv"))
write.csv2(AR_PFT_1081_perc[7,], paste0(out_folder,"/AR18x18/AR_PFT_1081_perc_excl.csv"))
# excluding the excluded areas and recalculating proportions
AR_PFT_1081_perc6PFT <- AR_PFT_1081_perc[1:6,]
ARtemp_df <- list()
ARtemp <- sapply(names(AR_PFT_1081_perc6PFT), function(x) {
  ARtemp_df[paste0(x)] <<- (AR_PFT_1081_perc6PFT[x] / sum(AR_PFT_1081_perc6PFT[x]))*100
})
AR_PFT_1081_perc6PFT<-as.data.frame(ARtemp)
rownames(AR_PFT_1081_perc6PFT) <- rownames(AR_PFT_1081)[1:6]
AR_PFT_1081_perc6PFT <- round(AR_PFT_1081_perc6PFT, 2)
names(AR_PFT_1081_perc6PFT) <- names(AR_PFT_1081)

write.csv2(AR_PFT_1081_perc6PFT, paste0(out_folder,"/single_cells_1081_compare/AR_PFT_sc1081.csv"))


AR_PFT_1081_freq <- rowSums(AR_PFT_1081/sum(AR_PFT_1081)*100)
write.csv2(AR_PFT_1081_freq, paste0(out_folder,"/AR18x18/AR_PFT_1081_freq.csv"))

#### test representativness of PFT 20 plots ####
#sum areas
AR_PFT_1081_freq <- rowSums(AR_PFT_1081/sum(AR_PFT_1081)*100) #
write.csv2(AR_PFT_1081_freq, paste0(out_folder,"/AR18x18/AR_PFT_1081_freq.csv"))
AR_PFT_20_freq <- rowSums(AR_PFT_20/sum(AR_PFT_20)*100)
write.csv2(AR_PFT_20_freq, paste0(out_folder,"/AR18x18/AR_PFT_20_freq.csv"))


#AR_stats_20/AR_stats_1081
#AR_perc_20/AR_perc_1081
AR_PFT_repre <- cbind(as.data.frame(AR_PFT_1081_freq), as.data.frame(AR_PFT_20_freq))
write.csv2(AR_PFT_repre, paste0(out_folder,"/AR18x18/AR_PFT_representative1.csv"))
plot(AR_PFT_1081_freq~AR_PFT_20_freq)
var.test(AR_PFT_1081_freq, AR_PFT_20_freq)
t.test(AR_PFT_1081_freq, AR_PFT_20_freq)
ks.test(AR_PFT_1081_freq, AR_PFT_20_freq)
cor.test(AR_PFT_1081_freq, AR_PFT_20_freq)
chisq.test(AR_PFT_1081_freq, AR_PFT_20_freq)
prop.test(AR_PFT_1081_freq, AR_PFT_20_freq)

vegdist(t(cbind(AR_PFT_1081_freq, AR_PFT_20_freq)))

#************

pft_data <- PFT_df # PFT_df is created within 01_comparison_vegan.R
path="E:/Project_2_RUNS/OUTPUT/results"
pft_data <- read.csv2(file = paste0(path, "/PFT_df.csv"))
pft_data <- read.csv2(file = paste0(path, "/PFT_df_adj.csv"), dec = ".") # adjusted values for precipitation seasonality (bioclim_15) from DGVM sensitivity analysis

######### GRAPHICS #############
# One specific example of GGPLOT
g <- ggplot(pft_data, aes(x = Method, y = pft_data$X622, fill = PFT_name)) +
  colScale

g + geom_bar(position = "fill",stat = "identity") + # position (stack, dodge, fill)
    labs(title="PFT coverage", 
       subtitle= paste(names(pft_data)[3]),
       caption="Source: NCAR, SatVeg",
       x="Modelling Method",
       y="plot 622")

g <- ggplot(pft_data, aes(x = Method, y = X513, fill = PFT_code)) + scale_fill_brewer(palette = "Spectral") 
g 
 # scale_y_continuous(labels = percent_format())

# One specific example of GGPLOT split with GGarrange
g1 <- ggplot(pft_data[pft_data$Method %in% c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM"),], aes(x = Method, y = X405, fill = PFT_code)) +
  scale_fill_brewer(palette = "Spectral") 
g1 <- g1 + geom_bar(position = "fill",stat = "identity") + # position (stack, dodge, fill)
  labs(title="PFT coverage", 
       #caption="Source: NCAR, SatVeg",
       x="Modelling Method",
       y="plot 222") +
  theme(legend.position = "none") # hide legend
g2  <- ggplot(pft_data[pft_data$Method %in% c("AR18x18"),], aes(x = Method, y = X405, fill = PFT_code)) +
  scale_fill_brewer(palette = "Spectral") 
g2 <- g2 + geom_bar(position = "fill",stat = "identity") + # position (stack, dodge, fill)
  labs(title=" ", 
       x="Ground Survey", 
       y= NULL) + 
  scale_y_continuous(breaks = NULL)
#arrange next to each other
grid.arrange(g1,g2, widths = c(3,1))

g <- ggplot(pft_data, aes(x = Method, y = X513, fill = PFT_code)) + scale_fill_brewer(palette = "Spectral") 
g + geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())


# * generic loop for all ####
colNames <- names(pft_data)[5:NCOL(pft_data)] 
#colNames <- names(pft_data)[7:NCOL(pft_data)] # could be also 7:12 in adjusted version of DGVM
for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  g <- ggplot(pft_data, aes_string(x = "Method", y = i, fill = "PFT_name")) + scale_fill_brewer(palette = "Spectral") +
    geom_bar(position = "fill",stat = "identity") +
    scale_y_continuous(labels = percent_format()) + 
    labs(title="PFT coverage", 
         subtitle= paste(i),
         caption="Source: NCAR, SatVeg , NIBIO",
         x="Modelling Method",
         y= paste("plot",i))
  # print(g)
  # output png
  ggsave(filename = paste0(i, "DM3_DGVM_RS_AR.png"), device = "png", path = paste0(out_ggplots, "/plots_6way/"),dpi = 300 )
 }
# * select only maxval method for DM ####
for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  g <- ggplot(pft_data[pft_data$Method %in% c("AR18x18", "DM_maxval", "RS", "DGVM"),], aes_string(x = "Method", y = i, fill = "PFT_name"))+
    scale_fill_brewer(palette = "Spectral") +
    geom_bar(position = "fill",stat = "identity") +
    scale_y_continuous(labels = percent_format()) + 
    labs(title="PFT coverage", 
         subtitle= paste(i),
         caption="Source: NCAR, SatVeg , NIBIO",
         x="Modelling Method",
         y= paste("plot",i))
  # print(g)
  # output png
  ggsave(filename = paste0(i, "DM_maxval_DGVM_RS_AR.png"), device = "png", path = paste0(out_ggplots, "/plots_4way/"),dpi = 300 )
}

# * select also DGVM_adj method for comparison ####
for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  g <- ggplot(pft_data[pft_data$Method %in% c("DGVM", "DGVM_adj","RS", "DM_maxval", "AR18x18"),], aes_string(x = "Method", y = i, fill = "PFT_code"))+
    colScale + 
    # scale_fill_brewer(palette = "Spectral") +
    geom_bar(position = "fill",stat = "identity") +
    #scale_y_continuous(labels = percent_format()) + 
    labs(title="PFT coverage", 
         subtitle= paste(i),
         #caption="Source: NCAR, SatVeg , NIBIO",
         x="Modelling Method",
         y= paste("plot",i))
  # print(g)
  # output png
  ggsave(filename = paste0(i, "DM_maxval_DGVM_RS_AR.png"), device = "png", path = paste0(out_ggplots, "/DGVM_adj/"),dpi = 300 )
}


# * Maxval method vs AR for DM ####
# and split AR from the rest
for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  g1 <- ggplot(pft_data[pft_data$Method %in% c("DM_maxval", "RS", "DGVM"),],
                aes_string(x = "Method", y = i, fill = "PFT_name")) + 
    scale_fill_brewer(palette = "Spectral") +
    geom_bar(position = "fill",stat = "identity") +
    scale_y_continuous(labels = percent_format()) + 
    labs(title="PFT coverage", 
         #caption="Source: NCAR, SatVeg",
         x="Modelling Method",
         y=paste("plot ",i)) +
    theme(legend.position = "none")# hide legend
  g2 <- ggplot(pft_data[pft_data$Method %in% c("AR18x18"),],
               aes_string(x = "Method", y = i, fill = "PFT_name")) + 
    scale_fill_brewer(palette = "Spectral") +
    geom_bar(position = "fill",stat = "identity") +
    labs(title=" ", 
         #caption="Source: NCAR, SatVeg",
         x="Ground Survey",
         y=NULL) +
           scale_y_continuous(breaks = NULL)
  
   #arrange next to each other
  g12 <- arrangeGrob(g1,g2, widths = c(1,1))
  
  ggsave(plot = g12, filename = paste0("AR_vs_Methods_DM_maxval",i,".png"), device = "png", width = 10, height = 7, path = paste0(out_ggplots, "/plots_4way/"),dpi = 300 )
}

# * graphs FOR PUBLICATION ON A MAP IN QGIS ####
# no legend, no labels
colNames<-colnames(pft_data)[7:ncol(pft_data)] #[5:24] #
colNames_number<-gsub("X", "#", colNames)

#Create a custom color scale (https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin)
library(RColorBrewer)
myColors <- brewer.pal(6,"Spectral")
names(myColors) <- c("Broadleaf deciduous boreal shrub",
                     "Broadleaf deciduous boreal tree",
                     "Broadleaf deciduous temperate tree",
                     "C3 grass",
                     "Needleleaf evergreen boreal tree",
                     "Bare ground / Not vegetated")
#names(myColors) <- order(levels(pft_data$PFT_code), c(6, 3,2,5,1,4))
#colScale <- scale_colour_manual(name = "PFT_code",values = myColors)
colScale <- scale_fill_manual(name = "PFT_name",values = myColors)


for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  # need to get an increasing integer from 1:20  to replace for colNames_number
  # j <- colNames_number[i] #
  g1 <- ggplot(pft_data[pft_data$Method %in% c("DM_maxval", "RS", "DGVM", "AR18x18"),],
               aes_string(x = "Method", y = i, fill = "PFT_name")) + 
    # scale_fill_brewer(palette = "Spectral") +
    scale_x_discrete(limits= c("DGVM","RS","DM_maxval", "AR18x18")) + # changing the order and the label names ##labels=c("DGVM","RS","DM")
    geom_bar(position = "fill",stat = "identity") +
    colScale + # adjusting colors to match
    #scale_y_continuous(labels = percent_format()) + 
    labs(#title=paste("plot ",i), 
         #caption="Source: NCAR, SatVeg",
         x="Modelling Method",
         y="PFT profile") +
    theme_void() + theme(legend.position="none")    # hide legend and labels on X axis  # theme_void()
  
  # g2 <- ggplot(pft_data[pft_data$Method %in% c("AR18x18"),],
  #              aes_string(x = "Method", y = i, fill = "PFT_name")) + 
  #   scale_fill_brewer(palette = "Spectral") +
  #   geom_bar(position = "fill",stat = "identity") +
  #   labs(title=" ", 
  #        #caption="Source: NCAR, SatVeg",
  #        x="Reference Dataset",
  #        y=NULL) +
  #   scale_y_continuous(breaks = NULL) +
  #   theme_void() + theme(legend.position="none")# hide legend
  # g2
  # #arrange next to each other
  # g12 <- arrangeGrob(g1,g2, widths = c(2,1))
  # #plot(g12)
  ggsave(plot = g1, filename = paste0("AR_vs_Methods_DM_maxval",i,".png"), device = "png", width = 10, height = 7, path = paste0(out_ggplots, "/plots_4way/"),dpi = 300 ) 
  #("AR_vs_Methods_DM_maxval",i,".png")
}

#* DGVM adjusted ####
for (i in colNames){
  print(paste0('creating graphs DGVM vs RS for plot #', i ))
  # need to get an increasing integer from 1:20  to replace for colNames_number
  # j <- colNames_number[i] #
  g1 <- ggplot(pft_data[pft_data$Method %in% c("DGVM", "DGVM_adj", "AR18x18"),],
               aes_string(x = "Method", y = i, fill = "PFT_name")) + 
    # scale_fill_brewer(palette = "Spectral") +
    scale_x_discrete(limits= c("DGVM","DGVM_adj", "AR18x18")) + # changing the order and the label names ##labels=c("DGVM","RS","DM")
    geom_bar(position = "fill",stat = "identity") +
    colScale + # adjusting colors to match
    #scale_y_continuous(labels = percent_format()) + 
    labs(#title=paste("plot ",i), 
      #caption="Source: NCAR, SatVeg",
      x="Modelling Method",
      y="PFT profile")  + 
    theme(legend.position="none")
    #theme_void()     # hide legend and labels on X axis  # theme_void()
  ggsave(plot = g1, filename = paste0("AR_vs_DGVM_adj",i,".png"), device = "png", width = 10, height = 7, path = paste0(out_ggplots, "/plots_4way/"),dpi = 300 ) 
  #("AR_vs_Methods_DM_maxval",i,".png")
}


### ....####
#### COMPARISON Scheme####
# in file 01_comparison_vegan.R