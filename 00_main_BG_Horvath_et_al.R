#### Improving the representation of high-latitude vegetation in Dynamic Global Vegetation Models ####
# Dynamical Vegetation model version CLM4.5-BGCDV
# updated: 21.10.2020
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
library(data.table) 
library(dplyr)
library(vegan)

#ploting
library(ggplot2)
library(scales)
library(gridExtra)

#### Define paths ####
# adjust paths to point to specific files
proj_folder <- #"E:/Project_2"
in_NOR_data <- #"E:/Project_2/HUI_CLM_layers/Norway_borders_N50"
in_veg_data <- #"E:/Project_2/HUI_CLM_layers/AR18x18_plots"
in_pts_data <- #"E:/Project_2/Single point plots CLM4.5"
in_RS_data <- #"E:/Project_2/LAND COVER DATA/Satveg_deling_nd_partnere_09_12_2009/tiff/"
in_binary_raster <- #"E:/Project_1_RUNS/new_MODEL_RUN/07_proportion_raster/Binary_raster/"
in_probab_raster <- #"E:/Project_1_RUNS/new_MODEL_RUN/04_predict_raster/"
in_raster <- #"E:/Project_1_FINAL_layers/MASKED_TIFF/"
in_dm_data <- #"E:/Project_2/OUTPUT/"
in_csv_files <- #"C:/Users/peterhor/Documents/GitHub/Project_2/"
out_folder <- #"E:/Project_2_RUNS/OUTPUT"
out_ggplots <- #"E:/Project_2_RUNS/OUTPUT/results/ggplots"
dir.create(file.path(out_folder),showWarnings=F)

#### CRS:  WGS84 UTM zone 33N ####
project_crs <- crs("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")



#### Create 1km boxes 20x####
veg_plots <- readOGR(in_veg_data, "Hui_dissolved_plots")
veg_plots@proj4string
#subset
# plot ID numbers to be subset (as in supplement S1)
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
#here, since DGVM is not available we can just crop with "Hui_dissolved_plots" and create centroids
veg_plots <- readOGR(in_veg_data, "Hui_dissolved_plots")
center_1081 <- gCentroid(veg_plots,byid=TRUE, id = veg_plots@data$FLATE_NR)

## frequency TEMP & PRECIP distribution diagrams ####
#20 vs 1081 plots along temp and precip gradient 
# load in rasters for temp and precip
bioclim_1 <- raster(paste0(in_raster,"bioclim_1.tif"))
bioclim_12 <- raster(paste0(in_raster,"bioclim_12.tif"))
bioclim_15 <- raster(paste0(in_raster,"bioclim_15.tif"))
# check plotting
plot(bioclim_15)
plot(center_20_test, add=TRUE)

c20_bioclim_1 <- extract(bioclim_1, center_20_test)
c1081_bioclim_1 <- extract(bioclim_1, center_1081)
c20_bioclim_12 <- extract(bioclim_12, center_20_test)
c1081_bioclim_12 <- extract(bioclim_12, center_1081)
c20_bioclim_15 <- extract(bioclim_15, center_20_test)
c1081_bioclim_15 <- extract(bioclim_15, center_1081)

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

var.test(c20_bioclim_15, c1081_bioclim_15)
t.test(c20_bioclim_15, c1081_bioclim_15)
ks.test(c20_bioclim_15, c1081_bioclim_15)

# publication plots in ggplot
temperature_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_bioclim_1,c1081_bioclim_1))
precipitation_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_bioclim_12,c1081_bioclim_12))
bio15_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_bioclim_15,c1081_bioclim_15))
swe_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_swe_10,c1081_swe_10))
tmin_data <- data.frame(type=rep(c("c_20", "c_1081"), c(20,1081)), data= c(c20_tmin_5,c1081_tmin_5))
library(ggplot2)
ggplot(temperature_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
ggplot(precipitation_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
ggplot(bio15_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
ggplot(swe_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
ggplot(tmin_data, aes(x=data, fill=type)) + geom_density(alpha=.3)
library(plyr)
t_mean <- ddply(temperature_data, "type", summarise, rating.mean=mean(data, na.rm=T))
t_mean
p_mean <- ddply(precipitation_data, "type", summarise, rating.mean=mean(data, na.rm=T))
p_mean
bio15_mean <- ddply(bio15_data, "type", summarise, rating.mean=mean(data, na.rm=T))
bio15_mean

t_plot <- ggplot(temperature_data, aes(x=data, fill=type)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_density(alpha=.3) + 
  geom_vline(data = t_mean, aes(xintercept = rating.mean, colour = type),
          linetype = "longdash", size=1) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 17), axis.text = element_text(size = 17)) +
  labs(title="frequency distribution", 
       #caption="Source: NCAR, SatVeg",
       x="Annual Mean Temperature (°C)",
       y="Density")

p_plot <- ggplot(precipitation_data, aes(x=data, fill=type)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_density(alpha=.3) + 
  geom_vline(data = p_mean, aes(xintercept = rating.mean, colour = type),
             linetype = "longdash", size=1) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 17), axis.text = element_text(size = 17)) +
  labs(title="frequency distribution", 
       #caption="Source: NCAR, SatVeg",
       x="Annual Precipitation (mm)",
       y="Density")

bio15_plot <- ggplot(bio15_data, aes(x=data, fill=type)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_density(alpha=.3) + 
  geom_vline(data = bio15_mean, aes(xintercept = rating.mean, colour = type),
             linetype = "longdash", size=1) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 17), axis.text = element_text(size = 17)) +
  labs(title="frequency distribution", 
       #caption="Source: NCAR, SatVeg",
       x="Precipitation Seasonality",
       y="Density")


#export
ggsave(plot = t_plot, filename = paste0("Temp_distribution_20vs1081",".png"), device = "png", path = paste0(out_ggplots),dpi = 300 )
ggsave(plot = p_plot, filename = paste0("Precip_distribution_20vs1081",".png"), device = "png", path = paste0(out_ggplots),dpi = 300 )
ggsave(plot = bio15_plot, filename = paste0("Precip_seasonality_20vs1081",".png"), device = "png", path = paste0(out_ggplots),dpi = 300 )

#### climate data test####
# cordex vs seNorge
# read in climate data from table S1 in the supplement, where centervalues for each plot are shown
cordex_check <-read.csv2("C:/Users/peterhor/Google Drive/2015 - University of Oslo/Project 2/20_plots_lat_long_alt_temp_precip_CORDEX.csv", sep = ";", dec = ".")
plot(cordex_check$SeNorge_precip~cordex_check$CORDEX_precip)
abline(a=1, b=1)
plot(cordex_check$SeNorge_temp~cordex_check$CORDEX_temp)
abline(a=1, b=1)

library(ggplot2)

climate_precip <- ggplot(cordex_check, aes(x=SeNorge_precip, y=CORDEX_precip)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_point(data = cordex_check, size=2) +
  geom_abline(intercept =  , slope = 1, color="black", 
              linetype="dashed", size=0.5) +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted", size=0.5, color="blue") +
  labs(title="CORDEX vs SeNorge", 
       #caption="Source: NCAR, SatVeg",
       x="SeNorge - Annual Precipitation (mm)",
       y="CORDEX - Annual Precipitation (mm)")
climate_precip
ggsave(plot = climate_precip, filename = paste0("CORDEX_vs_seNorge_precip.png"), device = "png", width = 10, height = 7, path = paste0(path),dpi = 300 )

climate_temp <- ggplot(cordex_check, aes(x=SeNorge_temp, y=CORDEX_temp)) + # scale_fill_discrete(name="Dataset", labels=c("1081 plots", "20 test plots")) +
  geom_point(data = cordex_check, size=2) +
  geom_abline(intercept =  , slope = 1, color="black", 
              linetype="dashed", size=0.5) +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted", size=0.5, color="red") +
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
# save output
write.csv2(AR_VT_20_perc, paste0(out_folder,"/AR18x18/AR_VT_20.csv"))

#### test representativness of VT coverage on the 20 plots ####
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


#### .... #####################
#### ** RS dataset ####
# Remote sensing rasters from SatVeg (https://kartkatalog.miljodirektoratet.no/MapService/Details/satveg)
RS_Norway33 <- raster(paste0(in_RS_data,"RS_Norway_utm33.tif"))
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

# reclassify VT to PFT ####
RS_Norway_PFT <- reclassify(RS_Norway33, transl_RS, include.lowest=TRUE, overwrite=TRUE,
                            filename = paste0(out_folder,"/RS/RS_PFT.tif")) 

RS_freq <- freq(RS_Norway33)
write.csv2(RS_freq, paste0(out_folder,"/RS/RS_freq.csv"))
RS_PFT_freq <- freq(RS_Norway_PFT)
write.csv2(RS_PFT_freq, file = paste0(out_folder,"RS/RS_PFT_freq.csv"))
#writeRaster(RS_Norway_PFT, filename = paste0(out_folder,"/RS/RS_Norway_PFT.tif"), format="GTiff", overwrite=TRUE)

# then mask out only 1x1km from raster (this takes 15min)
RS_PFT_single <- mask(RS_Norway_PFT, single_cell_1km_df, filename=paste0(out_folder,"/RS_PFT_single_cells_R"), format="GTiff", overwrite=TRUE)

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
# DM results from Horvath et al. 2019 (Dryad repository)

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

# all Norway
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
# reclassify VT to PFT ####
DM_PFT <- reclassify(max_val, transl_DM, include.lowest=TRUE, overwrite=TRUE,
                            filename = paste0(out_folder,"/DM/DM_PFT.tif")) 
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
#### TRANSLATION scheme ####
#### ******PFT*****#####
# prepare or load data 
# Translation schemes are based on the table S5 in supplement
#load in translation rules
transl_rule_DM <- read.csv(paste0(in_csv_files, "02_Translation_schemes/VEG_to_PFT_final.csv"), sep = ";", dec = ".")
transl_rule_RS <- read.csv(paste0(in_csv_files, "02_Translation_schemes/RS_to_PFT_final.csv"), sep = ";", dec = ".")
transl_rule_AR <- read.csv(paste0(in_csv_files, "02_Translation_schemes/AR18x18_to_PFT_final.csv"), sep = ";", dec = ".")



#### .... #####################     
## function RASTER to PFT table ####

# this function converts a raster (either DM or RS) on specified plots (20 selected plots)
# into a PFT profile - table with percentage representation of each PFT
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
  # for ease of saving 
  subfolder=20
  subtype="single_cells_"
  filetype="_sc"
  colnames(rast_mask_table) <- paste0("X", plots@data$id)
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
DM_maxval_PFT <- raster(paste0(out_folder,"/DM/DM_maxval_PFT.tif"))

# execute function 20squares
rast_to_table(RS_PFT, single_cell_1km_df)
rast_to_table(DM_maxval_PFT, single_cell_1km_df)

#### test representativness of PFT coverage on the 20 plots ####
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

#************####

pft_data <- PFT_df # PFT_df is created within 01_comparison_***.R script
path="E:/Project_2_RUNS/OUTPUT/results"
pft_data <- read.csv2(file = paste0(path, "/PFT_df.csv"))
colNames <- names(pft_data)[7:NCOL(pft_data)] 
#colNames <- names(pft_data)[7:NCOL(pft_data)] # could be also 7:12 in adjusted version of DGVM
pft_data <- read.csv2(file = paste0(path, "/PFT_df_adj.csv"), dec = ".") # adjusted values for precipitation seasonality (bioclim_15) from DGVM sensitivity analysis
colNames <- names(pft_data)[7:NCOL(pft_data)] 

######### GRAPHICS #############
# custom color scale
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


# * graphs FOR PUBLICATION ON A MAP IN QGIS ####
# no legend, no labels
colNames<-colnames(pft_data)[7:ncol(pft_data)] #[5:24] #
colNames_number<-gsub("X", "#", colNames)

#Create a custom color scale (https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin)

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
     theme_void() + theme(legend.position="none")      # hide legend and labels on X axis  # theme_void()
  ggsave(plot = g1, filename = paste0("AR_vs_Methods_DM_maxval",i,".png"), device = "png", width = 10, height = 7, path = paste0(out_ggplots, "/plots_3way/"),dpi = 300 ) #DGVM_adj
  #("AR_vs_Methods_DM_maxval",i,".png")
}

### ....####
### ....####

#### Improving the representation of high-latitude vegetation in Dynamic Global Vegetation Models ####
# Dynamical Vegetation model version CLM-DV4.5
# updated: 20.10.2020
# author: Peter Horvath

#### COMPARISON Scheme####
# of different modeling approaches 

library(vegan)
library(ggplot2)
# read in data
in_csv_files <- "C:/Users/peterhor/Documents/GitHub/Project_2/01_Comparison_schemes"
#path="C:/Users/peterhor/Documents/GitHub/Project_2/output"
path="E:/Project_2_RUNS/OUTPUT/results"
csv_list <- list.files(in_csv_files, pattern="csv", full.names=FALSE)
# read in with column and row names

#read in from GIThub original data 1x1km
in_csv_files <- "C:/Users/peterhor/Documents/GitHub/GitHub_Enterprise/Project_2/01_Comparison_schemes/"
comp_DM_maxval <- read.csv2(paste0(in_csv_files,"DM_maxval_PFT.csv"), sep = ";", dec = ".")
comp_RS <- read.csv2(paste0(in_csv_files,"RS_PFT.csv"), sep = ";", dec = ".")
comp_DGVM <- read.csv2(paste0(in_csv_files,"DGVM_PFT.csv"), sep = ";", dec = ".")
comp_AR18x18 <- read.csv2(paste0(in_csv_files,"AR_PFT.csv"), sep = ";", dec = ".")

########***********##########
# testing 20 plots####


PFT_method_names <-list(comp_DM_maxval, comp_RS, comp_DGVM, comp_AR18x18) #comp_DM_maxval, comp_DM_preval,
names(PFT_method_names) <- c("DM_maxval", "RS", "DGVM", "AR18x18")
# only data with no labels
PFT_method_raw <- lapply(PFT_method_names, "[", 5:24)
#PFT_method_raw[[1]][,2]
names(PFT_method_raw) <- c("DM_maxval", "RS", "DGVM", "AR18x18")
saveRDS(PFT_method_raw, file = paste0(path, "/PFT_profiles20.RDS"))
PFT_method_raw <- readRDS(file = paste0(path, "/PFT_profiles20.RDS"))
#str(pft_data)

# 0 ####
# before we merge and compare all the methods, let's explore each method and the sites
# explore list
lapply(PFT_method_raw, colSums)
lapply(PFT_method_raw, rowSums)
# rename rows in list of DFs
PFT_method_raw<-lapply(PFT_method_raw, function(x){
  rownames(x) <- t(comp_DM_maxval[,2])
  return(x)
})
norway_PFT <- as.data.frame(lapply(PFT_method_raw, rowSums))
#norway_PFT <- format(norway_PFT, digits=1)
colSums(norway_PFT)
saveRDS(norway_PFT, file = paste0(path, "/Norway_PFTs_ar20.RDS"))
write.csv2(norway_PFT, file = paste0(path, "/Norway_PFT_ar20.csv"))

norway_PFT_perc <- norway_PFT/20

cumsum(norway_PFT_perc)
colsum(norway_PFT_perc)
colSums(norway_PFT_perc)
saveRDS(norway_PFT_perc, file = paste0(path, "/Norway_PFTs_ar20_perc.RDS"))
write.csv2(norway_PFT_perc, file = paste0(path, "/Norway_PFT_ar20_perc.csv"))
write.csv2(norway_PFT_perc[,c(5,4,2,6)], file = paste0(path, "/Norway_PFT_ar20_perc_TABLE3results.csv"))
norway_PFT_perc <- readRDS(file = paste0(path, "/Norway_PFTs_ar20_perc.RDS"))
# area statistics NORWAY testing method against reference AR 
# CHI
chisq.test(cbind(norway_PFT_perc$DM_maxval, norway_PFT_perc$AR18x18))
chisq.test(cbind(norway_PFT_perc$RS, norway_PFT_perc$AR18x18))
chisq.test(cbind(norway_PFT_perc$DGVM, norway_PFT_perc$AR18x18))
# FISHER
fisher.test(cbind(norway_PFT_perc$DM_maxval, norway_PFT_perc$AR18x18))
fisher.test(cbind(norway_PFT_perc$RS, norway_PFT_perc$AR18x18))
fisher.test(cbind(norway_PFT_perc$DGVM, norway_PFT_perc$AR18x18))
#KOLMOGOROV-Smirnov test
ks.test(norway_PFT_perc$DM_maxval, norway_PFT_perc$AR18x18)

#testing vegdist across 20 sites merged
vegdist(t(norway_PFT))
vegdist(t(PFT_method_raw[[1]]))
hist(vegdist(t(PFT_method_raw[[1]])), main = names(PFT_method_raw)[[1]])
lapply(PFT_method_raw, function(x){vegdist(t(x))})
lapply(PFT_method_raw, function(x){hist(vegdist(t(x)))})


# 
bray_sites <- lapply(PFT_method_raw, function(x){
  res <- vegdist(t(x))
  return(round(res,2))
})
lapply(bray_sites,summary)
# plot histogram of similarity between sites
lapply(bray_sites, function(x) hist(x))
# library(ggplot2)
# ggplot(bray_sites[[1]], aes(x=DM_maxval) + geom_histogram())
# lapply(bray_sites, function(g) ggplot(g, aes=(x=)))

#PFT_method_raw[["DGVM"]][,2]
lapply(PFT_method_raw, "[", 4) #list of dataframes on the fourth place of parent list
lapply(PFT_method_raw, "[[", 4) #list of vectors on the fourth place of parent list

# we can compute Bray-Curtis dissimilarities for each site individually 

#In bray-curtis dissimilarity calculations, the rows should represent sites (in our case methods),
#while the colums are species (PFTs)
# thus need to transpose data

#1 Bray-Curtis dissimilarities for each site ####
# packed into a function
bray_curtis <- function(site){
  site_1_20<- t(as.data.frame(lapply(PFT_method_raw, "[", site)))
  rownames(site_1_20) <- c( "DM_maxval", "RS", "DGVM", "AR18x18")
  colnames(site_1_20) <- t(comp_DM_maxval[,2])
  bray_site_1_20<-vegdist(site_1_20, method = "bray")
  return(bray_site_1_20)
}
# note site locations are on columns 4:23
#ncol(PFT_method_raw[[1]])
bray_curtis(5)  

a <- vector('list', length(2))
for(i in 1:20){
  a[[i]] <- bray_curtis(i)
  names(a)[i] <- colnames(PFT_method_raw[[1]])[i]
}
# preparing the output
a
a[[4]]
a[[20]][c(5,9,12,14,15)]
mypos <- c(5,9,12,14,15) #positions for BC dissimilarity of AR18x18 to each of the other methods
#extracting positions from each list
sapply(a, "[", mypos)
lapply(a, "[", mypos)
#getting means for each Method compared to AR18x18
mean(sapply(a, "[", mypos[1]))
mean(sapply(a, "[", mypos[2]))
mean(sapply(a, "[", mypos[3]))

#plotting mean BC Dissimilarity for method across sites
# alternatively load Braycurtis_3methods.csv from file #
bc_dis <- read.csv2(paste0(path, "/braycurtis_3methods.csv"), sep = ";", dec=".", header = TRUE)
bc_dis <- as.data.frame(t(sapply(a, "[", mypos)))
bc_dis <- cbind(rownames(bc_dis),bc_dis)
colnames(bc_dis) <- c("plot_nr", "DM_maxval", "RS", "DGVM")
summary(bc_dis)
write.csv2(bc_dis, file = paste0(path, "/braycurtis_ar20.csv"))
write.csv2(bc_dis[,c("DM_maxval", "RS", "DGVM")], file = paste0(path, "/braycurtis_3methods_ar20.csv"))
#rowMeans(bc_dis)
#colMeans(bc_dis)
#
bc_dis <- read.csv2(file = paste0(path, "/braycurtis_3methods_ar20.csv"))
wilcox.test(bc_dis$DM_maxval, bc_dis$DGVM)
wilcox.test(bc_dis$RS, bc_dis$DGVM)
wilcox.test(bc_dis$RS, bc_dis$DM_maxval)

# plot BC####
library(reshape2)
bc_dis_l <- melt(bc_dis)
colnames(bc_dis_l) <- c("plot_nr","Method", "Bray_Curtis")
saveRDS(bc_dis_l, file = paste0(path, "/braycurtis_long.RDS"))
#bc_dis_l <- readRDS(file = paste0(path, "/braycurtis_long.RDS"))
write.csv2(bc_dis_l, file = paste0(path, "/braycurtis_long.csv"))
## all methods including "DM_maxval", "DM_preval"
bc_plot1 <- ggplot(bc_dis_l, aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Modelling method",
       y="Proportional dissimilarity") 
bc_plot1
ggsave(filename = "Project_2_bc_plot1.png",plot = bc_plot1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

#only three relevant methods
bc_plot3 <- ggplot(bc_dis_l[bc_dis_l$Method %in% c("DM_maxval", "RS", "DGVM"),], aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Modelling method",
       y="Proportional dissimilarity") 
bc_plot3
ggsave(filename = "Project_2_bc_methods3.png",plot = bc_plot3, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

# with corresponding plot numbers connected with grey line
bc_plot3.1 <- ggplot(bc_dis_l[bc_dis_l$Method %in% c("DM_maxval", "RS", "DGVM"),], aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Modelling method",
       y="Proportional dissimilarity") +
  geom_line(aes(group=plot_nr), colour="grey",  alpha=0.5)
bc_plot3.1
ggsave(filename = "Project_2_bc_methods3.1.png",plot = bc_plot3.1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

bc_plot6.1 <- ggplot(bc_dis_l, aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Modelling method",
       y="Proportional dissimilarity") +
  geom_line(aes(group=plot_nr), colour="grey",  alpha=0.5)
bc_plot6.1
ggsave(filename = "Project_2_bc_methods6.1.png",plot = bc_plot6.1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

# DM_maxval with corresponding plot numbers connected with grey line
# modorder <- c("DGVM","RS","DM")
bc_plot3.11 <- ggplot(bc_dis_l[bc_dis_l$Method %in% c("DGVM","RS","DM_maxval"),], aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  #scale_color_manual(labels=c("DGVM","RS","DM")) +
  scale_x_discrete(limits= c("DGVM","RS","DM_maxval"),labels=c("DGVM","RS","DM")) + # changing the order and the label names
  #scale_fill_discrete(name = "Method", guide=c("DGVM","RS","DM_maxval"),  labels=c("DGVM","RS","DM")) + # change legend # eventually + 
  theme(legend.position = "none") +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Modelling method",
       y="Proportional dissimilarity") +
  geom_line(aes(group=plot_nr), colour="grey",  alpha=0.5)
bc_plot3.11
ggsave(filename = "Project_2_bc_methods3.11_ar20_new2.png",plot = bc_plot3.11, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )



# 1.1. extension on spatial patterns ####
#categorize plots into South-North, and Coast-Inland.
#assess influence of latitude (select north of 62degrees)
#assess influence of coastalness ()
library(rgdal)
center_20plots <- readOGR(in_pts_data, "centerpoint_21test")
#the lat long is swapped in the dataset
colnames(center_20plots@data)[5] <- "LONG"
colnames(center_20plots@data)[6] <- "LAT"
proxy_coast <- raster("D:/Project_1_FINAL_layers/MASKED_TIFF_uncorrelated/CONTINUOUS/Proxy_coast.tif")
# values to points
proxy_vals <- extract(proxy_coast, center_20plots)
# add data to SpatialPointsDataFrame
center_20plots@data$proxy_coast <- proxy_vals

sp_pattern <- cbind(center_20plots@data, bc_dis)
sp_pattern<-sp_pattern[,-c(1:4)]
# 1.12 statistical tests####
summary(lm(LAT~RS, data = sp_pattern))
summary(lm(LAT~DM_maxval, data = sp_pattern))
summary(lm(LAT~DGVM, data = sp_pattern))
summary(lm(proxy_coast~RS, data = sp_pattern))
summary(lm(proxy_coast~DM_maxval, data = sp_pattern))
summary(lm(proxy_coast~DGVM, data = sp_pattern))

# 1.13 plot tests####

library(reshape2)
sp_pattern_l <- melt(sp_pattern, id = c("LAT", "LONG", "plot_nr", "proxy_coast")) #melt with ID also on LAT & Long
colnames(sp_pattern_l) <- c("LAT","LONG","plot_nr","proxy_coast","Method", "Bray_Curtis")

# look at latitude influence 
sp_plot1 <- ggplot(sp_pattern_l, aes(x=Bray_Curtis, y=LONG, colour=Method, fill=Method)) + 
  geom_point(alpha=0.8, size = 3, shape = 21) + 
  facet_wrap(~Method, scales = "free") +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted") + #fill="grey",linetype="dotted",
  labs(x="Proportional dissimilarity index across methods",
       y="Latitude") 
sp_plot1
ggsave(filename = "Project_2_sp_plot1.png",plot = sp_plot1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )
# look at coast proximity influence
sp_plot2 <- ggplot(sp_pattern_l, aes(x=Bray_Curtis, y=proxy_coast,colour=Method, fill=Method)) + 
  geom_point(alpha=0.8, size = 3, shape = 21) +
  facet_wrap(~Method, scales = "free") +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted") + #fill="grey",
  labs(x="Proportional dissimilarity index along methods", 
       y="proximity to coast") 
sp_plot2
ggsave(filename = "Project_2_sp_plot2.png",plot = sp_plot2, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

# 3final methods look at latitude influence 
sp_plot1 <- ggplot(sp_pattern_l[sp_pattern_l$Method %in% c("DM_maxval", "RS", "DGVM"),], aes(x=Bray_Curtis, y=LONG, colour=Method, fill=Method)) + 
  geom_point(alpha=0.8, size = 3, shape = 21) + 
  facet_wrap(~Method, scales = "free") +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted") + #fill="grey",linetype="dotted",
  labs(x="Proportional dissimilarity index across methods",
       y="Latitude") 
sp_plot1
ggsave(filename = "Project_2_sp_plot1_3methods.png",plot = sp_plot1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )
# 3final methods look at coast proximity influence
sp_plot2 <- ggplot(sp_pattern_l[sp_pattern_l$Method %in% c("DM_maxval", "RS", "DGVM"),], aes(x=Bray_Curtis, y=proxy_coast,colour=Method, fill=Method)) + 
  geom_point(alpha=0.8, size = 3, shape = 21) +
  facet_wrap(~Method, scales = "free") +
  stat_smooth(formula=y~x, method = lm, alpha=0.2, linetype="dotted") + #fill="grey",
  labs(x="Proportional dissimilarity index along methods", 
       y="proximity to coast") 
sp_plot2
ggsave(filename = "Project_2_sp_plot2_3methods.png",plot = sp_plot2, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )


#### 2 working with data frame#### 
#MELT DATA FROM LIST INTO LONG FORM FOR PLOTTING IN GGPLOT
library(reshape2)
library(data.table)
PFT_df <- rbindlist(PFT_method_names, idcol = TRUE)
colnames(PFT_df)[1] <- "Method"
colnames(PFT_df)[4] <- "PFT_code"
saveRDS(PFT_df, file = paste0(path, "/PFT_df.RDS"))
write.csv2(PFT_df, file = paste0(path, "/PFT_df.csv"))
PFT_melt <- melt(PFT_df[,-3], variable.name = "site", value.name = "percentage") # exclude 4th column as it includes a integer name for PFTs - and is restraining from proper use of melt function
saveRDS(PFT_melt, file = paste0(path, "/PFT_melt.RDS"))
write.csv2(PFT_melt, file = paste0(path, "/PFT_melt.csv"))
PFT_melt<- read.csv2(file = paste0(path, "/PFT_melt.csv"))

#Pft representation per site with regard to method
pft <- ggplot(PFT_melt, aes(x=PFT_code, y=percentage, fill=Method)) +
  geom_bar(position = "stack",stat = "identity") + 
  facet_wrap(~site, scales = "free")+
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + 
  labs(title="TRALALA", 
       fill="FILL",
       x="PFT",
       y=expression ("percent i"~m^2))
pft
ggsave(filename = "PFT_nonsense.png",plot = pft, width = 15, height = 12, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )


#Pft representation with regard to method 
# we can observe that there is a general lack of representation of PFT 8 in DGVM and a total overrepresentation of PFT2
pft1 <- ggplot(PFT_melt[PFT_melt$Method %in% c("AR18x18","DM_maxval", "RS", "DGVM"),], aes(x=PFT_shortcut, y=percentage, fill = PFT_name)) +
  geom_bar(position = "stack",stat = "identity") + 
  facet_wrap(~Method)+ #, scales = "free", ncol = 1
  colScale + #scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + 
  labs(title="Coverage of PFTs per Method", 
       fill="Legend",
       x="PFT",
       y=expression ("cumulative percent "))
pft1
ggsave(filename = "PFT_method_3.1_DM_maxval.png",plot = pft1, width = 15, height = 12, device = "png", path = paste0(path, "/ggplots/"),dpi = 300 )
# position dodge
pft1.1 <- ggplot(PFT_melt[PFT_melt$Method %in% c("AR18x18","DM_maxval", "RS", "DGVM"),], aes(x=PFT_shortcut, y=percentage, fill=Method)) +
  geom_bar(position = "dodge",stat = "identity") + 
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + 
  labs(title="Coverage of PFTs per Method", 
       fill="Legend",
       x="PFT",
       y=expression ("cumulative percent "))
pft1.1
ggsave(filename = "PFT_method_3.1_DM_maxval.png",plot = pft1.1, width = 15, height = 12, device = "png", path = paste0(path, "/ggplots/"),dpi = 300 )

pft1.2 <- ggplot(PFT_melt, aes(x=PFT_code, y=percentage, fill = PFT_name)) +
  geom_bar(position = "stack",stat = "identity") + 
  facet_wrap(~Method)+ #, scales = "free", ncol = 1
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + 
  labs(title="PFT per Method", 
       fill="FILL",
       x="PFT",
       y=expression ("cumulative percent "))
pft1.2
ggsave(filename = "PFT_method_6.1.png",plot = pft1.2, width = 15, height = 12, device = "png", path = paste0(path, "/ggplots/"),dpi = 300 )



