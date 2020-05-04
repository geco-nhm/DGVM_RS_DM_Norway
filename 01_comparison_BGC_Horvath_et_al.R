#### Improving the representation of high-latitude vegetation in Dynamic Global Vegetation Models ####
# Dynamical Vegetation model version CLM-DV4.5
# updated: 30.4.2020
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
comp_DM_preval <- read.csv2(paste0(in_csv_files,"DM_preval_PFT.csv"), sep = ";", dec = ".")
comp_DM_AUC <- read.csv2(paste0(in_csv_files,"DM_AUC_PFT.csv"), sep = ";", dec = ".")
comp_RS <- read.csv2(paste0(in_csv_files,"RS_PFT.csv"), sep = ";", dec = ".")
comp_DGVM <- read.csv2(paste0(in_csv_files,"DGVM_PFT.csv"), sep = ";", dec = ".")
comp_AR18x18 <- read.csv2(paste0(in_csv_files,"AR_PFT.csv"), sep = ";", dec = ".")
#PFT_all <- read.csv(paste0(in_csv_files,"PFT_descr.csv"), sep = ";", dec = ".")

#adjust rows with descriptions
PFT_descr <- read.csv("C:/Users/peterhor/Documents/GitHub/GitHub_Enterprise/Project_2/02_Translation_schemes/PFT_descr.csv", sep = ";", dec = ".")
comp_DM_maxval<- cbind(PFT_descr[1:6,],comp_DM_maxval[,2:ncol(comp_DM_maxval)])
comp_DM_preval<- cbind(PFT_descr[1:6,],comp_DM_preval[,2:ncol(comp_DM_preval)])
comp_DM_AUC<- cbind(PFT_descr[1:6,],comp_DM_AUC[,2:ncol(comp_DM_AUC)])
comp_RS<- cbind(PFT_descr[1:6,],comp_RS[,2:ncol(comp_RS)])
comp_DGVM<- cbind(PFT_descr[1:6,],comp_DGVM[,4:ncol(comp_DGVM)])
comp_AR18x18<- cbind(PFT_descr[1:6,],comp_AR18x18[,4:ncol(comp_AR18x18)])


########***********##########
# testing 20 plots####


PFT_method_names <-list(comp_DM_AUC, comp_DM_maxval, comp_DM_preval, comp_RS, comp_DGVM, comp_AR18x18) #comp_DM_maxval, comp_DM_preval,
names(PFT_method_names) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
# only data with no labels
PFT_method_raw <- lapply(PFT_method_names, "[", 5:24)
#PFT_method_raw[[1]][,2]
names(PFT_method_raw) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
saveRDS(PFT_method_raw, file = paste0(path, "/PFT_profiles20x6.RDS"))
PFT_method_raw <- readRDS(file = paste0(path, "/PFT_profiles20x6.RDS"))
#str(pft_data)

# 0 ####
# before we merge and compare all the methods, let's explore each method and the sites
# explore list
lapply(PFT_method_raw, colSums)
lapply(PFT_method_raw, rowSums)
# rename rows in list of DFs
PFT_method_raw<-lapply(PFT_method_raw, function(x){
  rownames(x) <- t(comp_DM_AUC[,2])
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
# ggplot(bray_sites[[1]], aes(x=DM_AUC) + geom_histogram())
# lapply(bray_sites, function(g) ggplot(g, aes=(x=)))

#PFT_method_raw[["DGVM"]][,2]
lapply(PFT_method_raw, "[", 4) #list of dataframes on the fourth place of parent list
lapply(PFT_method_raw, "[[", 4) #list of vectors on the fourth place of parent list

# we can compute Bray-Curtis dissimilarities for each site individually 

#In bray-curtis dissimilarity calculations, the rows should represent sites (in our case methods),
#while the colums are species (PFTs)
# thus need to transpose data

X405<- t(as.data.frame(lapply(PFT_method_raw, "[", 4)))
rownames(X405) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
colnames(X405) <- t(comp_DM_AUC[,2])
bray_X405<-vegdist(X405)
bray_X405


#0.1 chi square test for each combination of the methods ####
chi_trial <- t(as.data.frame(lapply(PFT_method_raw, "[", 5)))
rownames(chi_trial) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
colnames(chi_trial) <- t(comp_DM_AUC[,2])
fisher.test(chi_trial[c(1,6),c(1,5,6)])
lapply(PFT_method_raw, "[", c(4,6))


site=4
fish <- function(site){
  site_1_20<- t(as.data.frame(lapply(PFT_method_raw, "[", site)))
  rownames(site_1_20) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
  colnames(site_1_20) <- t(comp_DM_AUC[,2])
  f_site_1_20<-fisher.test(site_1_20[c(2,6),])
  return(f_site_1_20)
}

PFT_methods_only = PFT_method_raw[c(2,4,5)]
PFT_reference = PFT_method_raw[[6]]

for(method in PFT_methods_only){
  for (name in names(method)){
    site = method[name]
    reference = PFT_reference[name]
    print(site)
    print(reference)
  }
  
}

typeof(unlist(PFT_method_raw[[5]]))
a = PFT_reference[1]
names(a)[[1]]
PFT_method_raw[[2]][1]

a = cbind(PFT_method_raw[[4]][1],PFT_method_raw[[6]][1])
fisher.test(a, workspace = 200000000)
chisq.test(a)
t.test(PFT_method_raw[[5]][[1]][3],v, paired = TRUE)
PFT_method_raw[[5]]


#1 Bray-Curtis dissimilarities for each site ####
# packed into a function
bray_curtis <- function(site){
  site_1_20<- t(as.data.frame(lapply(PFT_method_raw, "[", site)))
  rownames(site_1_20) <- c("DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM", "AR18x18")
  colnames(site_1_20) <- t(comp_DM_AUC[,2])
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
colnames(bc_dis) <- c("plot_nr","DM_AUC", "DM_maxval", "DM_preval", "RS", "DGVM")
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
  labs(x="Vegetation Modelling Method",
       y="Bray Curtis Dissimilarity index") 
bc_plot1
ggsave(filename = "Project_2_bc_plot1.png",plot = bc_plot1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

#only three relevant methods
bc_plot3 <- ggplot(bc_dis_l[bc_dis_l$Method %in% c("DM_AUC", "RS", "DGVM"),], aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Vegetation Modelling Method",
       y="Bray Curtis Dissimilarity index") 
bc_plot3
ggsave(filename = "Project_2_bc_methods3.png",plot = bc_plot3, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

# with corresponding plot numbers connected with grey line
bc_plot3.1 <- ggplot(bc_dis_l[bc_dis_l$Method %in% c("DM_AUC", "RS", "DGVM"),], aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Vegetation Modelling Method",
       y="Bray Curtis Dissimilarity index") +
  geom_line(aes(group=plot_nr), colour="grey",  alpha=0.5)
bc_plot3.1
ggsave(filename = "Project_2_bc_methods3.1.png",plot = bc_plot3.1, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )

bc_plot6.1 <- ggplot(bc_dis_l, aes(x=Method, y=Bray_Curtis, fill=Method)) + 
  geom_boxplot(alpha=0.8) +
  geom_point(aes(fill = Method), alpha=0.8, size = 3, shape = 21, position = position_jitterdodge()) +
  labs(x="Vegetation Modelling Method",
       y="Bray Curtis Dissimilarity index") +
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
  labs(x="Vegetation Modelling Method",
       y="Bray Curtis Dissimilarity index") +
  geom_line(aes(group=plot_nr), colour="grey",  alpha=0.5)
bc_plot3.11
ggsave(filename = "Project_2_bc_methods3.11_ar20_new.png",plot = bc_plot3.11, device = "png", path = paste0(path, "/ggplots"),dpi = 300 )



#

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

#Pft representation with regard to method + along individual SITES
# 
pft2 <- ggplot(subset(PFT_melt, site=="X405"), aes(x=PFT_code, y=percentage, fill = PFT_name)) 
pft2 <- pft2 + geom_bar(position = "stack",stat = "identity") + 
  facet_wrap(~Method, scales = "free")+
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + 
  labs(title="PFT per Method at site X405", 
       fill="FILL",
       x="PFT",
       y=expression ("cumulative percent "))
pft2
ggsave(filename = "PFT_method_X405.png",plot = pft2, width = 15, height = 12, device = "png", path = paste0(path, "/ggplots/"),dpi = 300 )


