#### Set up ####
#Clear environment
rm(list=ls())

#Set working directory
setwd("/Users/ab61/Documents/Papers/Cholera/Nature Submission/Revisions Round 2/Reproducibility")

#Activate packages
library(readxl)
library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(ggstream)
library(RColorBrewer)
library(tidytree)
library(treeio)
library(bangladesh)
library(plyr)
library(sf)
library(gridExtra)
library(phytools)
library(CoordinateCleaner)
library(rnaturalearth)
library(ggrepel)
library(viridis)
library(ggridges)
library(reshape2)
library(scatterpie)
library(gggenes)
library(ggnewscale)
library(colorspace)

#Vibrio cholerae metadata
metadata_all = as.data.frame(read.delim("SupplementaryDocument1.txt"), header = T, fill = T)
metadata_all$ID = gsub("#", "_", metadata_all$ID)

#Create dates
metadata_all$Date_Exact = paste(metadata_all$Year, sprintf("%02d", metadata_all$Month), sprintf("%02d", metadata_all$Day), sep = "-")
metadata_all[grep("NA", metadata_all$Date_Exact),"Date_Exact"] = NA
metadata_all$Date_Exact = as.Date(metadata_all$Date_Exact)

metadata_all$Date_Month = paste(metadata_all$Year, sprintf("%02d", metadata_all$Month),"01", sep = "-")
metadata_all[grep("NA", metadata_all$Date_Month ),"Date_Month"] = NA
metadata_all$Date_Month = as.Date(metadata_all$Date_Month )

#Subset to those with over 90% reads attributed to V. cholerae
metadata = metadata_all[which(metadata_all$Excluded == "No"),]

#SNP distances
snpdists = as.matrix(read.delim("snpdists.txt", row.names = 1))
colnames(snpdists) = rownames(snpdists)
colnames(snpdists) = gsub("#", "_", colnames(snpdists))
rownames(snpdists) = gsub("#", "_", rownames(snpdists))

#Tree of all 7PET Vibrio Cholerae
vc = read.tree("Trees/vibriocholerae.tree")
vc_annotations = read.table("Trees/vibriocholerae.annotations", header = T)
vc = list(tree = vc, annotation = vc_annotations)
rm(vc_annotations)

#Tree of all 7PET Vibrio Cholerae, with outliers removed so can plot a timed phylogeny
vc_clean = read.tree("Trees/vibriocholerae_clean.tree")
vc_clean_annotations = read.table("Trees/vibriocholerae_clean.annotations", header = T)
vc_clean = list(tree = vc_clean, annotation = vc_clean_annotations)
rm(vc_clean_annotations)

#Tree of BD1 in surveillance study only
bd1 = read.tree("Trees/bd1.tree")
bd1_annotations = read.table("Trees/bd1.annotations", header = T)
bd1 = list(tree = bd1, annotation = bd1_annotations)
rm(bd1_annotations)

#Tree of BD2 in surveillance study only
bd2 = read.tree("Trees/bd2.tree")
bd2_annotations = read.table("Trees/bd2.annotations", header = T)
bd2 = list(tree = bd2, annotation = bd2_annotations)
rm(bd2_annotations)

#Function to remove missing data
no_na = function(x){x[which(!is.na(x))]}

#### Define sublineages ####
for(lineage in c("sBD1", "BD2", "Other")){
  ids = metadata_all[which(metadata_all$Lineage == lineage),"ID"]
  lineage_dist_matrix = snpdists[which(rownames(snpdists) %in% ids), which(colnames(snpdists) %in% ids)]
  hc = hclust(as.dist(lineage_dist_matrix))
  Sublineages = cutree(hc, h = 20)
  Sublineages = data.frame(ID = names(Sublineages), Sublineages = Sublineages)
  Sublineages = merge(Sublineages, vc$annotation[,c("label", "y")], by.x = "ID", by.y = "label")
  #Identify singletons
  singletons = names(which(table(Sublineages$Sublineages) == 1))
  singletons = Sublineages[which(Sublineages$Sublineages %in% singletons),"ID"]
  Sublineages = Sublineages[-which(Sublineages$ID %in% singletons),]
  #Name sublineages by position in the tree. To order by when the sublineage first appears in the tree but also exclude outliers, look for 5% percentile
  percentile_5 = function(x){quantile(x, probs = 0.05)}
  order_in_tree = aggregate(Sublineages$y, by = list(Sublineages = Sublineages$Sublineages), FUN = "percentile_5")
  order_in_tree = order_in_tree[order(order_in_tree$x),]
  order_in_tree$SublineageName = paste(lineage, sprintf("%03d", as.numeric(factor(order_in_tree$Sublineages, levels = unique(order_in_tree$Sublineages)))), sep = ".")
  Sublineages = merge(Sublineages, order_in_tree, by = "Sublineages")
  #Create a dataframe with sublineages for all lineages
  if(lineage == "sBD1"){SublineageNames = Sublineages[,c("ID", "SublineageName")]}else{SublineageNames = rbind(SublineageNames, Sublineages[,c("ID", "SublineageName")])}
}
colnames(SublineageNames)[2] = "Sublineage"

metadata = merge(metadata, SublineageNames, by = "ID", all.x = T)
metadata_all = merge(metadata_all, SublineageNames, by = "ID", all.x = T)
#Label Singletons as "Singleton"
metadata[which(is.na(metadata$Sublineage)),"Sublineage"] = paste(metadata[which(is.na(metadata$Sublineage)),"Lineage"], "Singleton", sep = ".")

save(metadata, file = "metadata_sublineages.Rdata")


#### Sublineage dynamics (Fig. 1a) ####
load(file = "metadata_sublineages.Rdata")

metadata$DecimalDate = decimal_date(metadata$Date_Month)
dynamics = as.data.frame(table(metadata[which(metadata$Country == "Bangladesh" & metadata$Reference == "This manuscript"),c("DecimalDate", "Sublineage")]))
dynamics$DecimalDate = as.numeric(as.character(dynamics$DecimalDate))

#Create colour scheme for sublineages
colorscheme = data.frame(Sublineage = sort(unique(dynamics$Sublineage)), color = NA)
colorscheme[grep("BD2", colorscheme$Sublineage),"color"] = (colorRampPalette(brewer.pal(9, "GnBu"))(length(colorscheme[grep("BD2", colorscheme$Sublineage),"color"])))
colorscheme[grep("sBD1", colorscheme$Sublineage),"color"] = (colorRampPalette(brewer.pal(9, "RdPu"))(length(colorscheme[grep("sBD1", colorscheme$Sublineage),"color"])))

ggplot(dynamics, aes(x = DecimalDate, y = Freq, fill = Sublineage)) + 
  geom_stream(type = "ridge", color = "black", bw = 0.4) +
  scale_fill_manual(values = colorscheme[which(colorscheme$Sublineage %in% dynamics$Sublineage),"color"]) +  
  theme_bw() + 
  theme(legend.position = "right", panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10)) +
  annotate(geom = "rect", xmin = 2016, xmax = 2016.333, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") +
  annotate(geom = "rect", xmin = 2018.5, xmax = 2022, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") +
  labs(x = "Year", y = "Frequency", fill = "Sublineage") +
  annotate(geom = "text", x = 2015, y = 50, label = "Nationwide Surveillance\n10 sites") +
  annotate(geom = "text", x = 2017.3, y = 120, label = "Nationwide Surveillance\n22 sites") +
  annotate(geom = "text", x = 2022.6, y = 30, label = "2% study\nDhaka") 

ggsave("1a Sublineage Dynamics.svg", width = 10, height = 4)

#### MGE/Gene dynamics (Fig. 1b and Extended data Fig. 3) ####
dynamics = as.data.frame(table(metadata[which(metadata$Country == "Bangladesh" & metadata$Reference == "This manuscript"),c("DecimalDate", "Sublineage")]))
dynamics$DecimalDate = as.numeric(as.character(dynamics$DecimalDate))

#Categorise neatly by MGE/Gene Profile
study = metadata[which(metadata$Country == "Bangladesh" & metadata$Reference == "This manuscript"),]
study$MGE_Replacement = ""
#VSP-II
study[grep("K147fs", study$VC_0490.var),"MGE_Replacement"] = "VC0490-K147fs"
study[grep("VC_0491-ins", study$VSP_II),"MGE_Replacement"] = paste(study[grep("VC_0491-ins", study$VSP_II),"MGE_Replacement"], "VC0491-ins")
study[grep("VC_0491-VC_0498", study$VSP_II),"MGE_Replacement"] = paste(study[grep("VC_0491-VC_0498", study$VSP_II),"MGE_Replacement"], "ΔVC0491-VC0494")
#Ogawa/Inaba
study[which(study$VC_0258.assembled != "yes"),"MGE_Replacement"] = paste(study[which(study$VC_0258.assembled != "yes"),"MGE_Replacement"], "wbeT-ins")
#PLE
study[which(study$PLE == "PLE1"),"MGE_Replacement"] = paste(study[which(study$PLE == "PLE1"),"MGE_Replacement"] , "PLE1")
study[which(study$PLE == "PLE11"),"MGE_Replacement"] = paste(study[which(study$PLE == "PLE11"),"MGE_Replacement"] , "PLE11")
#SXT
study[which(study$SXT == "ICE-GEN"),"MGE_Replacement"] = paste(study[which(study$SXT == "ICE-GEN"),"MGE_Replacement"] , "ICE-GEN")
study[which(study$SXT == "ICE-TET"),"MGE_Replacement"] = paste(study[which(study$SXT == "ICE-TET"),"MGE_Replacement"] , "ICE-TET")
#Summarise
study$MGE_Replacement = trimws(study$MGE_Replacement)
levels = c("VC0490-K147fs PLE1 ICE-TET", "VC0490-K147fs VC0491-ins PLE1 ICE-TET", "VC0490-K147fs ΔVC0491-VC0494 PLE1 ICE-TET","VC0490-K147fs ΔVC0491-VC0494 wbeT-ins PLE1 ICE-TET", "VC0490-K147fs ΔVC0491-VC0494 wbeT-ins ICE-TET", "VC0490-K147fs ΔVC0491-VC0494 wbeT-ins", "ICE-GEN", "PLE11 ICE-GEN", "wbeT-ins PLE11 ICE-GEN")
study[which(!study$MGE_Replacement %in% levels),"MGE_Replacement"] = "Other"


dyn = as.data.frame(table(study[,c("DecimalDate", "MGE_Replacement")]))
dyn$DecimalDate = as.numeric(as.character(dyn$DecimalDate))
dyn$MGE_Replacement = factor(dyn$MGE_Replacement, levels = c(levels, "Other"))



ggplot(dyn, aes(x = DecimalDate, y = Freq, fill = MGE_Replacement)) + 
  geom_stream(type = "ridge", color = "black", bw = 0.4) +
  scale_fill_manual(values = c(viridis(6), magma(6)[3:5], "grey")) +  
  theme_bw() + 
  theme(legend.position = "right", panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10)) +
  annotate(geom = "rect", xmin = 2016, xmax = 2016.333, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") +
  annotate(geom = "rect", xmin = 2018.5, xmax = 2022, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") +
  labs(x = "Year", y = "Frequency", fill = "Mobile genetic elements") +
  annotate(geom = "text", x = 2015.1, y = 50, label = "Nationwide Surveillance\n10 sites") +
  annotate(geom = "text", x = 2017, y = 120, label = "Nationwide Surveillance\n22 sites") +
  annotate(geom = "text", x = 2022.6, y = 30, label = "2% study\nDhaka") 

ggsave("1b MGEdynamics.svg", width = 10, height = 4)



MGE = study[,c("Sublineage", "MGE_Replacement", "DecimalDate")]
colnames(MGE)[2] = "Gene & MGE presence/absence"
MGE[which(!MGE$`Gene & MGE presence/absence` %in% levels),"Gene & MGE presence/absence"] = "Other"
MGE$`Gene & MGE presence/absence`  = factor(MGE$`Gene & MGE presence/absence`, levels = c(levels, "Other"))
sublineage_levels = c(sort(unique(MGE$Sublineage))[grep("BD2", sort(unique(MGE$Sublineage)))], sort(unique(MGE$Sublineage))[grep("BD1", sort(unique(MGE$Sublineage)))])
MGE$Sublineage = factor(MGE$Sublineage, levels = rev(sublineage_levels))
labels = as.data.frame(table(MGE$Sublineage))
MGE = merge(MGE,labels,by.x= "Sublineage", by.y = "Var1", all.x = T)
MGE$Freq = paste("n =", MGE$Freq)
MGE[which(duplicated(MGE$Sublineage)),"Freq"] = NA

ggplot(MGE, aes(x = DecimalDate, y = Sublineage,  fill = `Gene & MGE presence/absence`, color = `Gene & MGE presence/absence`)) +  geom_density_ridges(bandwidth = 0.05,alpha = 0.8) + 
  theme_bw() + 
  annotate(geom = "segment", y = "BD2.030", yend = "BD2.030", x = min(MGE$DecimalDate, na.rm = T), xend = max(MGE$DecimalDate, na.rm = T)) +
  scale_fill_manual(values = c(viridis(6), magma(6)[3:5], "grey")) + 
  scale_color_manual(values = colorspace::darken(c(viridis(6), magma(6)[3:5], "grey"))) + 
  labs(y = "Sublineage", x = "Year") +
  theme(legend.position = "bottom", panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10)) + guides(color = F, point_color = F) +
  geom_text(aes(label = Freq),x = 2023.9, color = "black", nudge_y = 0.5, hjust = "right") +
  geom_point(data = MGE[which(MGE$Sublineage == "BD2.030"),]) +
  annotate(geom = "rect", xmin = 2016, xmax = 2016.333, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") +
  annotate(geom = "rect", xmin = 2018.5, xmax = 2022, ymin = -Inf, ymax = Inf, alpha = 0.9, color = NA, fill = "white") 

table(metadata[which(metadata$Country == "Bangladesh" & metadata$Reference == "This manuscript" & metadata$Superintegron == " ΔVC_A0455-VC_A0459"),"Sublineage"])

ggsave("Figure S6 MGE Dynamics.svg", width = 10, height = 8)


#### Dissemination from Dhaka and Chittagong (Fig. 1c) ####

#For each of the BD1 and BD2 trees, find the location of each node's parent
for(i in 1:nrow(bd1$annotation)){
  bd1$annotation[i,"Parent_Division"] = bd1$annotation[which(bd1$annotation[,"node"]  == bd1$annotation[i,"parent"]),"Division"]}
for(i in 1:nrow(bd2$annotation)){
  bd2$annotation[i,"Parent_Division"] = bd2$annotation[which(bd2$annotation[,"node"]  == bd2$annotation[i,"parent"]),"Division"]}
bd1$annotation$Lineage = "sBD1"
bd2$annotation$Lineage = "BD2"

#Subtree for each population of interest
#1) sBD1.070
sBD1.070_MRCA = getMRCA(bd1$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "sBD1.070"),"ID"])
sBD1.070 = tree_subset(bd1$tree, sBD1.070_MRCA, levels_back = 0)
sBD1.070 = merge(as.data.frame(fortify(sBD1.070)),bd1$annotation[,c(1,10:ncol(bd1$annotation))], by = "label")
sBD1.070$Category = "sBD1.070"

#2) PLE- BD2.032
PLEneg_BD2.032_MRCA =  getMRCA(bd2$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.032" & metadata$PLE == "None"),"ID"])
PLEneg_BD2.032 = tree_subset(bd2$tree, PLEneg_BD2.032_MRCA, levels_back = 0)
#Drop tips which have also lost ICE-TET
PLEneg_BD2.032 = drop.tip(PLEneg_BD2.032, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.032" & metadata$PLE == "None" & metadata$SXT == "None"),"ID"])
PLEneg_BD2.032 = merge(as.data.frame(fortify(PLEneg_BD2.032)),bd2$annotation[,c(1,10:ncol(bd2$annotation))], by = "label")
PLEneg_BD2.032$Category = "PLE- BD2.032"

#3) SXT- BD2.032
SXTneg_BD2.032_MRCA =  getMRCA(bd2$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.032" & metadata$PLE == "None" & metadata$SXT == "None"),"ID"])
SXTneg_BD2.032 = tree_subset(bd2$tree, SXTneg_BD2.032_MRCA, levels_back = 0)
SXTneg_BD2.032 = merge(as.data.frame(fortify(SXTneg_BD2.032)),bd2$annotation[,c(1,10:ncol(bd2$annotation))], by = "label")
SXTneg_BD2.032$Category = "SXT- BD2.032"

#4) PLE1- in BD2.029
PLEneg_BD2.029_MRCA =  getMRCA(bd2$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.029" & metadata$PLE == "None"),"ID"])
PLEneg_BD2.029 = tree_subset(bd2$tree, PLEneg_BD2.029_MRCA, levels_back = 0)
PLEneg_BD2.029 = merge(as.data.frame(fortify(PLEneg_BD2.029)),bd2$annotation[,c(1,10:ncol(bd2$annotation))], by = "label")
PLEneg_BD2.029$Category = "PLE- BD2.029"

#5) ddmAB in BD2.026
ddmABneg_BD2.026_MRCA = getMRCA(bd2$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.026" & metadata$VSP_II == "VC_0490-interrupted ΔVC_0491-VC_0498"),"ID"])
ddmABneg_BD2.026 = tree_subset(bd2$tree, ddmABneg_BD2.026_MRCA, levels_back = 0)
#Drop tips belonging to other descendent sublineages
ddmABneg_BD2.026 = drop.tip(ddmABneg_BD2.026, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage != "BD2.026"),"ID"])
ddmABneg_BD2.026 = merge(as.data.frame(fortify(ddmABneg_BD2.026)),bd2$annotation[,c(1,10:ncol(bd2$annotation))], by = "label")
ddmABneg_BD2.026$Category = "ΔVC0491-VC0498 BD2.026"

#6) VCA0455-VCA0459 in BD2.028
superintegron_BD2.028_MRCA = getMRCA(bd2$tree, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage == "BD2.028"  & metadata$Superintegron == " ΔVC_A0455-VC_A0459"),"ID"])
superintegron_BD2.028 = tree_subset(bd2$tree, superintegron_BD2.028_MRCA, levels_back = 0)
#Drop tips belonging to other descendent sublineages
superintegron_BD2.028 = drop.tip(superintegron_BD2.028, metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018" & metadata$Sublineage != "BD2.028"),"ID"])
superintegron_BD2.028 = merge(as.data.frame(fortify(superintegron_BD2.028)),bd2$annotation[,c(1,10:ncol(bd2$annotation))], by = "label")
superintegron_BD2.028$Category =  "ΔVCA0455-VCA0459 BD2.028"

#Combine into one data frame
dissemination = rbind.fill(sBD1.070, PLEneg_BD2.032, SXTneg_BD2.032, PLEneg_BD2.029, ddmABneg_BD2.026, superintegron_BD2.028)
dissemination = dissemination[which(dissemination$Parent_Division != dissemination$Division & dissemination$Division %in% metadata$Location & dissemination$Parent_Division %in% metadata$Location),]
dissemination = dissemination[order(dissemination$numeric.date),]

#Number and first date of transitions
dissemination$numeric.date = as.numeric(dissemination$numeric.date)
transitions = merge(dissemination[which(!duplicated(paste(dissemination$Division, dissemination$Category))),], as.data.frame(table(dissemination[,c("Category", "Parent_Division", "Division")])), by = c("Category","Division", "Parent_Division") )

#Get x and y co-ordinates for Bangladesh divisions
division = bangladesh::get_map("division")
centroid = st_drop_geometry(st_centroid(division) %>%dplyr::mutate(lon = sf::st_coordinates(.)[,1],lat = sf::st_coordinates(.)[,2]))
transitions = merge(transitions, centroid, by = "Division")
transitions = merge(transitions, centroid, by.x = "Parent_Division", by.y = "Division")
transitions$Transition = paste(transitions$Parent_Division, transitions$Division, sep = "-> ")
transitions = transitions[order(transitions$numeric.date),]

#Time relative to first transition
for(i in unique(transitions$Category)){transitions[which(transitions$Category == i),"relative.date"] = transitions[which(transitions$Category == i),"numeric.date"] - min(transitions[which(transitions$Category == i),"numeric.date"])}
transitions$relative.date = transitions$relative.date*12

#Scale number of transmission events
for(i in unique(transitions$Category)){transitions[which(transitions$Category == i),"Proportion"] = transitions[which(transitions$Category == i),"Freq"]/sum(transitions[which(transitions$Category == i),"Freq"])}


#Categorise by transmission pattern
transitions$Source = "From Dhaka"
transitions[which(transitions$Category %in% c("PLE- BD2.029","ΔVCA0455-VCA0459 BD2.028")),"Source"] = "From Chittagong"


ggplot(transitions) + 
  geom_sf(data = division) + 
  geom_sf_text(data = division, aes(label = Division), size = 3, vjust = -0.5) + 
  geom_curve(data = transitions[which(transitions$Category %in% unique(transitions$Category)[1:2]),],aes(yend = lat.y, xend = lon.y, y = lat.x, x = lon.x, size = Proportion, color = Category, alpha = relative.date), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), curvature = 0.1) +
  geom_curve(data = transitions[which(transitions$Category %in% unique(transitions$Category)[c(3,6)]),],aes(yend = lat.y, xend = lon.y, y = lat.x, x = lon.x, size = Proportion, color = Category, alpha = relative.date), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), curvature = 0.2) +
  geom_curve(data = transitions[which(transitions$Category %in% unique(transitions$Category)[4]),],aes(yend = lat.y, xend = lon.y, y = lat.x, x = lon.x, size = Proportion, color = Category, alpha = relative.date), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), curvature = 0.3) +
  geom_curve(data = transitions[which(transitions$Category %in% unique(transitions$Category)[5]),],aes(yend = lat.y, xend = lon.y, y = lat.x, x = lon.x, size = Proportion, color = Category, alpha = relative.date), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), curvature = 0.4) +
  facet_wrap(~Source)+ 
  scale_alpha(range =c(1, 0.1)) + 
  scale_size_continuous(range= c(0.2,2)) + 
  theme_void() + 
  scale_color_brewer(palette= "Set2") +
  labs(color = "Date of first transition", size = "Number of transitions") + 
  labs(alpha = "Months relative to first transmission event", size = "Relative number of transmission events", color = "Subpopulation") 

ggsave("Figure 1c.svg", width = 14, height = 5)



#### Sublineage dynamics compared to India (Fig. 2a) ####

#only display Sublineages with 10 samples or more to avoid too many categories
show_lineages = names(which(table(metadata[which(((metadata$Site_y > 25 & metadata$Country == "India")|metadata$Country == "Bangladesh" | metadata$Location == "Kolkata") & (metadata$Year > 2003)),"Sublineage"]) >= 10))
show_lineages = show_lineages[-grep("Singleton", show_lineages)]

northindia = as.data.frame(table(metadata[which(metadata$Site_y > 25 & metadata$Country == "India" & metadata$Year > 2003 & metadata$Sublineage %in% show_lineages),c("Year", "Sublineage")]))
kolkata = as.data.frame(table(metadata[which(metadata$Location == "Kolkata"  & metadata$Year > 2003  & metadata$Sublineage %in% show_lineages),c("Year", "Sublineage")]))
bangladesh = as.data.frame(table(metadata[which(metadata$Country == "Bangladesh" & metadata$Year > 2003  & metadata$Sublineage %in% show_lineages),c("Year", "Sublineage")]))

#Create colour scheme
colorscheme = data.frame(Sublineage = sort(show_lineages), color = NA)
colorscheme[grep("sBD1", colorscheme$Sublineage), "color"] = colorRampPalette(brewer.pal(9, "RdPu"))(length(colorscheme[grep("sBD1", colorscheme$Sublineage), "color"]))
colorscheme[grep("BD2", colorscheme$Sublineage), "color"] = colorRampPalette(brewer.pal(9, "GnBu"))(length(colorscheme[grep("BD2", colorscheme$Sublineage), "color"]))
colorscheme[grep("Other", colorscheme$Sublineage), "color"] = colorRampPalette(brewer.pal(9, "Greys"))(length(colorscheme[grep("Other", colorscheme$Sublineage), "color"]))

#Plot for each region
bangladesh$Year = as.numeric(as.character(bangladesh$Year))
a = ggplot(bangladesh, aes(x = Year, y = Freq, fill = Sublineage)) + 
  geom_stream(type = "proportional", color = "black") + scale_fill_manual(values = colorscheme[which(colorscheme$Sublineage %in% bangladesh$Sublineage), "color"]) + theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_blank()) + facet_grid(~"Bangladesh") + scale_y_continuous(breaks=c(0,0.5, 1)) + scale_x_continuous(breaks=c(2005, 2010, 2015, 2020), limits = c(2003, 2024)) + labs(y = "Proportion") 

northindia$Year = as.numeric(as.character(northindia$Year))
b = ggplot(northindia, aes(x = Year, y = Freq, fill = Sublineage)) + 
  geom_stream(type = "proportional", color = "black") + scale_fill_manual(values = colorscheme[which(colorscheme$Sublineage %in% northindia$Sublineage), "color"]) + theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_blank()) + facet_grid(~"North India") + scale_y_continuous(breaks=c(0,0.5, 1)) + scale_x_continuous(breaks=c(2005, 2010, 2015, 2020), limits = c(2003, 2024)) + labs(y = "Proportion") 

kolkata$Year = as.numeric(as.character(kolkata$Year))
c = ggplot(kolkata, aes(x = Year, y = Freq, fill = Sublineage)) + 
  geom_stream(type = "proportional", color = "black") + scale_fill_manual(values = colorscheme[which(colorscheme$Sublineage %in% kolkata$Sublineage), "color"]) + theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_blank()) + facet_grid(~"Kolkata") + scale_y_continuous(breaks=c(0,0.5, 1)) + scale_x_continuous(breaks=c(2005, 2010, 2015, 2020), limits = c(2003, 2024)) + labs(y = "Proportion") 

g = arrangeGrob(a,b,c, ncol = 1)
ggsave("Figure 2a.svg", g, height = 6, width = 7)

#Create legend
ggplot(colorscheme, aes(x = 1, y = 1, fill = Sublineage)) + geom_tile(color = "black") + scale_fill_manual(values = colorscheme$color) + theme_void() 
ggsave("Figure 2a Legend.svg",height = 6, width = 7)


#### Dissemination from India (Fig. 2b and 2c) ####

#For each node, find the inferred country of the parent node
vc_clean$annotation$parentcountry = NA
for(i in 1:nrow(vc_clean$annotation)){vc_clean$annotation[i,"parentcountry"] = vc_clean$annotation[which(vc_clean$annotation$node == vc_clean$annotation[i,"parent"]),"Country"]}

south_asia = c("Bangladesh", "India","India/Nepal","Nepal","Pakistan", "Sri Lanka")
not_south_asia = unique(metadata[-which(metadata$Country %in% south_asia),"Country"])

#Export events
export = vc_clean$annotation[which(vc_clean$annotation$parentcountry != vc_clean$annotation$Country),]
#For each export event, calculate the number of samples descended from that export event that fall in a) The country exported to b) Any country other than the country of origin
export$Downstream_country = NA
export$Downstream_total = NA
for(i in 1:nrow(export)){
  print(i)
  parent = export[i,"parentcountry"]
  child = export[i,"Country"]
  des = getDescendants(vc_clean$tree,export[i,"node"])
  des = vc_clean$annotation[which(vc_clean$annotation$node %in% des),]
  #If re-enters parent country and spreads from there, exclude from this
  if(length(des[which(des$Country == parent),"node"] > 0)){
    exclude = c()
    for(x in des[which(des$Country == parent),"node"]){
      exclude = unique(c(exclude, getDescendants(vc_clean$tree, x)))
    }
    des = des[-which(des$node %in% exclude),]}
  
  des = des[which(des$isTip == T),]
  
  export[i,"Downstream_total"] = nrow(des)
  export[i,"Downstream_country"] = nrow(des[which(des$Country == child),])
}

#Get GPS co-ordinates of centroid of each country
gps = CoordinateCleaner::countryref
gps = gps[which(!duplicated(gps$name) & gps$type == "country"),]
gps = gps[,c("name", "centroid.lon",  "centroid.lat")]
colnames(gps) = c("Location", "X", "Y")

#Make sure country names match up
export[which(export$Country == "Republic of South Africa"),"Country"] = "South Africa"
export[which(export$Country == "Kurdistan"),"Country"] = "Turkey"
export[which(export$Country %in% c("Gaza", "Palestine")),"Country"] = "Palestinian Territories"
export[which(export$Country == "Russian Federation"),"Country"] = "Russia"
export[which(export$parentcountry == "Kurdistan"),"parentCountry"] = "Turkey"

gps[which(gps$Location == "Guinea-Bissau"),"Location"] = "Guinea Bissau"
gps[which(gps$Location == "Congo - Kinshasa"),"Location"] = "Democratic Republic of the Congo"
gps[which(gps$Location == "Myanmar (Burma)"),"Location"] = "Myanmar" 

export = merge(export, gps, by.x = "parentcountry", by.y = "Location")
export = merge(export, gps, by.x = "Country", by.y = "Location")
colnames(export) = gsub("[.]x",".from", colnames(export))
colnames(export) = gsub("[.]y",".to", colnames(export))

for(i in 1:nrow(export)){
  export[i,"Pair"] = paste(sort(as.character(export[i,c("Country", "parentcountry")])), collapse = " ")
}

export$numeric.date = as.numeric(export$numeric.date)

nrow(export[which(export$parentcountry == "India" & export$numeric.date >= 2003),])
sum(export[which(export$parentcountry == "India" & export$numeric.date >= 2003),"Downstream_total"])
sum(export[which(export$parentcountry == "India" & export$numeric.date >= 2003 & export$Country %in% not_south_asia),"Downstream_total"])

nrow(export[which(export$parentcountry == "Bangladesh" & export$numeric.date >= 2003),])
sum(export[which(export$parentcountry == "Bangladesh" & export$numeric.date >= 2003),"Downstream_total"])
sum(export[which(export$parentcountry == "Bangladesh" & export$numeric.date >= 2003 & export$Country %in% not_south_asia),"Downstream_total"])


#Global source of cholera: number of transmission events
order = names(sort(table(vc_clean$annotation[which(vc_clean$annotation$Country != vc_clean$annotation$parentcountry),c("parentcountry")])))
source = as.data.frame((table(vc_clean$annotation[which(vc_clean$annotation$Country != vc_clean$annotation$parentcountry),c("parentcountry", "Sublineage")])))
source$parentcountry = factor(source$parentcountry, levels = order)
ggplot(source, aes(x = parentcountry, y= Freq, fill= Sublineage)) + geom_bar(stat = "identity", color = "black") + theme_minimal() + coord_flip()  + coord_flip() + labs(x = "Country", y = "Export events", fill = "HierBAPS Lineage") + theme(legend.position = "none")

#Global source of cholera: downstream cases
down_stream = aggregate(export$Downstream_total, by = list(parentcountry = export$parentcountry, Sublineage = export$Sublineage), FUN = sum)
colnames(down_stream)[3] = "Freq"

order =aggregate(export$Downstream_total, by = list(parentcountry = export$parentcountry), FUN = sum)
down_stream$parentcountry = factor(down_stream$parentcountry, levels = order[order(order$x),"parentcountry"])
ggplot(down_stream, aes(x = parentcountry, y= Freq, fill= Sublineage)) + geom_bar(stat = "identity", color = "black") + theme_minimal()  + coord_flip() + labs(x = "Country", y = "Tips descended from export event", fill = "HierBAPS Lineage") + theme(legend.position = "none")

#Global source of cholera: downstream cases, 2003-2023
ex_20 = export[which(export$numeric.date > 2003),]
down_stream = aggregate(ex_20$Downstream_total, by = list(parentcountry = ex_20$parentcountry, Sublineage = ex_20$Sublineage), FUN = sum)
colnames(down_stream)[3] = "Freq"

order =aggregate(ex_20$Downstream_total, by = list(parentcountry = ex_20$parentcountry), FUN = sum)
down_stream$parentcountry = factor(down_stream$parentcountry, levels = order[order(order$x),"parentcountry"])
ggplot(down_stream, aes(x = parentcountry, y= Sublineage, fill= Sublineage)) + geom_bar(stat = "identity", color = "black") + theme_minimal()  + coord_flip() + labs(x = "Country", y = "Tips descended from export event (2003-)", fill = "HierBAPS Lineage") + theme(legend.position = "none")

#Map of export events from India and Bangladesh
world = rnaturalearth::ne_countries()
sf::sf_use_s2(FALSE)
world = world[order(world$pop_est, decreasing = T),]
world = world[which(!duplicated(world$sovereignt)),]

bd1bd2 = export[which(export$Lineage %in% c("sBD1", "BD2")),]
bd1bd2 = bd1bd2[which(bd1bd2$numeric.date > 2003),]
bd1bd2 = bd1bd2[which(bd1bd2$parentcountry %in% south_asia & bd1bd2$Country %in% not_south_asia),]

#Change order so that smaller events are shown on top
bd1bd2 = bd1bd2[order(bd1bd2$Downstream_total, decreasing = T),]
bd1bd2$Group = factor(bd1bd2$Sublineage, levels = unique(bd1bd2$Sublineage))
#Fill in missing sublineage (is a singleton)
bd1bd2[which(is.na(bd1bd2$Sublineage)),"Sublineage"] = metadata[which(metadata$ID == bd1bd2[which(is.na(bd1bd2$Sublineage)),"label"]),"Sublineage"]



world_cropped = st_crop(world, y = c(xmin = -80, xmax = 110, ymin = -50, ymax = 90))

ggplot(world_cropped) + 
  geom_sf(fill= "grey", color = NA) + 
  theme_void() + 
  theme(text = element_text(size = 14), legend.position = "bottom")   + scale_size(range = c(0, 6)) +
  scale_color_manual(values = c(brewer.pal(3,"GnBu")[2:3],brewer.pal(9,"RdPu"))) + 
  labs(size = "No. samples descended\nfrom export event", color = "Sub-lineage") +
  geom_curve(data = bd1bd2, aes(x = X.from, xend = X.to, y = Y.from, yend = Y.to, size = Downstream_total, color = Sublineage), curvature = 0.15, arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = bd1bd2, aes(label = floor(bd1bd2$numeric.date), x = (X.from+X.to)/2, y =(Y.from+Y.to)/2))

ggsave("Figure 2b World map arrows.svg", width = 7, height = 9)


#Just South Asia
bd1bd2 = export[which(export$Lineage %in% c("sBD1", "BD2")),]
bd1bd2 = bd1bd2[which(bd1bd2$parentcountry %in% south_asia & bd1bd2$Country %in% south_asia),]
bd1bd2[which(is.na(bd1bd2$Sublineage)),"Sublineage"] = metadata[which(metadata$ID == bd1bd2[which(is.na(bd1bd2$Sublineage)),"label"]),"Sublineage"]

#Significant first instances of PLEs and ctxB7
events = vc_clean$annotation
events = events[order(events$numeric.date),]
events = events[which(!is.na(events$numeric.date)),]
events = events[which(events$numeric.date != "--"),]
events = rbind(events[which(events$PLE == "PLE1" & events$Lineage == "BD2"),][1,],events[which(events$PLE == "PLE11" & events$Lineage == "sBD1"),][1,],events[which(events$ctxB == "ctxB7" & events$Lineage == "sBD1"),][1,])
events = merge(events,gps, by.x = "Country", by.y = "Location")
events$Label = c("+PLE1", "+PLE11", "+ctxB7")
events$Label = paste0(events$Label,"\n", floor(as.numeric(events$numeric.date)))

events$Sublineage =factor(events$Sublineage, levels = unique(sort(c(bd1bd2$Sublineage, events$Sublineage))))
bd1bd2$Sublineage =factor(bd1bd2$Sublineage, levels = unique(sort(c(bd1bd2$Sublineage, events$Sublineage))))

south_asia_map = st_crop(world, y = c(xmin = 65, xmax = 95, ymin = 5, ymax = 33))

ggplot(south_asia_map) + geom_sf(fill= "grey", color = NA) + theme_void() + theme(text = element_text(size = 14), legend.position = "bottom") +
  geom_curve(alpha = 0, data = bd1bd2, aes(x = X.from, xend = X.to, y = Y.from, yend = Y.to, size = Downstream_total, color = Sublineage), curvature = 0.15, arrow = arrow(length = unit(0.01, "npc"), type="closed"), position = position_jitter(width = 1 ,height = 1), color="black")  + 
  scale_size(range = c(0, 4)) + 
  scale_color_manual(values = c(brewer.pal(8, "GnBu"), colorRampPalette(brewer.pal(9, "RdPu"))(15))) +
  labs(size = "No. samples descended\nfrom export event", color = "Sub-lineage") +
  geom_curve(data = bd1bd2, aes(x = X.from, xend = X.to, y = Y.from, yend = Y.to, size = Downstream_total, color = Sublineage), curvature = 0.15, arrow = arrow(length = unit(0.01, "npc"), type="closed"), position = position_jitter(width = 1 ,height = 1)) +
  guides(colour = guide_legend(ncol = 3), size = guide_legend(ncol = 1), fill = F) + scale_fill_manual(values = c(brewer.pal(8, "GnBu"), colorRampPalette(brewer.pal(9, "RdPu"))(15)), drop = F) + 
  geom_label_repel(data = events, aes(x = X, y = Y, label = Label, fill = Sublineage), size = 5)  

ggsave("Figure 2c South Asia.svg", width = 7, height = 9)




#### Association of different genomic features with disease severity (Fig. 3a) ####
metadata_all$DecimalDate = decimal_date(metadata_all$Date_Month)
for(i in 1:ncol(metadata_all)){
  metadata_all[which(metadata_all[,i] == ""),i] = NA
}

severity = metadata_all[which(metadata_all$Country == "Bangladesh" & metadata_all$Year > 2013 & !is.na(metadata_all$DecimalDate) & metadata_all$Lineage %in% c("BD2", "sBD1")),c("DecimalDate", "VC_0258.assembled", "Superintegron","Percent_ICP1", "Plasmid", "VSP_II", "K139", "Dehydration", "Stool.Nature", "SXT", "PLE", "Lineage", "Site", "Location", "VC_A0219.assembled", "cas1.assembled", "odn.assembled", "SXT_phenotype", "DO_phenotype", "TE_phenotype", "VC_0490.assembled", "VC_0491.assembled", "VC_0492.assembled", "VC_A0455.assembled", "Excluded")]

#Convert to binary
severity$SevereDehydration =  as.numeric(severity$Dehydration == "Severe")
severity$Ricewaterstool =  as.numeric(severity$Stool.Nature == "Rice watery")
severity$BD1 = as.numeric(severity$Lineage == "sBD1")
severity$DdmB = as.numeric(severity$VC_0491.assembled == "yes")
severity$DdmA = as.numeric(severity$VC_0492.assembled == "yes")
severity$pSA7G1 = as.numeric(severity$Plasmid == "pSA7G1")
severity$hlyA = 0
severity[grep("yes", severity$VC_A0219.assembled),"hlyA"] = 1
severity$VCA0455 = 0
severity$VCA0455[grep("yes", severity$VC_A0455.assembled)] = 1
severity$ICETET = as.numeric(severity$SXT == "ICE-TET")
severity$PLE1 = as.numeric(severity$PLE == "PLE1")
severity$PLE11 = as.numeric(severity$PLE == "PLE11")
severity$pSA7G1 = as.numeric(severity$Plasmid == "pSA7G1")
severity$K139 = as.numeric(severity$K139== "Yes")
severity$wbeT = 0
severity[grep("yes", severity$VC_0258.assembled),"wbeT"] = 1
severity$ICP1 = as.numeric(severity$Percent_ICP1 > 0)

#Run analysis controlling for date and site
logreg = as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(logreg) = c("Outcome","Factor", "pvalue", "OR", "UL", "LL")
for(outcome in c("Ricewaterstool", "SevereDehydration")){
  subset = severity[,c(outcome, "hlyA","VCA0455","ICETET","PLE1","PLE11","pSA7G1","K139","wbeT","BD1","DdmB","DdmA","DecimalDate", "Site")]
  colnames(subset)[1] = "outcome"
  model = glm(outcome ~., family = "binomial", data = subset)
  df = data.frame(Factor = row.names(summary(model)$coefficients), pvalue = summary(model)$coefficients[,"Pr(>|z|)"], estimate = summary(model)$coefficients[,"Estimate"])
  conf = confint(model)
  df = merge(df, conf, by = "row.names")
  df$OR = exp(df$estimate)
  df$UL = exp(df$`97.5 %`)
  df$LL = exp(df$`2.5 %`)
  df$Outcome = outcome
  df = df[,c("Outcome","Factor", "pvalue", "OR", "UL", "LL")]
  logreg = rbind(logreg, df)}

#Not enough room on the graph to show every site controlled for
logreg = logreg[-grep("Site|Intercept", logreg$Factor),]

#Neaten up labels
logreg$Factor = gsub("yes", "", logreg$Factor, ignore.case = T) %>% gsub("PLEPLE", "PLE", .) %>% gsub("[.]y", "", .) %>% gsub("SXT", "ICE-TET", .) %>% gsub("Lineage", "", .) %>% gsub("Decimal", "", .)  
logreg$Outcome = gsub( "Dehydration", "Severe dehydration", logreg$Outcome) %>% gsub("Stool.Nature", "Rice water stool", .)
logreg$Factor[which(logreg$Factor == "DdmA")] = "DdmA (VC0492)"
logreg$Factor[which(logreg$Factor == "DdmB")] = "DdmB (VC0491)"
logreg$Factor[which(logreg$Factor == "hlyA")] = "hlyA (VCA0219)"
logreg$Factor[which(logreg$Factor == "wbeT")] = "wbeT (VC0258)"
logreg$Outcome = gsub("Ricewaterstool", "Rice water stool", logreg$Outcome )
logreg$Outcome = gsub("SevereSevere dehydration","Severe dehydration", logreg$Outcome )

ggplot(logreg, aes(x = log(OR), y = Factor, fill = Outcome))  + 
  geom_vline(xintercept = 0, color = "grey20", linetype = "dotted") + 
  geom_errorbar(aes(xmin = log(LL), xmax = log(UL)), width = 0.2, position=position_dodge(width = 0.7), color = "black") + 
  geom_point(position=position_dodge(width = 0.7), size = 3, shape = 21, colour = "black", aes(fill = Outcome)) + 
  theme_bw() + scale_fill_brewer() + labs(x = "log(Odds ratio)") + scale_fill_manual(values = c("#1B998B", "#C5D86D")) +
  theme(legend.position = "bottom", text = element_text(size = 14), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_blank(),  panel.border = element_rect(fill = NA)) +
  annotate(geom = "label", y = 13, x = -2.5, label = "Reduced risk", fill = "grey") +
  annotate(geom = "label", y = 13, x = 2.5, label = "Increased risk", fill = "grey") +
  annotate(geom = "point", y = 13.6, x = 0, color = NA) + labs(fill = "Symptom")


#### Vibrio cholerae reads in metagenomes (Fig. 3b) ####
metagenomics = read.table("SupplementaryDocument3.txt", header = T)

#Limit to samples where there are reads mapping to the regions flanking PLE1 and SXT, and ctxB allele (and therefore likely lineage) has been determined
metagenomics_subset = metagenomics[which(metagenomics$PLE1_flanking_region_reads_minimap > 1 & metagenomics$SXT_flanking_region_reads_minimap > 0 & metagenomics$ctxB %in% c("ctxB1", "ctxB7")),]

#Categorise
metagenomics_subset[which(metagenomics_subset$ctxB == "ctxB7"),"Category"] = "ctxB7 (BD1)"
metagenomics_subset[which(metagenomics_subset$ctxB == "ctxB1" & metagenomics_subset$SXT == "ICE-TET" & metagenomics_subset$PLE == "PLE1"),"Category"] = "ctxB1+ (BD2) PLE1+ ICE-TET+"
metagenomics_subset[which(metagenomics_subset$ctxB == "ctxB1" & metagenomics_subset$SXT == "ICE-TET" & metagenomics_subset$PLE == "ND" ),"Category"] = "ctxB1+ (BD2) ICE-TET+"
metagenomics_subset[which(metagenomics_subset$ctxB == "ctxB1" & metagenomics_subset$SXT == "ND" & metagenomics_subset$PLE == "ND"),"Category"] = "ctxB1+ (BD2) ICE-TET-"

metagenomics_subset$Category = factor(metagenomics_subset$Category, levels = c("ctxB1+ (BD2) PLE1+ ICE-TET+", "ctxB1+ (BD2) ICE-TET+", "ctxB1+ (BD2) ICE-TET-", "ctxB7 (BD1)"))

metagenomics_subset$ICP1_present = c("ICP1 absent", "ICP1 present")[as.numeric(metagenomics_subset$ICP1_screen_kraken_viraldb_percent >= 0.1)+1]

metagenomics_subset$VCholerae_Percent = 100*metagenomics_subset$Vibrio_cholerae_kraken_fulldb_reads/metagenomics_subset$Total_kraken_fulldb_reads

aggregate(metagenomics_subset$VCholerae_Percent, by = list(metagenomics_subset$SXT, metagenomics_subset$ICP1_present, metagenomics_subset$ctxB), FUN = "median")

ggplot(metagenomics_subset[which(!is.na(metagenomics_subset$Category)),], aes(x =Category, y = VCholerae_Percent, fill = Category)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, shape = 21, size = 3) + 
  theme_bw()+
  theme(text = element_text(size = 14), panel.grid = element_blank(), panel.grid.minor.x = element_blank()) + labs(x = "Mobile genetic elements", y = "% V. cholerae reads")  + 
  facet_wrap(~ICP1_present, ncol = 2)+ 
  scale_fill_manual(values = c(brewer.pal(9,"Spectral")[4:8], "grey90", "white")) + guides(fill = F) + ylim(0,110) 
ggsave("FigureS Microbiome ICE.svg",  width = 8, height = 4)

#### ICP1 anti-defence switch (Fig. 3c) ####

icp1 = metadata_all[which(metadata_all$Percent_ICP1 > 0.1  & metadata_all$Country == "Bangladesh" & metadata_all$Year >= 2013 ),c("ctxB_type","VC_1451.assembled","VC_A0219.assembled","VC_0841.assembled","VC_1451.assembled","cas1.assembled","odn.assembled","cas3.assembled",colnames(metadata_all)[grep("csy.*assembled", colnames(metadata_all))],"SXT","Percent_ICP1", "Year", "Lineage")]

#Many ICP1+ positive samples were not included in the 7PET tree due to low V. cholerae reads, and therefore not assigned a lineage. However it is clear from their gene & MGE presence/absence profile that ctxB1 samples are BD2 and ctxB7 sBD1
icp1[which(icp1$ctxB_type == "ctxB1"),"Lineage"] = "BD2"
icp1[which(icp1$ctxB_type == "ctxB7"),"Lineage"] = "sBD1"

icp1[which(icp1$cas1.assembled %in% c("yes", "partial")), "Antidefence"] = "CRISPR"
icp1[which(icp1$odn.assembled %in% c("yes", "yes_nonunique")), "Antidefence"] = "Odn"
icp1[which(icp1$odn == "no" & icp1$cas1 == "no"), "Antidefence"] = "None detected"
icp1[which(icp1$odn.assembled == "interrupted"),"Antidefence"] = "Odn (interrupted)"
icp1$SXT  = factor(icp1$SXT, levels = c("ICE-TET", "None", "ICE-GEN"))
icp1$Antidefence = factor(icp1$Antidefence, levels = c("None detected", "CRISPR", "Odn", "Odn (interrupted)"))

ggplot(icp1, aes(x = Year, fill = Antidefence)) + geom_bar(stat = "count", color = "black") + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA)) + labs(x = "Year", y = "Frequency", fill = "Anti-defence system") + scale_fill_manual(values = c("white", brewer.pal(3,"RdPu")[2:3]))

ggplot(icp1[which(icp1$Lineage %in% c("sBD1", "BD2")),], aes(x = Lineage, fill = Antidefence)) + geom_bar(stat = "count", color = "black") + theme_minimal() + scale_fill_manual(values = (brewer.pal(4,"RdPu"))) + theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA)) + labs(x = "Lineage", y = "Frequency", fill = "Anti-defence system") + scale_fill_manual(values = c("white", brewer.pal(3,"RdPu")[2:3]))


#### K139 and DdmABC (Fig. 3d) ####

#Use same dataframe as for Figure 3a, but subset to BD2 only
subset = severity[which(severity$Lineage == "BD2"),c(outcome, "hlyA","VCA0455","ICETET","PLE1","PLE11","pSA7G1","K139","wbeT","DdmB","DdmA","DecimalDate", "Site")]
summary(glm(K139 ~., family = "binomial", data = subset))

bd_only = metadata[which(metadata$Lineage %in% c("sBD1", "BD2") & metadata$Country == "Bangladesh" & metadata$Year > 2013),]

bd_only$DdmABC = "WT"
bd_only[grep("K147fs", bd_only$VC_0490.var),"DdmABC"] = "VC0490-K147fs"
bd_only[grep("VC_0491-ins", bd_only$VSP_II),"DdmABC"] = paste(bd_only[grep("VC_0491-ins", bd_only$VSP_II),"DdmABC"], "VC0491-ins")
bd_only[grep("VC_0491-VC_0498", bd_only$VSP_II),"DdmABC"] = paste(bd_only[grep("VC_0491-VC_0498", bd_only$VSP_II),"DdmABC"], "ΔVC0491-VC0494")
k = table(bd_only[,c("DdmABC", "K139")])
k_percent = data.frame(DdmABC = rownames(k),percent = 100*k[,"Yes"]/(k[,"Yes"]+k[,"None"]))
k_percent$Proportion = paste(k[,2], k[,1]+k[,2], sep = "/")
k_percent = k_percent[which((k_percent$DdmABC != "NA")),]
k_percent$DdmABC = gsub("VC0490", "DdmC", k_percent$DdmABC) %>% gsub("VC0491", "DdmB", .)
k_percent$DdmABC  = gsub("DdmC-K147fs ΔDdmB-VC0494", "DdmC-K147fs ΔDdmAB", k_percent$DdmABC)


ggplot(k_percent, aes(x = DdmABC, y = percent, color = percent == 0)) + 
  geom_bar(stat = "identity", fill = "#1B998B") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), panel.border = element_rect(fill = NA)) + 
  labs(y = "% with K139") + 
  geom_text(aes(label = Proportion), position = position_stack(vjust = 0.6), color = "black") +
  scale_color_manual(values = c("black", "white")) +
  theme(legend.position = "none") +
  annotate(geom = "segment", x = 1.1, xend = 1.9, y = 70.5, yend = 70.5) +
  annotate(geom = "segment", x = 2.1, xend = 2.9, y = 20, yend = 20) + 
  annotate(geom = "text", x = mean(c(1.1,1.9)), y = 72.4, label = expression(paste("3.3x", 10^-6))) + 
  annotate(geom = "text", x = mean(c(2.1,2.9)), y = 22, label = 0.01)

#### Export events (Fig. 3e) ####

#Total number of samples with each PLE, Bangladesh 2003 onwards
total = as.data.frame(table(metadata[which(metadata$Country == "Bangladesh" & metadata$Year >= 2003),"PLE"]))
colnames(total) = c("PLE", "Freq")
total$Category = "Total"

#Export events from Bangladesh (make sure to run code for Figure 2b first)
bd = export[which(export$parentcountry == "Bangladesh" & export$numeric.date >= 2003),]
export_events = as.data.frame(table(bd$PLE))
colnames(export_events) = c("PLE", "Freq")
export_events$Category = "Export events"

#Downstream samples
downstream = aggregate(bd$Downstream_total, by = list(c(PLE = bd$PLE)), FUN = sum)
colnames(downstream ) = c("PLE", "Freq")
downstream$Category = "Samples descended from export events"

#Combine into one dataframe
PLE_export = rbind(total, export_events, downstream)
rm(total, export_events, downstream)
PLE_export$PLE =factor(PLE_export$PLE, levels = c("None", paste0("PLE", 1:11)))
PLE_export$Category = factor(PLE_export$Category , levels = c("Total", "Export events", "Samples descended from export events"))

#Fisher test
#Group different PLEs together
fisher = PLE_export
fisher$PLE = as.character(fisher$PLE)
for(i in unique(fisher$Category)){
  fisher[which(fisher$Category == i & fisher$PLE %in% paste0("PLE", 1:11)),"Freq"] = sum(fisher[which(fisher$Category == i & fisher$PLE %in% paste0("PLE", 1:11)),"Freq"])
  fisher[which(fisher$Category == i & fisher$PLE %in% paste0("PLE", 1:11)),"PLE"] = "PLE"}
fisher = unique(fisher)
fisher = fisher[,c("PLE", "Category", "Freq")]
#Create a matrix 
fisher = dcast(fisher, Category~PLE, value.var = "Freq")
#Total samples vs export events
fisher.test(as.matrix(fisher[1:2,2:3]))$p.value
#Total samples vs samples descended from export events
fisher.test(as.matrix(fisher[c(1,3),2:3]), simulate.p.value = T)$p.value
#Export events vs samples descended from export events
fisher.test(as.matrix(fisher[c(2,3),2:3]))$p.value

ggplot(PLE_export, aes(x = Category, fill = PLE, y = Freq)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_fill_manual(values = c("white", brewer.pal(11, "Spectral")[c(1:4, 9:11)])) +
  geom_text(aes(label = Freq),  position = position_fill(vjust = 0.5)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "bottom") +
  labs(y = "Proportion", x = "Population") + guides(fill=guide_legend(nrow = 2))

#Repeat for wbeT: total number of samples with each PLE, Bangladesh 2003 onwards
metadata$VC_0258 = "Present"
metadata[which(metadata$VC_0258.assembled %in% c("fragmented", "partial")),"VC_0258"] = c("Fragmented/partial")
metadata[which(metadata$VC_0258.assembled == "no"),"VC_0258"] = c("Not present")
total = as.data.frame(table(metadata[which(metadata$Country == "Bangladesh" & metadata$Year >= 2003),"VC_0258"]))
colnames(total) = c("VC_0258", "Freq")
total$Category = "Total"

#Export events
bd = export[which(export$parentcountry == "Bangladesh" & export$numeric.date >= 2003),]
bd$VC0258[which(bd$VC0258 == "fragmented")] = "Fragmented/partial"
bd$VC0258[which(bd$VC0258 == "yes")] = "Present"
export_events = as.data.frame(table(bd$VC0258))
colnames(export_events) = c("VC_0258", "Freq")
export_events$Category = "Export events"

#Downstream samples
downstream = aggregate(bd$Downstream_total, by = list(c(VC_0258 = bd$VC0258)), FUN = sum)
colnames(downstream ) = c("VC_0258", "Freq")
downstream$Category = "Samples descended from export events"

#Combine into one dataframe
wbeT_export = rbind(total, export_events, downstream)
rm(total, export_events, downstream)
wbeT_export$VC_0258 =factor(wbeT_export$VC_0258, levels = c("Not present", "Fragmented/partial", "Present"))
wbeT_export$Category = factor(wbeT_export$Category , levels = c("Total", "Export events", "Samples descended from export events"))

#Fisher tests
fisher = wbeT_export
fisher$VC_0258 = as.character(fisher$VC_0258)
for(i in unique(fisher$Category)){
  fisher[which(fisher$Category == i & fisher$VC_0258 %in% c("Fragmented/partial", "Not present")),"Freq"] = sum(fisher[which(fisher$Category == i & fisher$VC_0258 %in% c("Fragmented/partial", "Not present")),"Freq"])
  fisher[which(fisher$Category == i & fisher$VC_0258 %in% c("Fragmented/partial", "Not present")),"VC_0258"] = "Not complete"}
fisher = unique(fisher)
fisher = dcast(fisher, Category~VC_0258, value.var = "Freq")
#Total vs export events
fisher.test(as.matrix(fisher[1:2,2:3]))$p.value
#Export events vs samples descended from export events
fisher.test(as.matrix(fisher[c(1,3),2:3]), simulate.p.value = T)$p.value
#Total vs samples descended from export events
fisher.test(as.matrix(fisher[c(2,3),2:3]))$p.value

ggplot(wbeT_export, aes(x = Category, fill = VC_0258, y = Freq))  +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = c("white", "grey85", brewer.pal(8, "Dark2")[1])) +
  geom_text(aes(label = Freq),  position = position_fill(vjust = 0.5)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "bottom") + labs(y = "Proportion", x = "Population", fill = "wbeT (VC0258)")+ guides(fill=guide_legend(nrow = 2))

g = arrangeGrob(a,b, ncol = 2)
ggsave("PLE_overrep.svg",g, width = 5, height = 3)






#### Global 7PET tree (Extended data Fig. 1) ####

#Colour branches by whether inferred to occur in India, Bangladesh or other
vc_clean$annotation$GlobalRegion = "Other"
vc_clean$annotation[which(vc_clean$annotation$Country == "India"),"GlobalRegion"] = "India"
vc_clean$annotation[which(vc_clean$annotation$Country == "Bangladesh"),"GlobalRegion"] = "Bangladesh"

#Create and clean a heatmap
metadata$ICP1 = "None"
metadata[which(metadata$Percent_ICP1 > 0.1),"ICP1"] = "ICP1"
metadata[which(metadata$cas1.assembled == "yes"),"ICP1"] = "ICP1: CRISPR+"
metadata[which(metadata$odn.assembled == "yes"),"ICP1"] = "ICP1: odn+"
metadata$Superintegron = trimws(metadata$Superintegron)
metadata$Superintegron = gsub(" _VC"," ΔVC",gsub("^_VC", "ΔVC", metadata$Superintegron))

heatmap = metadata[which(metadata$ID %in% vc_clean$tree$tip.label),c("ID","Lineage","Sublineage","ctxB_type","PLE","Superintegron","SXT","VSP_II", "Plasmid", "K139","ICP1","VC_0258.assembled","VC_0841.assembled","VC_1451.assembled","VC_A0219.assembled")]
rownames(heatmap) = heatmap$ID
heatmap = heatmap[,-1]
#Reformat PLE column so that will be shown in order (PLE1-11)
heatmap[which(heatmap$PLE %in% paste0("PLE", 1:9)),"PLE"] = paste(" ", heatmap[which(heatmap$PLE %in% paste0("PLE", 1:9)),"PLE"])

#Rename categories for each to make more immediately understandable
heatmap[which(!heatmap$Superintegron %in% names(which(table(heatmap$Superintegron)> 100))),"Superintegron"] = "Other"
heatmap$Superintegron = gsub("WT", "Wild-type", heatmap$Superintegron)

heatmap[which(heatmap$VSP_II %in% names(which(table(heatmap$VSP_II) < 100))),"VSP_II"] = "Other"

colnames(heatmap)[1:3] = c("Lineage","Sublineage","ctxB")
colnames(heatmap) = gsub(".assembled", "", colnames(heatmap))

colnames(heatmap)[9] = "K139"


for(i in c("Lineage","Sublineage","ctxB","PLE","Superintegron","SXT", "VSP_II", "Plasmid", "K139")){
  heatmap[,i][which(!is.na(heatmap[,i]))]   =  paste0(i, ": ", heatmap[,i][which(!is.na(heatmap[,i]))])}

for(i in c("VC_0258","VC_0841","VC_1451","VC_A0219")){
  heatmap[,i][which(heatmap[,i] %in% c("partial", "disrupted", "fragmented"))] = c("partial/disrupted/fragmented")
  heatmap[,i][which(!is.na(heatmap[,i]))]   =  paste0("Gene: ", heatmap[,i][which(!is.na(heatmap[,i]))])}

for(i in 1:ncol(heatmap)){
  heatmap[,i][grep("no|None", heatmap[,i])] = "None"}

#Create colour scheme
color_scheme = data.frame(value = sort(unique(unlist(heatmap))), color = NA)
color_scheme[grep("ctxB", color_scheme$value),"color"] = brewer.pal(7, "Set1")[c(1:3,5:7)]
color_scheme[grep("Gene", color_scheme$value),"color"] = c("grey95", "grey85", brewer.pal(8, "Dark2")[1])
color_scheme[grep("ICP1", color_scheme$value),"color"] = brewer.pal(3, "RdPu")
color_scheme[grep("K139", color_scheme$value),"color"] = brewer.pal(8, "Dark2")[1]
color_scheme[grep("Lineage", color_scheme$value),"color"] = c(brewer.pal(9, "GnBu")[5], "Grey", brewer.pal(9, "RdPu")[5])
color_scheme[which(color_scheme$value == "None"),"color"] = "white"
color_scheme[grep("Plasmid", color_scheme$value),"color"] =  brewer.pal(3, "Set2")[1:2]
color_scheme[grep("PLE", color_scheme$value),"color"] = brewer.pal(10, "Spectral")[c(1:2,4,6:10)]
color_scheme[grep("Superintegron", color_scheme$value),"color"] = brewer.pal(5, "Purples")

color_scheme[grep("SXT", color_scheme$value),"color"] = brewer.pal(6, "Set2")
color_scheme[grep("VSP_II", color_scheme$value),"color"] = c("grey", colorRampPalette(brewer.pal(9, "GnBu"))(9))

color_scheme[grep("Sublineage: BD2", color_scheme$value),"color"] = colorRampPalette(brewer.pal(9, "GnBu")[2:9])(length(color_scheme[grep("Sublineage: BD2", color_scheme$value),"color"]))
color_scheme[grep("Sublineage: sBD1", color_scheme$value),"color"] = colorRampPalette(brewer.pal(9, "RdPu")[2:9])(length(color_scheme[grep("Sublineage: sBD1", color_scheme$value),"color"]))
color_scheme[grep("Sublineage: Other", color_scheme$value),"color"] = colorRampPalette(brewer.pal(9, "Greys")[2:9])(length(color_scheme[grep("Sublineage: Other", color_scheme$value),"color"]))
color_scheme[grep("Sublineage:.*S", color_scheme$value),"color"] = "white"

recent = max(as.Date(vc_clean$annotation$date), na.rm = T)

tree = ggtree(vc_clean$tree, aes(color = GlobalRegion), mrsd=recent) %<+% vc_clean$annotation + scale_color_manual(values = c(brewer.pal(9, "GnBu")[c(5,8)], "grey")) + theme(legend.position = "left") + theme_tree2()
gheatmap(tree, heatmap, color = NA, colnames_angle = 90, colnames_offset_y = -450,   font.size = 4, width = 0.5, offset = 10) + theme(legend.position = "none") + scale_fill_manual(
  values = color_scheme$color) 
ggsave("Global.tree.png", height = 10, width = 15)



#### Sublineage persistence (Extended data Fig. 2) ####
ridge = metadata[which(metadata$SurveillanceStudy == "Bangladesh Surveillance 2014-2018"),c("Location", "DecimalDate", "Sublineage")]
ridge$Sublineage = factor(ridge$Sublineage, levels = rev(sort(unique(ridge$Sublineage))))
ridge = ridge[which(!is.na(ridge$Location)),]
#If more than 5 in sub-lineage show density, else show individual points


a = ggplot(ridge[which(ridge$Sublineage %in% names(which(table(ridge$Sublineage)> 5))),], aes(x = DecimalDate, y = Sublineage, point_color = Location, fill =Location, color = Location)) + 
  scale_y_discrete(drop=FALSE) +
  geom_density_ridges(bandwidth = 0.05,alpha = 0.8, rel_min_height = 0.01) + 
  theme_minimal() + 
  geom_jitter(data = ridge[which(ridge$Sublineage %in% names(which(table(ridge$Sublineage) < 6))),], shape = 21, alpha = 0.5, width = 0, height = 0.05) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) + 
  scale_color_manual(values = colorspace::darken(c(brewer.pal(8, "Set2"), "grey"), amount = 0.5)) + labs(y = "Sublineage", x = "Year")  + ggtitle("a.")


#Calculate each sub-lineage lasted for
length_lasted = metadata[which(metadata$Country == "Bangladesh" & metadata$Reference == "This manuscript"),c("Location", "DecimalDate", "Sublineage")]
min = aggregate(length_lasted$DecimalDate, by = list(Division = length_lasted$Location, Sublineage = length_lasted$Sublineage), min)
max = aggregate(length_lasted$DecimalDate, by = list(Division = length_lasted$Location, Sublineage = length_lasted$Sublineage), max)
colnames(min)[3] = "min"
colnames(max)[3] = "max"
divisions = merge(min, max, by=c("Division", "Sublineage"))
#Adjust to when first detected in country
divisions$min_adj = NA
divisions$order = NA
divisions = divisions[order(divisions$Sublineage, divisions$min),]

for(i in unique(divisions$Sublineage)){
  divisions[which(divisions$Sublineage == i),"min_adj"] = divisions[which(divisions$Sublineage == i),"min"]-min(divisions[which(divisions$Sublineage == i),"min"])
  divisions[which(divisions$Sublineage == i),"order"] = 1:nrow(divisions[which(divisions$Sublineage == i),])}

model = lm(min_adj~order, divisions)
newdata = divisions[,-which(colnames(divisions) %in% "min_adj")]
divisions[,c("fit", "lwr", "upr")] = predict(model, newdata = divisions, interval = "confidence")
#Predicted time (in months) to reach all 8 divisions?
divisions[which(divisions$order == 8),c("fit", "lwr", "upr")]*12


b = ggplot(divisions, aes(x = order, y = min_adj,fill = Division, color = Division)) + geom_jitter(shape = 21, size = 3, width = 0.1)   + labs(x = "Number of divisions sublineage detected in", y = "Time since sublineage first detected (years)") + theme_minimal()  + scale_fill_brewer(palette="Set2") + scale_color_manual(values =  colorspace::darken(brewer.pal(8, "Set2"), amount = 0.5)) + ggtitle("b.") + geom_abline(slope = 0.15695, intercept = -0.07049, size = 2)

#How long sub-lineages last in Bangladesh overall?
length_lasted = length_lasted[which(!is.na(length_lasted$DecimalDate)),]
min = aggregate(length_lasted$DecimalDate, by = list(Sublineages = length_lasted$Sublineage), min)
max = aggregate(length_lasted$DecimalDate, by = list(Sublineages = length_lasted$Sublineage), max)
colnames(min)[2] = "min"
colnames(max)[2] = "max"
nationwide = merge(min, max, by=c("Sublineages"))
quantile(12*(nationwide$max - nationwide$min))

#Plot Figure
g = arrangeGrob(a, b)
ggsave("Figure S4.png", g, height = 9, width = 7)


#### Genomic context of genetic changes (Extended data Fig. 4) ####
annotation = as.data.frame(rtracklayer::import("Vibrio_cholerae_O1_biovar_eltor_str_N16961_v2.gff"))
annotation$forward = annotation$strand == "+"
#As two chromosomes, create columns with cumulative position
annotation$MergedStart = annotation$start
annotation$MergedEnd = annotation$end
annotation[which(annotation$seqnames == "AE003853"),"MergedStart"] = annotation[which(annotation$seqnames == "AE003853"),"MergedStart"] + 2961182
annotation[which(annotation$seqnames == "AE003853"),"MergedEnd"] = annotation[which(annotation$seqnames == "AE003853"),"MergedEnd"] + 2961182
annotation$locus_tag = gsub("VC_", "VC", annotation$locus_tag)

#hlyA
hlyA = data.frame(x = c(238320, 238337), y = "", label = c("17bp deletion", ""))

a = ggplot(annotation[which(annotation$start > 235585 & annotation$end < 241810 & annotation$seqnames == "AE003853"),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 238320) + 
  geom_feature(y = "", x = 238337) + 
  geom_feature_label(data = hlyA, aes(x =x, y = y, label = label)) + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. II)") + ggtitle("a. hlyA (VCA0219)")

#VSP-II
minidel = data.frame(x =524676, y = "", label = "1. K147fs")
annotation[which(annotation$locus_tag == "VC0491"),]
transposon = data.frame(x =  525235, y = "", label = "2. Transposon insertion")
#496-498 always deleted in BD2

b = ggplot(annotation[which(annotation$locus_tag %in% paste0("VC0", 489:502)),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag))  + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  labs(x = "Position (N16961 chr. I)") + 
  geom_feature(y = "", x = 524676) + 
  geom_feature_label(data = minidel, aes(x =x-1000, y = y, label = label))  + geom_feature(y = "", x =525235) + 
  geom_feature_label(data = transposon, aes(x =x+1500, y = y, label = label)) +
  annotate(geom = "segment", x = 529100, xend = 532440, y = 1.15, yend = 1.15, color = "black") +
  annotate(geom = "text", x = mean(c(529100,532440)), y = 1.25, label = "Deleted in all BD2", color = "black", size = 3)  +
  annotate(geom = "segment", x = 525300, xend = 534100, y = 0.9, yend = 0.9, color = "black") +
  annotate(geom = "text", x = mean(c(525300,534100)), y = 0.8, label = "3. Deletion", color = "black", size = 3) + ggtitle("b. VSP-II")

b

#wbeT
insertions = data.frame(x =  c(263687, 263685), y = "", label = c("Insertion site (BD2)", "Insertion site (BD1)"))

c = ggplot(annotation[which(annotation$locus_tag %in% paste0("VC0", 257:259)),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() + 
  geom_feature(y = "", x = 263687)+ 
  geom_feature(y = "", x = 263685) + 
  annotate(geom = "text", label = "Insertion site", x =263685, y = 1.4, size = 3) +
  geom_gene_label(align = "left", aes(label = locus_tag)) + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  labs(x = "Position (N16961 chr. I)") + ggtitle("c. wbeT (VC0258)")
c

#Superintegron deletion
d = ggplot(annotation[which(annotation$locus_tag %in% paste0("VCA0", 453:460)),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 520264) +
  annotate(geom = "segment", x = 3368780-2961182, xend = 3369690-2961182, y = 1.15, yend = 1.15) +
  annotate(geom = "segment", x = 3367500-2961182, xend = 3371300-2961182, y = 0.8, yend = 0.8) +
  annotate(geom = "text", x = 408053, y = 1.3, label = "Deleted in BD2")  +
  annotate(geom = "text", x = 408053, y = 0.63, label = "Deleted in BD2 subclade") +
  annotate(geom = "segment", x = 3367444-2961182, xend = 3367567-2961182, y = 1, yend = 1, color = "deeppink3", size = 2) +
  annotate(geom = "segment", x = 3371305-2961182, xend = 3371428-2961182, y = 1, yend = 1, color = "deeppink3", size = 2) + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. II)") +
  annotate(geom = "segment", x = 3367444-2961182, xend = 3367444-2961182, y = 1, yend = 1.5, color = "deeppink3") +
  annotate(geom = "segment", x = 3371428-2961182, xend = 3371428-2961182, y = 1, yend = 1.5, color = "deeppink3")+
  annotate(geom = "segment", x = 3371428-2961182, xend = 3367444-2961182, y = 1.5, yend = 1.5, color = "deeppink3") +
  annotate(geom = "label", x = 408254, y = 1.5, label = "Homologous regions", color = "deeppink3")+ ggtitle("d. Superintegron deletion")
d

#acfC
e = ggplot(annotation[which(annotation$locus_tag %in% paste0("VC0", 840:842)),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 903585) + geom_feature_label(y = "", x = 903585, aes(label = "1bp deletion"))  + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. I)")  + ggtitle("e. acfC (VC0841)")

#rtxA
f = ggplot(annotation[which(annotation$locus_tag %in% paste0("VC", 1450:1452)),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 1550155) + geom_feature_label(y = "", x = 1550155, aes(label = "TGG>TGA(Stop)"))  + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. I)")  + ggtitle("f. rtxA (VC1451)")
f


PLE = read.table("PLE_genes.txt", header = T, fill = T)

#Assign colour by category
colors = unique(PLE[,c("Gene", "Category")])
colors = colors[order(colors$Gene),]
colors$color = NA
colors[which(colors$Category == "Gene found in PLE1"),"color"] = colorRampPalette(brewer.pal(9, "GnBu"))(length(colors[which(colors$Category == "Gene found in PLE1"),"color"]))
colors[which(colors$Category == "Gene found in PLE2-PLE10"),"color"] = colorRampPalette(brewer.pal(9, "RdPu"))(length(colors[which(colors$Category == "Gene found in PLE2-PLE10"),"color"]))
colors[which(colors$Category == "Gene unique to PLE11"),"color"] = "grey"

#Plot genes in PLEs
ggplot(PLE, aes(xmin = start, xmax = end, y = PLE, fill = Gene, forward = orientation)) +
  geom_gene_arrow()  +
  theme_genes() + scale_fill_manual(values = colors$color) + guides(fill = F) + theme(text = element_text(size = 14), axis.title.y = element_blank())

#Plot genomic context
ple1 = data.frame(x =  411334, y = "", label = "PLE1 insertion site")
ggplot(annotation[which(annotation$start > 410334 & annotation$end < 413334 & annotation$seqnames == "AE003853"),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 411334) + 
  geom_feature_label(data = ple1, aes(x =x, y = y, label = label)) + scale_fill_brewer(palette = "Greys") + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. II)")

ple11 = data.frame(x =  520264, y = "", label = "PLE11 insertion site")
ggplot(gannotation[which(gannotation$start > 518264 & gannotation$end < 522264 & gannotation$seqnames == "AE003853"),],  aes(xmin = start, xmax = end, fill = locus_tag, y= "", forward = forward)) + geom_gene_arrow() +
  theme_genes() +
  geom_gene_label(align = "left", aes(label = locus_tag)) + geom_feature(y = "", x = 520264) + 
  geom_feature_label(data = ple11, aes(x =x, y = y, label = label)) + scale_fill_brewer(palette = "Greys") + theme(legend.position = "none", axis.title.y = element_blank()) + labs(x = "Position (N16961 chr. II)")




#### Distribution of sporadic mobile genetic elements (Extended data Fig. 5) ####
division = bangladesh::get_map("division")

pSA7G1 = metadata_all[which(metadata_all$Country == "Bangladesh" & !is.na(metadata_all$SurveillanceStudy) & metadata_all$Plasmid == "pSA7G1"),]
ggplot(pSA7G1) + 
  geom_sf(data = division) + theme_void() +  
  geom_point(aes(x = Site_x, y = Site_y), position = position_jitter(width = 0.1, height = 0.1)) +  ggtitle("a. pSA7G1") + theme(text = element_text(size = 18))

K139 = as.data.frame(table(metadata_all[which(metadata_all$Country == "Bangladesh" & !is.na(metadata_all$SurveillanceStudy)),c("Site", "K139")]))
K139 = dcast(K139, Site~K139)
K139  = merge(K139, unique(metadata_all[,c("Site", "Site_x", "Site_y")]), by = "Site")
colnames(K139)[3] = c("Yes")
ggplot(K139) +
  geom_sf(data = division) + geom_scatterpie(aes(x=Site_x, y=Site_y, group=Site), data=K139, cols=c("None", "Yes"), size = 0.1) + scale_fill_manual(values = c("white", brewer.pal(8, "Dark2")[1])) + theme_void() + labs(fill = "K139 in genome") + ggtitle("b. K139") + theme(legend.position = "bottom") + theme(text = element_text(size = 18))

ICP1 = metadata_all[which(metadata_all$Country == "Bangladesh" & !is.na(metadata_all$SurveillanceStudy) & metadata_all$Percent_ICP1 > 0.1),]
ggplot(ICP1) + 
  geom_sf(data = division) + theme_void() +  
  geom_point(aes(x = Site_x, y = Site_y, shape = Excluded), position = position_jitter(width = 0.1, height = 0.15)) + ggtitle("c. ICP1") + theme(text = element_text(size = 18))

#### Location of samples (Extended data Fig. 6) ####

#Ganges basin
url = "https://datacatalogfiles.worldbank.org/ddh-published/0041426/1/DR0051689/major_basins_of_the_world_0_0_0.zip"
download.file(url, "basins.zip")
unzip("basins.zip")
basin = st_read("Major_Basins_of_the_World.shp")
basin = basin[which(basin$ID == 406),]

#Rivers
url = "https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_as_shp.zip"
download.file(url, "rivers.zip")
unzip("rivers.zip")
rivers= st_read("HydroRIVERS_v10_as_shp")
rivers= rivers[which(rivers$ORD_FLOW %in% 2:4),]
rivers= st_crop(rivers, basin)
rivers= st_simplify(rivers)
rivers$ORD_CLAS = as.numeric(rivers$ORD_CLAS)

#National boundaries
world = ne_countries(scale = "medium", returnclass = "sf")
world = world[which(world$sovereignt %in% c("India", "Nepal", "Myanmar", "Bangladesh", "Bhutan", "China", "Pakistan")),]

#Split up Angermeyer samples inside and outside of Dhaka so show up as different labels
metadata_split = metadata[which(metadata$Country %in% c("India", "Nepal", "Bangladesh")),]
metadata_split[which(metadata_split$Reference == "Angermeyer et al, 2022 mBio" & metadata_split$Site == "Mathbaria Upozilla Health Complex, Pirojpur"),"Country"] = "Bangladesh (Mathbaria)"
#Show different sub-studies of this study separately
metadata_split[which((metadata_split$SurveillanceStudy != "")),"Reference"] = metadata_split[which((metadata_split$SurveillanceStudy != "")),"SurveillanceStudy"] 
sites = unique(metadata_split[which(metadata_split$Country %in% c("Bangladesh", "India", "Nepal","Bangladesh (Mathbaria)")),c("Country", "Reference", "Site_x", "Site_y")])

#Limit to those with more than 10 samples
metadata_split$Category = paste(metadata_split$Country, metadata_split$Reference)
sites = sites[which(paste(sites$Country, sites$Reference) %in% names(which(table(metadata_split$Category) >= 10))),]

#Aggregate label positions by study and country
for(i in unique(paste(sites$Country, sites$Reference))){
  sites[which(paste(sites$Country, sites$Reference) == i),"Label_x"] = median(sites[which(paste(sites$Country, sites$Reference) == i),"Site_x"], na.rm = T)
  sites[which(paste(sites$Country, sites$Reference) == i),"Label_y"] = median(sites[which(paste(sites$Country, sites$Reference) == i),"Site_y"], na.rm = T)}

#Sample size and year range
sites$N = NA
for(i in unique(sites$Reference)){
  for(j in unique(sites[which(sites$Reference == i),"Country"])){
    sites[which(sites$Reference == i & sites$Country == j),"Start"] = min(metadata_split[which(metadata_split$Reference == i & metadata_split$Country == j),"Year"], na.rm = T)
    sites[which(sites$Reference == i  & sites$Country == j),"End"] = max(metadata_split[which(metadata_split$Reference == i & metadata_split$Country == j),"Year"], na.rm = T)
    sites[which(sites$Reference == i & sites$Country == j),"N"] = nrow(metadata_split[which(metadata_split$Reference == i & metadata_split$Country == j),])
  }}

#Standardise reference format
sites[intersect(grep("et al", sites$Reference), grep("2[0-9][0-9][0-9]", sites$Reference)),"Year"] = gsub(",", "", gsub(" .*", "", gsub(".* 20", "20", sites[intersect(grep("et al", sites$Reference), grep("2[0-9][0-9][0-9]", sites$Reference)),"Reference"])))
sites[grep("et al", sites$Reference),"Author"] = trimws(gsub("et al.*", "", sites[grep("et al", sites$Reference),"Reference"]))
sites[grep("et al", sites$Reference),"Journal"] = trimws(gsub("[^[:alnum:] ]", "", gsub('[0-9]+', "", gsub(".*et al", "", sites[grep("et al", sites$Reference),"Reference"]))))
sites[grep("Comm", sites$Journal),"Journal"] = "Nat. Commun."
sites[grep("NatGen", sites$Journal),"Journal"] = "Nat. Genet."
sites[grep("PLOS|pNTD", sites$Journal, ignore.case = T),"Journal"] = "PLOS NTDs"
sites$Reference_reformatted = sites$Reference
sites[grep("et al", sites$Reference),"Reference_reformatted"] = paste(sites[grep("et al", sites$Reference),"Author"], "et al.", sites[grep("et al", sites$Reference),"Year"], sites[grep("et al", sites$Reference),"Journal"])

#Labels
sites$Label = sites$Reference_reformatted
sites[which(sites$Start == sites$End),"Label"] = paste(sites[which(sites$Start == sites$End),"Label"] , sites[which(sites$Start == sites$End),"Start"], sep = "\n")
sites[which(sites$Start != sites$End),"Label"] = paste(sites[which(sites$Start != sites$End),"Label"] , paste(sites[which(sites$Start != sites$End),"Start"], sites[which(sites$Start != sites$End),"End"], sep = "-"), sep = "\n")
sites = sites[order(sites$N, decreasing = T),]
sites$Label = factor(sites$Label, levels = unique(sites$Label))
sites$Label = paste0(sites$Label, " n=", sites$N, sep = " ")

#Combine each of Kolkata and icddr,b into one label
sites[which(sites$Reference %in% unique(metadata[which(metadata$Location == "Kolkata"),"Reference"]) & sites$Country == "India"),"Label"] = paste(unique(sites[which(sites$Reference %in% unique(metadata[which(metadata$Location == "Kolkata"),"Reference"]) & sites$Country == "India"),"Label"]), collapse = "\n")
sites[which(sites$Reference %in% unique(metadata_split[which(metadata_split$Site == "icddr,b"),"Reference"]) & sites$Country == "Bangladesh"),"Label"]  = paste(sites[which(sites$Reference %in% unique(metadata_split[which(metadata_split$Site == "icddr,b"),"Reference"]) & sites$Country == "Bangladesh"),"Label"], collapse = "\n")

map = ggplot(world)  + 
  geom_sf(fill = "white", linewidth = 0) + 
  geom_sf(data = basin, linewidth = 0,  fill = brewer.pal(9, "Blues")[1]) + 
  geom_sf(fill = NA, linewidth = 0.5) + 
  geom_sf(data =rivers, aes(linewidth = ORD_FLOW, color = ORD_FLOW)) + 
  scale_linewidth(range = c(1,0.2)) + 
  scale_color_gradient(high = brewer.pal(9, "Blues")[5], low = brewer.pal(9, "Blues")[9]) + theme_void() + xlim(74,96) + ylim(20,32) + 
  theme(panel.background = element_rect(fill = brewer.pal(9, "Blues")[6])) + 
  new_scale_color() + 
  geom_jitter(data = sites, aes(x = Site_x, y=  Site_y, fill = Label), shape = 21, width = 0.03, height = 0.02, size = 2, alpha = 0.8) + 
  theme(legend.position = "bottom") + scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(13)) + theme(legend.position = "none")

map  + geom_label(data = sites[which(!duplicated(sites$Label)),], aes(x = Label_x, y=  Label_y, label = Label, fill = Label), size = 2, color = "black", alpha = 0.8)

#Frequency over time
metadata_all$Category = "Rest of world"
metadata_all[which(metadata_all$Country == "Bangladesh"),"Category"] = "Bangladesh"
metadata_all[which(metadata_all$Country == "India"), "Category"] = "India"
metadata_all[which(metadata_all$Country == "Nepal"), "Category"] = "Nepal"
metadata_all[which(metadata_all$City %in% c("Dindigul")), "Category"] = "South India"
metadata_all[which(metadata_all$Country == "India" & metadata_all$SurveillanceStudy == "North India study"),"Category"] = "North India"
metadata_all[which(metadata_all$City %in% c("Kolkata", "Assam")), "Category"] = "East India"
metadata_all[which(metadata_all$City == "Delhi"),"Category"] = "North India"
metadata_all[which(metadata_all$Location %in% c("Yammuna Nagar", "Patiala", "Panchkula", "Morinda", "Mohali", "Hallomajra", "Derabassi", "Chandigarh", "Ambala")),"Category"] = "North India"
metadata_all[which(metadata_all$Category %in% c("India", "South India")),"Category"] = "Rest of world"
samples = as.data.frame(table(metadata_all[,c("Year", "Category")]))
samples = samples[which(samples$Freq > 0),]
samples$Category = factor(samples$Category, levels = rev(sort(unique(samples$Category))))
ggplot(samples, aes(x = as.numeric(as.character(Year)), y = Category)) + geom_point(aes(size = Freq), alpha = 0.8) + theme_bw() + labs(x = "Year", y = "Region")




#### Sharing between neighbouring regions in different countries (Extended data Fig. 7) ####
#Khulna
table(metadata[which(metadata$Location == "Khulna"),"Sublineage"] %in% metadata[which(metadata$Location == "Kolkata"),"Sublineage"])
table(metadata[which(metadata$Location == "Khulna"),"Sublineage"] %in% metadata[which(metadata$Location != "Kolkata" & metadata$Country == "Bangladesh"),"Sublineage"])

#Assam
metadata[which(metadata$Location == "Assam"),c("ID","Sublineage")]
metadata[which(metadata$Sublineage == "Other.159"),c("Country", "Location", "Year")]

#Nepal
table(metadata[which(metadata$Country == "Nepal" & metadata$Year == 2010),"Sublineage"] %in% metadata[which(metadata$Country == "Bangladesh"),"Sublineage"])
table(metadata[which(metadata$Country == "Nepal" & metadata$Year == 2010),"Sublineage"] %in% metadata[which(metadata$Country == "India"),"Sublineage"])


gps = data.frame(location = c("Nepal","Bangladesh","India","Kolkata","Assam","Khulna"), x = c(84,90,77.1,88.4,93.2,89.3), y = c(28,24,28.7,22.6,26.4,22.9))

#Nepal in 2010 vs India and Bangladesh
df_a = metadata[which(metadata$Year == 2010 & metadata$Country %in% c("Nepal","India", "Bangladesh")),]
df_a$Country = factor(df_a$Country, levels = c( "India", "Nepal","Bangladesh"))
df_a = dcast(as.data.frame(table(df_a[,c("Country", "Sublineage")])), Country~Sublineage)
df_a = merge(df_a, gps, by.x = "Country", by.y = "location")

#Add segments_a for shared sublineages
sublin = unique(metadata[which(metadata$Year == 2010 & metadata$Country %in% c("Nepal","India", "Bangladesh")), c("Sublineage", "Country")])
sublin = sublin[-grep("Singleton", sublin$Sublineage),]

segments_a = data.frame(Country1 = c("Nepal", "India", "Bangladesh"), Country2 = c("India", "Bangladesh", "Nepal"), Shared = NA)
for(i in 1:nrow(segments_a)){
  segments_a[i,"Shared"] = length(intersect(sublin[which(sublin$Country == segments_a[i,"Country1"]), "Sublineage"], sublin[which(sublin$Country == segments_a[i,"Country2"]), "Sublineage"])) }
segments_a = merge(segments_a, gps, by.x = "Country1", by.y = "location")
segments_a = merge(segments_a, gps, by.x = "Country2", by.y = "location")

world = ne_countries(scale = "medium", returnclass = "sf")
world = world[which(world$sovereignt %in% c("India", "Nepal", "Myanmar", "Bangladesh", "Bhutan", "China", "Pakistan")),]

a = ggplot(world)  +
  geom_sf(data = world)  + 
  theme_void()  + 
  theme(text = element_text(size = 18)) + 
  xlim(74,96) + ylim(20,32) +
  geom_segment(data = segments_a, aes(x = x.x, xend = x.y, y = y.x, yend = y.y, linewidth = Shared)) + 
  geom_scatterpie(aes(x=x, y=y, group=Country), data=df_a, cols=colnames(df_a)[2:12], size = 0.2, pie_scale =7) + 
  scale_linewidth(range = c(0,4), limits = c(0,9)) + 
  labs(linewidth = "No. sub-lineages shared", fill = "Sublineage") 
a
ggsave(a, file = "Sublineagesharedlegend.svg")



assam_subtree = tree_subset(vc$tree, getMRCA(vc$tree, metadata[which(metadata$Location == "Assam"),"ID"]), levels_back = 1)
ggtree(assam_subtree)
heatmap = metadata[which(metadata$ID %in% assam_subtree$tip.label),c("ID","Country", "Location")]
rownames(heatmap) = heatmap$ID
heatmap[which(heatmap$Country == "India" & heatmap$Location != "Assam"),"Location"] = "India (Other)"
heatmap[which(heatmap$Country == "India" & heatmap$Location == "Assam"),"Location"] = "India (Assam)"
heatmap[which(heatmap$Country == "Nepal"),"Location"] = "Nepal"
hmap = data.frame(location = heatmap$Location)
rownames(hmap) = heatmap$ID
gheatmap(ggtree(assam_subtree), hmap, color = NA) + scale_fill_manual(values = c("#868686", "#0B68AC", "red"))
ggsave("assam_subtree.svg", height = 4, width = 4)


#Assam vs North India and Bangladesh
metadata$Category = NA
metadata[which(metadata$Country == "Bangladesh"), "Category"] = "Bangladesh"
metadata[which(metadata$Country == "India" & metadata$Site_x > 25), "Category"] = "India" #Specifically North India
metadata[which(metadata$Location == "Kolkata"), "Category"] = "Kolkata"
metadata[which(metadata$Location == "Assam"), "Category"] = "Assam"

df_b = metadata[which(metadata$Year %in%  2002:2005),]
df_b = dcast(as.data.frame(table(df_b[,c("Category", "Sublineage")])), Category~Sublineage)
df_b = df_b[,-which(colnames(df_b) %in% names(which(colSums(df_b[,2:ncol(df_b)]) == 0)))]

#Function for unique combinations
combinations = function(x){
  all = as.data.frame(matrix(ncol = 2, nrow = 0))
  for(i in x){all = rbind(all,  df_b = data.frame(V1 = i, V2 = x))}
  for(i in 1:nrow(all)){all[i,c("V1", "V2")] = sort(as.character(all[i,c("V1", "V2")]))}
  all = unique(all[which(all$V1 != all$V2),])
  return(all)}

#Merge
df_b = merge(df_b, gps, by.x = "Category", by.y  = "location")
#Shift GPS co-ordinates for Bangladesh up from the centre, so that can see the line joining Kolkata and Assam better
df_b[which(df_b$Category == "Bangladesh"),"y"] = df_b[which(df_b$Category == "Bangladesh"),"y"] + 1

#Add a segment for shared
sublin = unique(metadata[which(metadata$Year %in%  2002:2005 & !is.na(metadata$Category)),c("Sublineage", "Category")])
sublin = sublin[-grep("Singleton", sublin$Sublineage),]
segments_b = combinations(no_na(unique(metadata$Category)))
colnames(segments_b) = c("Category1", "Category2")
segments_b$Shared = NA  
for(i in 1:nrow(segments_b)){
  segments_b[i,"Shared"] = length(intersect(sublin[which(sublin$Category == segments_b[i,"Category1"]), "Sublineage"], sublin[which(sublin$Category == segments_b[i,"Category2"]), "Sublineage"])) }

segments_b = merge(segments_b, gps, by.x = "Category1", by.y = "location")
segments_b = merge(segments_b, gps, by.x = "Category2", by.y = "location")

b = ggplot(world)  +
  geom_sf(data = world)  + 
  theme_void()  + 
  theme(text = element_text(size = 18)) + 
  xlim(74,96) + ylim(20,32) +
  geom_segment(data = segments_b, aes(x = x.x, xend = x.y, y = y.x, yend = y.y, linewidth = Shared)) + 
  scale_linewidth(range = c(0,4), limits = c(0,9)) + 
  labs(linewidth = "No. sub-lineages shared", fill = "Sublineage") +
  geom_scatterpie(aes(x=x, y=y, group=Category), data=df_b, cols=colnames(df_b)[2:59], size = 0.2, pie_scale = 4) 
b

#Khulna vs Kolkata and rest of Bangladesh
metadata$Category = NA
metadata[which(metadata$Country == "Bangladesh"), "Category"] = "Bangladesh"
metadata[which(metadata$Location == "Khulna"), "Category"] = "Khulna"
metadata[which(metadata$Location == "Kolkata"), "Category"] = "Kolkata"

df_c = metadata[which(metadata$Year %in%  2015:2018),]
graph_c_sublineages = unique(df_c$Sublineage)
df_c = dcast(as.data.frame(table(df_c[,c("Category", "Sublineage")])), Category~Sublineage)
df_c = df_c[,-which(colnames(df_c) %in% names(which(colSums(df_c[,2:ncol(df_c)]) == 0)))]
df_c = merge(df_c, gps, by.x = "Category", by.y = "location")

#Add a segment for shared
sublin = unique(metadata[which(metadata$Year %in%  2015:2018 & metadata$Category %in% c("Khulna", "Bangladesh", "Kolkata")),c("Sublineage", "Category")])
sublin = sublin[-grep("Singleton", sublin$Sublineage),]
segments_c = data.frame(Category1 = c("Khulna", "Khulna"), Category2 = c("Bangladesh", "Kolkata"), Shared = NA)

for(i in 1:nrow(segments_c)){
  segments_c[i,"Shared"] = length(intersect(sublin[which(sublin$Category == segments_c[i,"Category1"]), "Sublineage"], sublin[which(sublin$Category == segments_c[i,"Category2"]), "Sublineage"])) }

segments_c = merge(segments_c, gps, by.x = "Category1", by.y = "location")
segments_c = merge(segments_c, gps, by.x = "Category2", by.y = "location")

c = ggplot(world)  +
  geom_sf(data = world)  + 
  theme_void()  + 
  theme(text = element_text(size = 18)) + 
  xlim(85, 92) + ylim(21, 25) + 
  geom_segment(data = segments_c, aes(x = x.x, xend = x.y, y = y.x, yend = y.y, linewidth = Shared)) + 
  scale_linewidth(range = c(0,4), limits = c(0,6))  + 
  labs(linewidth = "No. sub-lineages shared", fill = "Sublineage") +
  geom_scatterpie(aes(x=x, y=y, group=Category), data=df_c, cols=colnames(df_c)[2:24], size = 0.2, pie_scale = 15) 
c

# Unified sublineage color scheme
sublineages = unique(sort(c(colnames(df_a), colnames(df_b), colnames(df_c))))
sublineages = sublineages[grep("BD|Other", sublineages)]
sublin_col = data.frame(Sublineages = sublineages, col = NA)
sublin_col[grep("sBD1", sublin_col$Sublineages),"col"] = colorRampPalette(brewer.pal(9, "RdPu"))(length(sublin_col[grep("sBD1", sublin_col$Sublineages),"col"]))
sublin_col[grep("BD2", sublin_col$Sublineages),"col"] = colorRampPalette(brewer.pal(9, "GnBu"))(length(sublin_col[grep("BD2", sublin_col$Sublineages),"col"]))
sublin_col[grep("Other", sublin_col$Sublineages),"col"] = colorRampPalette(brewer.pal(9, "Greys"))(length(sublin_col[grep("Other", sublin_col$Sublineages),"col"]))

ggsave("Figure S13_scale.svg", a, height = 3, width = 9)

a = a + scale_fill_manual(values = sublin_col[which(sublin_col$Sublineages %in% colnames(df_a)),"col"])

b = b + scale_fill_manual(values = sublin_col[which(sublin_col$Sublineages %in% colnames(df_b)),"col"])

c = c + scale_fill_manual(values = sublin_col[which(sublin_col$Sublineages %in% colnames(df_c)),"col"])


g = arrangeGrob(c+ theme(legend.position = "none",panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("a. Khulna (2015-2018)"),
                b+ theme(legend.position = "none",panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("b. Assam (2002-2005)"),
                a + theme(legend.position = "none",panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("c. Nepal (2010)"), ncol = 3)
g_legend = arrangeGrob(c+ guides(fill ="none")+ theme(panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("a. Khulna (2015-2018)"),
                       b+ theme(legend.position = "none",panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("b. Assam (2002-2005)"),
                       a + theme(legend.position = "none",panel.border = element_rect(colour = "white", fill=NA, linewidth=5), title = element_text(size = 12)) + ggtitle("c. Nepal (2010)"), ncol = 3)
ggsave("Figure S13.png", g, height = 3, width = 9)
ggsave("Figure S13_legend.svg", g_legend, height = 3, width = 9)

#Legend
ggplot(sublin_col, aes(x = 1, y = 1, fill = Sublineages)) + geom_tile(color = "black") + guides(fill=guide_legend(ncol=8)) + scale_fill_manual(values = sublin_col$col) + geom_line(data = data.frame(x = 1:9, y = 1:9, Sublineages = "BD2.002"), aes(x = x, y = y, linewidth = y)) + theme_void() + labs(linewidth = "No. shared sublineages")
ggsave("Figure S13 legend.svg", height = 3, width = 9)




supp1 = metadata
supp2 = vc_clean$annotation[,-c(25:26)]
write.csv(supp1, file = "../Supplementary1.csv",row.names = F)
write.csv(supp2, file = "../Supplementary2.csv", row.names = F)


#### ICP1 prevalence and tree (Extended data Fig. 8) ####
metagenomics = read.table("SupplementaryDocument3.txt", header = T)

#ICP1 prevalence
metagenomics$Category = ""
metagenomics[which(metagenomics$Vibrio_cholerae_screen_sylph == "Cholera" & metagenomics$ICP1_screen_kraken_viraldb_percent >= 0.01), "Category"] = "Vibrio cholerae and ICP1"
metagenomics[which(metagenomics$Vibrio_cholerae_screen_sylph == "No cholera" & metagenomics$ICP1_screen_kraken_viraldb_percent >= 0.01), "Category"] = "ICP1"
metagenomics[which(metagenomics$Vibrio_cholerae_screen_sylph == "Cholera" & metagenomics$ICP1_screen_kraken_viraldb_percent < 0.01), "Category"] = "Vibrio cholerae"

frequency = merge(as.data.frame(table(metagenomics$Country)), as.data.frame(table(metagenomics[which(!is.na(metagenomics$Category)), c("Country", "Category")])), by.x = "Var1", by.y = "Country", all.x = T)
colnames(frequency) = c("Country", "Total genomes", "Category", "Freq")
frequency = frequency[which(frequency$Category != ""),]
frequency$Percent = 100*frequency$Freq/frequency$`Total genomes`
frequency$Label = paste0(frequency$Country, " (n=", frequency$`Total genomes`, ")")
frequency$Label = factor(frequency$Label, levels = rev(unique(sort(frequency$Label))))
frequency$Category = factor(frequency$Category, levels = c("ICP1", "Vibrio cholerae and ICP1", "Vibrio cholerae", ""))

ggplot(frequency, aes(x = Label, y = Freq, fill = Category)) + geom_bar(stat = "identity", color = "black", size = 0.2) + coord_flip() + theme_bw() + labs(x = "Country",y = "% total microbiomes") + geom_text(data = frequency[which(frequency$Freq > 0),],aes(label = Freq), position = position_stack(vjust = .5), size = 3) + theme(panel.grid = element_blank(), legend.position = "bottom") + scale_fill_manual(values = c("deeppink", "deeppink3", "grey")) + geom_hline(yintercept = 0, size =0.2)

#ICP1 phylogenetic tree
icp1_tree = read.tree("Trees/icp1.treefile")
icp1_tree = midpoint.root(icp1_tree)

tree_metadata = read.delim("SupplementaryDocument4.txt", row.names = 1)

heatmap = tree_metadata[,c("Country", "ctxB_type", "PLE", "SXT", "Source", "Antidefence", "Date_category", "Lineage"),]
colnames(heatmap) = c("Country", "ctxB", "PLE", "SXT", "Source", "Anti-defence", "Date",  "Lineage")
heatmap = heatmap[,c("Country", "Anti-defence", "Date", "Source", "ctxB", "PLE", "SXT",  "Lineage")]
heatmap[which(heatmap$PLE == "ND"),"PLE"] = "None"
heatmap[which(heatmap$SXT == "ND"),"SXT"] = "None"
for(i in 1:ncol(heatmap)){heatmap[which(!is.na(heatmap[,i])),i] = paste0(colnames(heatmap)[i], ": ", as.character(heatmap[which(!is.na(heatmap[,i])),i]))}
heatmap[which(heatmap$ctxB %in% c("ctxB: ND", "ctxB: None")),"ctxB"] = NA

gheatmap(ggtree(icp1_tree)+ geom_nodelab(), heatmap, color = NA, colnames_angle = 90, colnames_offset_y = -10, font.size = 4) + scale_fill_manual(values = c(
  "#761152", "#FA9FB5", "#C51B8A",#defence
  "#7BCCC4","#FFF2CC", "#0868AC", "#C73157", #countries
  "#E41A1C", "#FF7F00", #ctxB
  brewer.pal(3, "Blues")[c(2,3,1)], #years
  "#A8DDB5","#BDBDBD","#FA9FB5",#Lineage
  "grey90", "white", "#9E0142", #PLE
  brewer.pal(3, "Accent"), #Source
  "#FC8D62", "#8DA0CB", "#E78AC3", "grey90", "white"
),

) + guides(fill = F) 

ggsave("ICP1tree.svg", width = 14, height = 10)





#### Phenotypic resistance to antibiotics (Extended data Fig. 9) ####
sxt = as.data.frame(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$SXT_phenotype) & metadata$SXT_phenotype != ""),c("SXT","SXT_phenotype")]))
sxt$Ab = "SXT"
colnames(sxt)[2] = "Resistance"

tet = as.data.frame(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh"  & !is.na(metadata$TE_phenotype) & metadata$TE_phenotype != ""),c("SXT","TE_phenotype")]))
tet$Ab = "Tetracycline"
colnames(tet)[2] = "Resistance"

dox = as.data.frame(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$DO_phenotype) & metadata$DO_phenotype != ""),c("SXT","DO_phenotype")]))
dox$Ab = "Doxycycline"
colnames(dox)[2] = "Resistance"

ab = rbind(sxt, tet, dox)
table(ab$SXT)
ab = ab[which(ab$SXT %in% c("ICE-TET", "ICE-GEN", "None")),]
ab$SXT  = factor(ab$SXT, levels = c("ICE-TET", "ICE-GEN", "None"))
ab$Resistance = gsub("R", "Resistant", ab$Resistance) %>% gsub("S", "Susceptible", .) %>% gsub("I", "Intermediate", .) %>% factor(., levels = c("Resistant", "Intermediate", "Susceptible"))
ab$Ab = factor(ab$Ab, levels = c("SXT", "Tetracycline", "Doxycycline"))
for(i in unique(ab$SXT)){
  for(j in unique(ab$Ab)){
    ab[which(ab$Ab == j & ab$SXT == i),"Percent"] = 100*ab[which(ab$Ab == j & ab$SXT == i),"Freq"]/sum(ab[which(ab$Ab == j & ab$SXT == i),"Freq"])
  }}
ab=ab[-which(ab$Percent == 0),]

for(i in c("SXT_phenotype", "DO_phenotype", "TE_phenotype")){
  metadata[which(metadata[,i] == ""),i] = NA
}

ab$p = NA
for(comparison in c("ICE-GEN", "None")){
  if(length(unique(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$SXT_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT")])) > 1){
    ab[which(ab$SXT == comparison & ab$Ab == "SXT"),"p"][1] = fisher.test(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$SXT_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT","SXT_phenotype")]))$p.value}
  
  
  if(length(unique(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$TE_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT")])) > 1){ 
    ab[which(ab$SXT == comparison & ab$Ab == "Tetracycline"),"p"][1] = fisher.test(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$TE_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT","TE_phenotype")]))$p.value}
  
  if(length(unique(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$DO_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT")])) > 1){ 
    ab[which(ab$SXT == comparison & ab$Ab == "Doxycycline"),"p"][1] = fisher.test(table(metadata[which(metadata$Year > 2013 & metadata$Country == "Bangladesh" & !is.na(metadata$DO_phenotype) & metadata$SXT %in% c("ICE-TET", comparison)),c("SXT","DO_phenotype")]))$p.value
  }
}

ab$label_y = 106
ab[which(ab$SXT == "ICE-GEN"), "label_y"] = 114
ab$label_x = 1.5
ab[which(ab$SXT == "ICE-GEN"), "label_x"] = 2
ab$SXT = factor(ab$SXT, levels = c("ICE-TET", "None", "ICE-GEN"))
ab$Ab = gsub("SXT", "Trimethoprim/sulfamethoxazole", ab$Ab)

ggplot(ab, aes(x = SXT, y = Percent, fill = Resistance)) + geom_bar(stat = "identity", color = "black") + facet_wrap(~Ab) + 
  theme_bw() + 
  scale_fill_manual(values = colorRampPalette(c("#1B998B", "white"))(3)) + 
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text_phenotype = element_text(angle = 90)) + 
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.6), color = "black") + labs(y = "%", x = "SXT-ICE") + annotate(geom = "segment", x = 1, xend = 2, y = 102, yend = 102) + 
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.6), color = "black") + labs(y = "%", x = "SXT-ICE") + annotate(geom = "segment", x = 1, xend = 3, y = 110, yend = 110) +
  geom_text(aes(label = signif(p, 2), x= label_x, y = label_y))





#### BD1/BD2 subtrees (Supplementary Figure 2) ####

#Subset tree
bd1bd2 = tree_subset(vc_clean$tree, 5490, levels_back = 0)
bd1bd2_annotation = merge(as.data.frame(fortify(bd1bd2)), vc_clean$annotation[,c(1,10:ncol(vc_clean$annotation))], by = "label", all.x = T)
bd1bd2 = list(tree = bd1bd2, annotation = bd1bd2_annotation)
recent = max(as.Date(bd1bd2$annotation$date), na.rm = T)

#Color by sublineage
sublineagecolors = color_scheme[grep("Sublineage", color_scheme$value),]
sublineagecolors$value = gsub("Sublineage: ", "", sublineagecolors$value)
sublineagecolors = sublineagecolors[which(sublineagecolors$value %in% bd1bd2$annotation$Sublineage),]

#Create legend
bd1bd2$annotation$labels = ""
bd1bd2$annotation = bd1bd2$annotation[order(bd1bd2$annotation$y),]
#Only label middle of sublineages with more than 50 members
for(i in names(which(table(no_na((bd1bd2$annotation$Sublineage))) >= 50))){
  rows = bd1bd2$annotation[which(bd1bd2$annotation$Sublineage == i & bd1bd2$annotation$isTip == T),"y"]
  middle = rows[floor(length(rows)/2)]
  bd1bd2$annotation[which(bd1bd2$annotation$y == middle),"labels"] = i}
#Create bounding box which excludes outliers
boundingbox = data.frame(Sublineage = unique(bd1bd2$annotation$labels)[which(unique(bd1bd2$annotation$labels)!="")], ymin = NA, ymax = NA)
for(i in boundingbox$Sublineage){
  
  lower = quantile(bd1bd2$annotation[which(bd1bd2$annotation$Sublineage == i &  bd1bd2$annotation$isTip == T),"y"], 0.01)
  upper = quantile(bd1bd2$annotation[which(bd1bd2$annotation$Sublineage == i &  bd1bd2$annotation$isTip == T),"y"], 0.99)
  
  boundingbox[which(boundingbox$Sublineage == i),"ymin"] = bd1bd2$annotation[which(abs(bd1bd2$annotation$y - lower) == min(abs(bd1bd2$annotation$y - lower))),"y"][1]
  
  boundingbox[which(boundingbox$Sublineage == i),"ymax"] = bd1bd2$annotation[which(abs(bd1bd2$annotation$y - upper) == min(abs(bd1bd2$annotation$y - upper))),"y"][1]
  
}

ggtree(bd1bd2$tree, aes(color = Sublineage), mrsd=recent) %<+% bd1bd2$annotation + theme_tree2() + scale_color_manual(values = sublineagecolors$color)  + theme(legend.position = "none")
ggsave("SublineageTreeonly.svg", height = 10, width = 15)

ggplot(bd1bd2$annotation[which(bd1bd2$annotation$isTip == T),], aes(y = y, x= "1", fill = Sublineage)) + geom_tile() + scale_fill_manual(values = sublineagecolors$color) + geom_rect(data = boundingbox, aes(ymin = ymin, ymax= ymax, xmin = 0.5, xmax= 1.5, y = 1, x = 1), fill = NA, color = "black") + geom_text(aes(label = labels), color = "white") + annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 3552, fill = NA, color = "black") + theme_void() + theme(legend.position = "none") 
ggsave("SublineageTreelabels.svg", height = 10, width = 1)

#Color by MGE/gene profile
bd1bd2$annotation$Profile = ""
bd1bd2$annotation[grep("VC_0491-ins", (bd1bd2$annotation$VSP_II)),"Profile"] = "VC_0491-ins"
bd1bd2$annotation[grep("ΔVC_0491-VC_0498", (bd1bd2$annotation$VSP_II)),"Profile"] = paste(bd1bd2$annotation[grep("ΔVC_0491-VC_0498", (bd1bd2$annotation$VSP_II)),"Profile"] , "ΔVC_0491-VC_0498")
bd1bd2$annotation[which(bd1bd2$annotation$VC0258 %in% c("fragmented", "interrupted", "partial", "no")),"Profile"] = paste(bd1bd2$annotation[which(bd1bd2$annotation$VC0258 %in% c("fragmented", "interrupted", "partial", "no")),"Profile"], "wbeT-ins")
bd1bd2$annotation[which(bd1bd2$annotation$PLE != "None"),"Profile"] = paste(bd1bd2$annotation[which(bd1bd2$annotation$PLE != "None"),"Profile"], bd1bd2$annotation[which(bd1bd2$annotation$PLE != "None"),"PLE"])
bd1bd2$annotation[which(bd1bd2$annotation$SXT != "None"),"Profile"] = paste(bd1bd2$annotation[which(bd1bd2$annotation$SXT != "None"),"Profile"], bd1bd2$annotation[which(bd1bd2$annotation$SXT != "None"),"SXT"])
bd1bd2$annotation$Profile = trimws(bd1bd2$annotation$Profile)
sort(table(bd1bd2$annotation$Profile), decreasing = T)

levels = c("PLE1 ICE-TET", "VC_0491-ins PLE1 ICE-TET", "ΔVC_0491-VC_0498 PLE1 ICE-TET","ΔVC_0491-VC_0498 wbeT-ins PLE1 ICE-TET", "ΔVC_0491-VC_0498 wbeT-ins ICE-TET", "ΔVC_0491-VC_0498 wbeT-ins", "ICE-GEN", "PLE11 ICE-GEN", "wbeT-ins PLE11 ICE-GEN")
bd1bd2$annotation[which(!bd1bd2$annotation$Profile %in% levels),"Profile"] = "Other"
bd1bd2$annotation$Profile = factor(bd1bd2$annotation$Profile, levels = c(levels, "Other"))

ggtree(bd1bd2$tree, aes(color = Profile), mrsd=recent) %<+% bd1bd2$annotation  + scale_color_manual(values = c(viridis(6), magma(6)[3:5], "grey", "white")) + theme_tree2() + guides() + theme(legend.position = "none")

ggsave("MGETreeonly.svg", height = 10, width = 15)




#### Subsampled transmission events within Bangladesh (Supplementary Figure 3) ####

#Calculate how many need to sub-sample for each region to get even representation
representation = read.csv("representation.csv")
representation$variable = make.names(as.character(representation$variable))

#Remove Mymenshing and combine Rajshahi/Rangpur and Khulna/Barisal so that when subsample, the number of samples isn't too low
representation = representation[which(representation$Division != "Mymenshing"),]
representation$Division = as.character(representation$Division)
representation[which(representation$Division  %in% c("Rajshahi", "Rangpur")),"Division"] = "Northern Bengal"
representation[which(representation$Division  %in% c("Khulna", "Barisal")),"Division"] = "Southern Bengal"
representation = aggregate(representation$value, by = list(Division = representation$Division, Year = representation$Year, variable = representation$variable, label = representation$label), FUN = "sum")

for(i in 2014:2018){
  
  df = dcast(unique(representation[which(representation$Year == i & (!representation$Division %in% c("Mymenshing"))),c("Division", "variable","x")]), Division~variable, id.vars = "Division", value.var = "x")
  
  # diarrhoea:cholera ratio
  diarrhoea_cholera_ratio = sum(df$AWD)/sum(df$Cholera)
  
  # expected diarrhoea cases proportional to population
  df$expected_diarrhoea <- round(df$Population/sum(df$Population)*sum(df$AWD))
  
  # expected cholera cases adjusted to maintain ratio
  df$expected_cholera = round(df$expected_diarrhoea/diarrhoea_cholera_ratio)
  
  # cap expected values at observed
  df$diarrhoea_subsample = pmin(df$AWD, df$expected_diarrhoea)
  df$cholera_subsample <- pmin(df$Cholera, df$expected_cholera)
  
  # adjust as number included in final analysis may be lower than number of cholera cases
  df$Year = i
  df$Final_subsample_number = floor(df$cholera_subsample/max(df$cholera_subsample/df$Final.analysis))
  
  #Combine into one dataframe
  if(i == 2014){adjusted = df}else{adjusted = rbind(adjusted, df)}
}

#Import results of subsampling
subsample = read.table("BD1_BD2_subsampling.txt", header = T)
subsample[which(subsample$Division %in% c("Rajshahi", "Rangpur")),"Division"] = "Northern Bengal\n(Rajshahi/Rangpur)"
subsample[which(subsample$Division %in% c("Khulna", "Barisal")),"Division"] = "Southern Bengal\n(Khulna/Barisal)"
subsample$Lineage = gsub("_.*", "", subsample$Sample)
subsample$SampleNumber = gsub(".*_sample", "", subsample$Sample)
subsample[which(subsample$Lineage == "BD1"),"Lineage"] = "sBD1"

#Find out parent location of each node
for(i in 1:nrow(subsample)){
  subsample[i,"Parent_Division"] = subsample[which(subsample$node == subsample[i,"parent"] & subsample$Sample ==  subsample[i,"Sample"]),"Division"]}

#Subset to transitions between regions
transitions = as.data.frame(table(subsample[,c("Sample","Division", "Parent_Division")]))
transitions = transitions[-which(transitions$Division == transitions$Parent_Division),]

#Get centroids of each region
division = bangladesh::get_map("division")
centroid = st_drop_geometry(st_centroid(division) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1],lat = sf::st_coordinates(.)[,2]))
centroid$Division = as.character(centroid$Division)
centroid[which(centroid$Division %in% c("Rajshahi","Rangpur")),"Division"] = "Northern Bengal\n(Rajshahi/Rangpur)"
centroid[which(centroid$Division == "Northern Bengal\n(Rajshahi/Rangpur)"),"lon"] = median(centroid[which(centroid$Division == "Northern Bengal\n(Rajshahi/Rangpur)"),"lon"])
centroid[which(centroid$Division == "Northern Bengal\n(Rajshahi/Rangpur)"),"lat"] = median(centroid[which(centroid$Division == "Northern Bengal\n(Rajshahi/Rangpur)"),"lat"])

centroid[which(centroid$Division %in% c("Khulna","Barisal")),"Division"] = "Southern Bengal\n(Khulna/Barisal)"
centroid[which(centroid$Division == "Southern Bengal\n(Khulna/Barisal)"),"lon"] = median(centroid[which(centroid$Division == "Southern Bengal\n(Khulna/Barisal)"),"lon"])
centroid[which(centroid$Division == "Southern Bengal\n(Khulna/Barisal)"),"lat"] = median(centroid[which(centroid$Division == "Southern Bengal\n(Khulna/Barisal)"),"lat"])

centroid = unique(centroid[,c("Division", "lon", "lat")])

transitions = merge(transitions, centroid, by = "Division", all.x = T)
transitions = merge(transitions, centroid, by.x = "Parent_Division", by.y = "Division", all.x = T)
transitions$Transition = paste(transitions$Parent_Division, transitions$Division, sep = "-> ")

#Plot median for lineage
transitions[which(transitions$lin)]

ggplot(transitions[grep("BD2", transitions$Sample),], aes(fill = Division)) + 
  geom_sf(data = division) + 
  scale_fill_manual(values = brewer.pal(9, "Greys")[c(2,3,4,2,1,5,5,6, 1, 1)]) +
  geom_curve(aes(yend = lat.y, xend = lon.y, y = lat.x, x = lon.x, size = Freq, alpha = Freq, color = Freq), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), curvature = 0.2) +
  scale_alpha(range =c(0.1, 0.3)) + 
  scale_size_continuous(range= c(0.2,3)) + 
  theme_void() + 
  scale_color_gradientn(colors = brewer.pal(9, "GnBu")) + 
  guides(alpha = F) 

ggplot(trans[grep("BD1", trans$Sample),], aes(y = Parent_Division, x = Division, fill = Freq)) + facet_wrap(~Sample, ncol = 5) + geom_tile() + scale_fill_gradientn(colors = c(brewer.pal(9, "RdPu"), "white"))
ggplot(trans[grep("BD2", trans$Sample),], aes(y = Parent_Division, x = Division, fill = Freq)) + facet_wrap(~Sample, ncol = 5) + geom_tile() + scale_fill_gradientn(colors = c(brewer.pal(9, "RdPu"), "white"))

ggplot(as.data.frame(table(subsample[which(subsample$Division != subsample$Parent_Division),c("SampleNumber","Parent_Division", "BAPS1")])), aes(x = Parent_Division, y= Freq)) + geom_boxplot() + geom_jitter(width = 0.2) + facet_wrap(~BAPS1) + theme_bw() + labs(x = "Source", y = "Number of inter-regional transmission events") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Origin of BD1
subsample$Date = as.numeric(subsample$Date)
subsample = subsample[order(subsample$Date),]
bd1 = subsample[which(subsample$BAPS1 == "BD1"),]
bd1 = bd1[which(!is.na(bd1$Division)),]
bd1[!duplicated(bd1$Sample),]

sort(table(trans$Sample, paste(trans$Parent_Division, trans$Division)))
ggplot(trans, aes(x = Sample))


#### Nucleotide diversity and tajima's D over time (Supplementary Figure 4) ####
tajimasd = read.delim("tajimasd.txt")

#Diversity
ggplot(tajimasd, aes(x = Year, y = pi)) +
  geom_jitter(aes(color = Country, fill = Country), shape = 21, size = 2.5, height = 0, width = 0.1) + 
  scale_color_manual(values = darken(c(brewer.pal(9, "GnBu")[c(5,8)]), 0.2)) + 
  scale_fill_manual(values = c(brewer.pal(9, "GnBu")[c(5,8)]))+
  guides(color = F, fill = F)+
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c(brewer.pal(9, "GnBu")[c(5,8)]))+
  geom_smooth(aes(color = Country, fill = Country), span = 0.1) + 
  theme_bw() + geom_hline(yintercept = 0) + 
  theme(panel.grid = element_blank()) +
  labs(y = "Nucleotide diversity") 

#Tajima's D
ggplot(tajimasd, aes(x = Year, y = D)) +
  geom_jitter(aes(color = Country, fill = Country), shape = 21, size = 2.5, height = 0, width = 0.1) + 
  scale_color_manual(values = darken(c(brewer.pal(9, "GnBu")[c(5,8)]), 0.2)) + 
  scale_fill_manual(values = c(brewer.pal(9, "GnBu")[c(5,8)]))+
  guides(color = F, fill = F)+
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c(brewer.pal(9, "GnBu")[c(5,8)]))+
  geom_smooth(aes(color = Country, fill = Country), span = 0.1) + 
  theme_bw() + geom_hline(yintercept = 0) + 
  theme(panel.grid = element_blank()) +
  labs(y = "Tajima's D")

#### MGE dynamics in different countries (Supplementary Figure 5) ####
global_subregions = read.delim("global_subregions.txt")
MGE = merge(metadata, global_subregions, by = "Country", all.x = T)
MGE$MGE_Replacement = ""
for(i in c("VC_0490.assembled", "VC_0491.assembled", "VC_0492.assembled")){
  MGE[which(MGE[,i] == "yes"),"MGE_Replacement"] = paste(MGE[which(MGE[,i] == "yes"),"MGE_Replacement"],  gsub("_", "", gsub(".assembled","", i)))}
MGE$Year = as.numeric(as.character(MGE$Year))
MGE$Location = MGE$Subregion
MGE[which(MGE$Country %in% c("Bangladesh", "India", "Democratic Republic of the Congo", "Mexico")),"Location"] = MGE[which(MGE$Country %in% c("Bangladesh", "India", "Democratic Republic of the Congo", "Mexico")),"Country"]
other = names(which(table(MGE$MGE_Replacement ) < 20))
MGE[which(MGE$MGE_Replacement %in% other),"MGE_Replacement"] = "Other"
MGE = MGE[-which(MGE$Location %in% c("Europe", "Oceania", "South America", "Unknown", "Eastern Asia", "Southeast Asia")),]
MGE$Location = gsub("Democratic Republic of the Congo", "DRC", MGE$Location)
MGE$PLE = factor(MGE$PLE, levels = c("None", paste0("PLE", 1:11)))
MGE$MGE_Replacement = gsub(" VC0490 VC0491 VC0492", "DdmABC", MGE$MGE_Replacement) %>% gsub(" VC0491 VC0492" , "DdmAB", .) %>% gsub(" VC0492"  , "DdmA", .)
MGE$MGE_Replacement[which(MGE$MGE_Replacement == "")] = "None"
MGE[grep("yes", MGE$VC_0258.assembled),"Ogawa"] = "Present"
MGE[which(MGE$VC_0258.assembled %in% c("partial", "fragmented", "disrupted")),"Ogawa"] = "Partial/disrupted/fragmented"
MGE[which(MGE$VC_0258.assembled %in% c("interrupted")),"Ogawa"] = "Interrupted"
MGE[which(MGE$VC_0258.assembled == "no"),"Ogawa"] = "Not present"


ggplot(MGE, aes(x = Year, point_color = PLE, fill =PLE, color = MGE_Replacement)) + 
  facet_wrap(~Location, scale = "free_y", nrow = 2)+
  geom_bar(stat = "count", color = "black") + 
  theme_bw() + xlim(2003,2023) + scale_fill_manual(values = c("white", brewer.pal(10, "Spectral")[c(1:4,9:10)])) + labs(y = "Frequency") + ggtitle("PLEs") + theme(legend.title = element_blank(), legend.position = "bottom") + guides(fill=guide_legend(nrow =1 ))

ggplot(MGE, aes(x = Year, point_color = MGE_Replacement, fill =MGE_Replacement, color = MGE_Replacement)) + 
  facet_wrap(~Location, scale = "free_y", nrow = 2)+
  geom_bar(stat = "count", color = "black") + 
  theme_bw() + xlim(2003,2023) + scale_fill_manual(values = c(brewer.pal(3, "GnBu"),"white", "grey")) + ggtitle("DdmABC") + theme(legend.title = element_blank(), legend.position = "bottom")+ labs(y = "Frequency")

ggplot(MGE, aes(x = Year, point_color = Ogawa, fill =Ogawa, color = MGE_Replacement)) + facet_wrap(~Location, scale = "free_y", nrow = 2)+
  geom_bar(stat = "count", color = "black") + 
  theme_bw() + xlim(2003,2023) + scale_fill_manual(values = c( c("grey95","white", "grey85", brewer.pal(8, "Dark2")[1])))+ ggtitle("wbeT (VC0258)") + theme(legend.title = element_blank(), legend.position = "bottom")+ labs(y = "Frequency")


#### Inferred international transmission events following sub-sampling (Supplementary Figure 6) ####
transmission_all = read.delim("Global_transmission_subsampling.txt")
ggplot(as.data.frame(table(transmission_all[,c("Category","Decade", "Sample")])), aes(x = Category, y = Freq, fill = Category, color = Category)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + facet_wrap(~Decade) + 
  geom_jitter(shape = 21, size = 2.5) + 
  scale_fill_manual(values =  c(brewer.pal(9, "GnBu")[c(5,8)], "grey"  )) + 
  scale_color_manual(values =  darken(c(brewer.pal(9, "GnBu")[c(5,8)], "grey"  ), 0.5)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Source Country", y ="Inferred transmission events") +
  guides(fill = F, color = F)

#### Proportion of the population, acute watery diarrhoea cases, cholera cases, sequenced samples, and samples included in the final analysis, from each Division (Supplementary Figure 7) ####
representation = read.table("representation.txt", header = T)

ggplot(representation, aes(x = label, y = value, fill = Division)) + geom_bar(stat = "identity",position = "fill") + facet_wrap(~Year, ncol = 5)  + theme_bw() + scale_fill_brewer(palette = "Set2") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Stage of analysis", y ="Proportion")



