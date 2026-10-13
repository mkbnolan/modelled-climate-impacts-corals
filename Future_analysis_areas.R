## Code authored by Megan K. B. Nolan ##

## Associated to the manuscript in review: Projected climate change impacts on deep-sea corals of the Red Sea (in review) ##

## Analysis for present vs future conditions and suitable habitat from ENMs ##

## Calculating results from rasters
  #1. Set up ----
setwd("Z:/nolanmc/ENMs/R_BRTs")
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(rstatix)

  #2. Create groups and load rasters ----
speciesname <- list("Acanthogorgia","Bathypathes", "Madracis", "Stylaster")
thresholds <- data.frame(Acanthogorgia = 0.237, Bathypathes = 0.195, Madracis = 0.391, Stylaster = 0.232)
scenario <- list("Present", "SSP119", "SSP126", "SSP245", "SSP370", "SSP585")

  #Calculate area of Red Sea (below 0 m)
depth_rs <- rast("EnvAll/Cropped/depth.tif.asc")
crs(depth_rs) <- "EPSG:4326"
rs_binary <- ifel(depth_rs < 0, 1, 0) 
area_rs <- sum(values(rs_binary) == 1, na.rm = TRUE) * 400 * 400 / 1000000

terra::plot(rs_binary, col = map.pal("viridis"))  


  #3. Create empty lists and dataframes ----
raster_results_all <- NULL
depths_all <- NULL
plot_list <- list()
plot_list_temp <- list()
plot_list_all <- list()
plot_list_temp_all <- list()
all_areas <- data.frame(species = character(), scenario = character(), threshold = numeric(), run=numeric(), 
                        total_area=numeric(), area_OG = numeric(), overlap_area=numeric())
lat_change <- data.frame(species= character(), scenario = character(), latitude=numeric())

  #4. Start of loop to calculate areas and plot histograms ----
for (a in speciesname){

  #Create dataframes to store results
raster_results <- data.frame(species = character(), scenario = character(), threshold = numeric(), total_area=numeric(), 
                            area_OG = numeric(), overlap_area=numeric())
depth_results <- data.frame(species = character(), scenario = character(), max_depth = numeric(),min_depth = numeric(),
                            mean_depth = numeric())
raster_list <- list()

  #start loop for scenarios
for (i in scenario){
  #Define threshold to create binary raster
  threshold = thresholds[[a]] 
  #Load raster for species and scenario, convert to binary by threshold
  rast_original <- rast(paste0(i,"/Plots/", a,"/mean", a,".tif"))
  crs(rast_original) <- "EPSG:4326" #ensure correct CRS
  binary_raster <- ifel(rast_original > threshold, 1, 0)

  #Extract area for presence
  area_1 <- sum(values(binary_raster) == 1, na.rm = TRUE) *400 *400 / 1000000

  #Extract area that falls within original (overlap)
  OGrast <- rast(paste0("Present/Plots/", a,"/mean",a,".tif"))
  crs(OGrast) <- "EPSG:4326"
  OGbinary <- ifel(OGrast > threshold, 1, 0)  #make binary with same threshold
  area_2 <- sum(values(OGbinary) == 1, na.rm = TRUE) * 400 * 400 / 1000000 #convert to km^2

  #Find overlap for where OGbinary and binary_raster are both 1
  overlap_raster <- binary_raster * OGbinary
  overlap_area <- sum(values(overlap_raster) == 1, na.rm = TRUE) * 400 * 400 / 1000000 #convert to km^2

    #Add all results to dataframe
  raster_results <- rbind(raster_results, data.frame(species= a, scenario = i, threshold = threshold, total_area=area_1, 
                                      area_OG = area_2, overlap_area = overlap_area))
  
    for (b in 1:10){
      accuracy_df <- read.csv(paste0("Present/Plots/", a,"/accuracy_",a,".csv")) 
      threshold_b <- accuracy_df$threshold[[b]]
      
      rast_original_b <- rast(paste0(i,"/Plots/", a,"/map", a, b,".tif"))
      crs(rast_original_b) <- "EPSG:4326" #ensure correct CRS
      binary_raster_b <- ifel(rast_original_b > threshold_b, 1, 0)
      
      #Extract area for presence
      area_1b <- sum(values(binary_raster_b) == 1, na.rm = TRUE) *400 *400 / 1000000
      
      #Find overlap for where OGbinary and binary_raster_b are both 1 (use binary from loop above)
      overlap_raster_b <- binary_raster_b * OGbinary
      overlap_area_b <- sum(values(overlap_raster_b) == 1, na.rm = TRUE) * 400 * 400 / 1000000 #convert to km^2  
    
      all_areas <- rbind(all_areas, data.frame(species=a, scenario=i, threshold=threshold_b, run=b, total_area=area_1b,
                                               overlap_area=overlap_area_b))
        
    }
  
  #Write rasters comparing present to future scenario
  gained <- ifel(rast_original > threshold, 1, 0)  #gained: score 1
  lost <- ifel(OGrast > threshold, 2, 0)        #lost: score 2
  retained <- gained + lost           #retained: score 3
  writeRaster(retained, paste0("Future_Analysis/", a, i, "gain_retain.tif"), overwrite=TRUE)
  
  #Find depth range for each raster
  depth_range <- values(depth_rs)[values(binary_raster) == 1]
  depth_df <- data.frame(depth = depth_range)
  range(depth_df$depth, na.rm = TRUE)
  mean_depth <- round(mean(depth_df$depth, na.rm = TRUE),0)*-1
  
  depth_results <- rbind(depth_results,data.frame(species=a, scenario =i, max_depth = (min(depth_df$depth, na.rm=T)*-1),
             min_depth = (max(depth_df$depth, na.rm=T)*-1), mean_depth = mean_depth))

  hist <- ggplot(depth_df, aes(x = (depth*-1))) +
    geom_histogram(binwidth = 10, fill = "steelblue", color = NA) +
    labs(title = paste0(a,": ",i),
         x = "Depth (m)",
         y = "Frequency") +
    annotate("text", x=-Inf, y=Inf, label = paste0(mean_depth, " m"), hjust= -0.2, vjust = 1.5, size=6) +
    theme_classic()
    
  plot_list[[i]] = hist
  
  #Plot as a function of temperature
  if (i == "Present"){
    temp <- rast("EnvAll/Cropped/tempMean.tif.asc")
  } else {
    temp <- rast(paste0(i,"/Unused/tempMean.asc"))
  }
  
  temp_values <- values(temp)[values(binary_raster) == 1]
  temp_df <- data.frame(temp = temp_values)
  mean_temp <- round(mean(temp_df$temp, na.rm = TRUE),2)
  
  hist_temp <- ggplot(temp_df, aes(x = temp)) +
    geom_histogram(binwidth = 0.2, fill = "steelblue", color = NA) +
    labs(title = paste0(a,": ",i),
         x = expression("Temperature ("^o*"C)"),
         y = "Frequency") +
    annotate("text", x=-Inf, y=Inf, label = paste0(mean_temp, " °C"), hjust= -0.2, vjust = 1.5, size=6) +
    theme_classic()
  
  plot_list_temp[[i]] = hist_temp
  
    #Find average latitude of species for all scenarios
  match_r <- binary_raster == 1
  match_points <- as.points(match_r, na.rm = TRUE)
  match_points <- match_points[values(match_points)[,1] == 1, ]
  
  coords <- crds(match_points, df=T)
  mean_latitude <- mean(coords[, "y"])
  lat_change <- rbind(lat_change, data.frame(species=a, scenario=i, latitude=mean_latitude))
  
  }

raster_results$OG_percent <- 100
raster_results$percent_total <- raster_results$total_area / raster_results$area_OG * 100
raster_results$percent_spatial <- raster_results$overlap_area / raster_results$area_OG * 100

raster_results_all <- rbind(raster_results_all, raster_results)
write.csv(raster_results, paste0("Future_Analysis/",a,"areas.csv"))

depths_all <- rbind(depths_all, depth_results)
write.csv(depths_all, "Future_Analysis/depths.csv")

#Combined plots for depth histograms
plot_list_all <- c(plot_list_all, plot_list)

custom_labels <- c(paste0(letters[1:24], ")"))
combined_plot <- cowplot::plot_grid(plotlist = plot_list_all, ncol=6, nrow=4, labels=custom_labels, align="hv")  
plot(combined_plot)

ggsave(
  filename = paste0("Future_Analysis/histograms_depth.pdf"), 
  plot = combined_plot,   
  width = 20, height = 10)
  
#Combined plots for temperature histograms
plot_list_temp_all <- c(plot_list_temp_all, plot_list_temp)

combined_plot_temp <- cowplot::plot_grid(plotlist = plot_list_temp_all, ncol=6, nrow=4, labels=custom_labels, align="hv")  
plot(combined_plot_temp)

ggsave(
  filename = paste0("Future_Analysis/histograms_temp.pdf"), 
  plot = combined_plot_temp,   
  width = 20, height = 10)

}

  #5. Further analysis ----
all_areas_summary <- all_areas %>%
  group_by(species, scenario) %>%
  summarise(mean_absolute=mean(total_area),
            sd_absolute=sd(total_area),
            se_absolute=sd_absolute/sqrt(length(total_area)),
            mean_overlap=mean(overlap_area),
            sd_overlap=sd(overlap_area),
            se_overlap=sd_overlap/sqrt(length(overlap_area))
            )


#Reshape df to long
areas <- raster_results_all[,c(1:2,5:6)] %>%
  pivot_longer(cols = c(area_OG, overlap_area),
               names_to = "time",
               values_to = "area") %>%
  mutate(time = recode(time, "area_OG" = "Present", "overlap_area" = "Future"),
         time = factor(time, levels = c("Present", "Future")))

se_areas <- all_areas_summary[,c(1,2,8)]


areas_present <- subset(areas, time =="Present")
areas_future <- subset(areas, time == "Future")
se_areas_present <- subset(se_areas, scenario == "Present")

areas_future <- left_join(areas_future, se_areas, by=c("species", "scenario"))
areas_present <- left_join(areas_present, se_areas_present[,c(1,3)], by="species")

areas <- rbind(areas_present, areas_future)

areas <- subset(areas, scenario != "Present")

palette <- c("forestgreen", "gold2", "darkorange2", "red3", "red4")

overlap <- ggplot(areas, aes(x = time, y = area, group = scenario, color = scenario)) +
  geom_line(size = 1.5) +  # Line connecting Current to Future
  geom_point(size = 1) + # Add points for Current and Future
  geom_errorbar(aes(ymin = area - se_overlap, ymax = area + se_overlap), width = 0.2)+
  facet_wrap(~ species, scales = "free_y", nrow=1, labeller = labeller(
    species = c("Acanthogorgia" = "Acanthogorgia spp.", "Bathypathes" = "Bathypathes thermophila", 
                "Madracis" = "Madracis interjecta", "Stylaster" = "Stylaster tritoni"))) +  # Separate panels for species
  labs(x = "Time Period",
       y = expression("Area (km"^2*")"%+-% " SE"),
       color = "Scenario") +
  scale_color_manual(values = palette)+
  theme_classic()+
  theme(strip.text = element_text(size = 12, face="italic"), strip.background = element_blank())

ggsave("Future_Analysis/overlap_area.pdf", plot=overlap, height=3, width=12)


#Reshape df to long
areas_absolute <- raster_results_all[,c(1:2,5,4)] %>%
  pivot_longer(cols = c(area_OG, total_area),
               names_to = "time",
               values_to = "area") %>%
  mutate(time = recode(time, "area_OG" = "Present", "total_area" = "Future"),
         time = factor(time, levels = c("Present", "Future")))

se_areas_abs <- all_areas_summary[,c(1,2,5)]

#Find significance
normality_total <- all_areas %>%
  group_by(species, scenario) %>%
  shapiro_test(total_area)
normality_overlap <- all_areas %>%
  group_by(species, scenario) %>%
  shapiro_test(overlap_area)
## Some are normally distributed, but most are not so will use wilcoxon rank-sum (Mann Whitney U)

acanthogorgia_all <- subset(all_areas, species == "Acanthogorgia")
results_acan_total <- acanthogorgia_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(total_area ~ scenario, ref.group = "Present")
results_acan_overlap <- acanthogorgia_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(overlap_area ~ scenario, ref.group = "Present")

styl_all <- subset(all_areas, species == "Stylaster")
results_styl_total <- styl_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(total_area ~ scenario, ref.group = "Present")
results_styl_overlap <- styl_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(overlap_area ~ scenario, ref.group = "Present")

mad_all <- subset(all_areas, species == "Madracis")
results_mad_total <- mad_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(total_area ~ scenario, ref.group = "Present")
results_mad_overlap <- mad_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(overlap_area ~ scenario, ref.group = "Present")

bathy_all <- subset(all_areas, species == "Bathypathes")
results_bathy_total <- bathy_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(total_area ~ scenario, ref.group = "Present")
results_bathy_overlap <- bathy_all %>%
  mutate(scenario = as.factor(scenario)) %>%
  wilcox_test(overlap_area ~ scenario, ref.group = "Present")


areas_present <- subset(areas_absolute, time =="Present")
areas_future <- subset(areas_absolute, time == "Future")
se_areas_abs_present <- subset(se_areas_abs, scenario == "Present")

areas_future <- left_join(areas_future, se_areas_abs, by=c("species", "scenario"))
areas_present <- left_join(areas_present, se_areas_abs_present[,c(1,3)], by="species")

areas_absolute <- rbind(areas_present, areas_future)

areas_absolute <- subset(areas_absolute, scenario != "Present")

absolute <- ggplot(areas_absolute, aes(x = time, y = area, group = scenario, color = scenario)) +
  geom_line(size = 1.5) +  # Line connecting Current to Future
  geom_point(size = 1) + # Add points for Current and Future
  geom_errorbar(aes(ymin = area - se_absolute, ymax = area + se_absolute), width = 0.2)+
  facet_wrap(~ species, scales = "free_y", nrow=1, labeller = labeller(
    species = c("Acanthogorgia" = "", "Bathypathes" = "", 
                "Madracis" = "", "Stylaster" = ""))) +  # Separate panels for species
  labs(x = "Time Period",
       y = expression("Area (km"^2*")"%+-% " SE"),
       color = "Scenario") +
  scale_color_manual(values = palette)+
  theme_classic()+
  theme(strip.text = element_text(size = 12, face="italic"), strip.background = element_blank())

ggsave("Future_Analysis/absolute_area.pdf", plot=absolute, height=3, width=12)

#Standard deviaiton and error of AUC and TSS
accuracy_se <- data.frame(species = character(), mean_auc = numeric(), sd_auc=numeric(), 
                          se_auc=numeric(), mean_tss = numeric(), sd_tss=numeric(), se_tss=numeric(),
                          mean_thresh = numeric(), sd_thresh=numeric(), se_thresh=numeric())

for (a in speciesname) {
accuracy_df <- read.csv(paste0("Present/Plots/", a,"/accuracy_",a,".csv"))[1:10,] 
accuracy_summary <- accuracy_df %>%
  summarise(species=a,
            mean_auc=mean(auc),
            sd_auc=sd(auc),
            se_auc=sd_auc/sqrt(length(auc)),
            mean_tss=mean(tss),
            sd_tss=sd(tss),
            se_tss=sd_tss/sqrt(length(tss)),
            mean_thresh=mean(threshold),
            sd_thresh=sd(threshold),
            se_thresh=sd_thresh/sqrt(length(threshold)))

accuracy_se <- rbind(accuracy_se, accuracy_summary)
}

write.csv(accuracy_se, "Z:/nolanmc/Papers/ENMs/Analysis/Accuracy_se.csv")

##### Count species in each cell
biodiv <- data.frame(scenario= character(), area1=numeric(), area2=numeric(), area3=numeric(), area4=numeric())

for (i in scenario){
rast_original <- rast(paste0(i,"/Plots/Acanthogorgia/meanAcanthogorgia.tif"))
crs(rast_original) <- "EPSG:4326" #ensure correct CRS
acan_raster <- ifel(rast_original > 0.237, 1, 0)

rast_original <- rast(paste0(i,"/Plots/Madracis/meanMadracis.tif"))
crs(rast_original) <- "EPSG:4326" #ensure correct CRS
madracis_raster <- ifel(rast_original > 0.391, 1, 0)

rast_original <- rast(paste0(i,"/Plots/Stylaster/meanStylaster.tif"))
crs(rast_original) <- "EPSG:4326" #ensure correct CRS
styl_raster <- ifel(rast_original > 0.232, 1, 0)

rast_original <- rast(paste0(i,"/Plots/Bathypathes/meanBathypathes.tif"))
crs(rast_original) <- "EPSG:4326" #ensure correct CRS
bathy_raster <- ifel(rast_original > 0.195, 1, 0)

total <- acan_raster + madracis_raster + styl_raster + bathy_raster
writeRaster(total, paste0("Future_Analysis/numberSpecies_",i,".tif"), overwrite = TRUE)

area4 <- sum(values(total) == 4, na.rm = TRUE) * 400 * 400 / 1000000 #convert to km^2
area3 <- sum(values(total) == 3, na.rm = TRUE) * 400 * 400 / 1000000
area2 <- sum(values(total) == 2, na.rm = TRUE) * 400 * 400 / 1000000
area1 <- sum(values(total) == 1, na.rm = TRUE) * 400 * 400 / 1000000

values <- as.data.frame(cbind(i, area1, area2, area3, area4)) 
biodiv <- rbind(biodiv, values)

}

biodiv_long <- biodiv %>%
  pivot_longer(cols = c(area1, area2, area3, area4),
        names_to = "NumberSpecies", values_to = "area") 
biodiv_long$area <- as.numeric(biodiv_long$area)

palette6 <- c("darkgrey","forestgreen", "gold2", "darkorange2", "red3", "red4")


biodiv_plot <- ggplot(biodiv_long, aes(x = NumberSpecies, y=area, fill=i))+
  geom_bar(stat="identity", position="dodge", width=0.8) +  # Line connecting Current to Future
  labs(x = "Number of Species", y = expression("Area (km"^2*")"), fill = "Scenario") +
  scale_fill_manual(values = palette6)+
  scale_x_discrete(labels= c("1", "2", "3", "4"))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  theme(legend.position=c(1,1), legend.justification=c(1,1))
ggsave(plot=biodiv_plot, "Figures/biodiversity_hotspots.pdf", width=4, height=3)


lats <- ggplot(lat_change, aes(x=scenario, y=latitude, group=species, color=species))+
  geom_point( size=1.5)+
  geom_line(size=1)+
  theme_classic()+
  labs(x="Scenario", y="Average Latitude", color= "Species")+
  scale_x_discrete(labels=c("Present", "SSP1-1.9", "SSP1-2.6","SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))+
  scale_color_discrete(labels=c(expression(italic("Acanthogorgia") ~ "spp."), 
                                expression(italic("Bathypathes thermophila")), 
                                expression(italic("Madracis interjecta")),
                                expression(italic("Stylaster tritoni"))))

ggsave(plot=lats, "Future_Analysis/latitudes.pdf")                 

  #6. Venn diagram analysis ----
  ## Alternative code also provided in vennDiagram.R- this one was ultimately not used
acan_rast <- rast("Present/Plots/Acanthogorgia/meanAcanthogorgia.tif")
crs(acan_rast) <- "EPSG:4326" #ensure correct CRS
acan_present <- ifel(acan_rast > 0.237, 1, 0)
acan_pres_only <- mask(acan_present, acan_present == 1, maskvalues=F)
acan_poly <- as.polygons(acan_pres_only, na.rm=T)
a_acan <- expanse(acan_poly)/1e+6 

bathy_rast <- rast("Present/Plots/Bathypathes/meanBathypathes.tif")
crs(bathy_rast) <- "EPSG:4326" #ensure correct CRS
bathy_present <- ifel(bathy_rast > 0.195, 1, 0)
bathy_pres_only <- mask(bathy_present, bathy_present == 1, maskvalues=F)
bathy_poly <- as.polygons(bathy_pres_only, na.rm=T)
a_bathy <- expanse(bathy_poly)/1e+6

madr_rast <- rast("Present/Plots/Madracis/meanMadracis.tif")
crs(madr_rast) <- "EPSG:4326" #ensure correct CRS
madr_present <- ifel(madr_rast > 0.391, 1, 0)
madr_pres_only <- mask(madr_present, madr_present ==1, maskvalues =F)
madr_poly <- as.polygons(madr_pres_only, na.rm=T)
a_madr <- expanse(madr_poly)/1e+6

styl_rast <- rast("Present/Plots/Stylaster/meanStylaster.tif")
crs(styl_rast) <- "EPSG:4326" #ensure correct CRS
styl_present <- ifel(styl_rast > 0.232, 1, 0)
styl_pres_only <- mask(styl_present, styl_present ==1, maskvalues=F)
styl_poly <- as.polygons(styl_pres_only, na.rm=T)
a_styl <- expanse(styl_poly)/1e+6


# Pairwise intersections
acan_bathy <- intersect(acan_poly, bathy_poly)
a_acan_bathy <- expanse(acan_bathy)/1e+6
acan_madr <- intersect(acan_poly, madr_poly)
a_acan_madr <- expanse(acan_madr)/1e+6
acan_styl <- intersect(acan_poly, styl_poly)
a_acan_styl <- expanse(acan_styl)/1e+6
bathy_madr <- intersect(bathy_poly, madr_poly)
a_bathy_madr <- expanse(bathy_madr)/1e+6
bathy_styl <- intersect(bathy_poly, styl_poly)
a_bathy_styl <- expanse(bathy_styl)/1e+6
madr_styl <- intersect(madr_poly, styl_poly)
a_madr_styl <- expanse(madr_styl)/1e+6

# Three-way
abm3 <- intersect(acan_bathy, madr_poly)
a_abm3 <- expanse(abm3)/1e+6
abs3 <- intersect(acan_bathy, styl_poly)
a_abs3 <- expanse(abs3)/1e+6
ams3 <- intersect(acan_madr, styl_poly)
a_ams3 <- expanse(ams3)/1e+6
bms3 <- intersect(bathy_madr, styl_poly)
a_bms3 <- expanse(bms3)/1e+6

# All four
all4 <- intersect(abm3, styl_poly)
a_all4 <- expanse(all4)/1e+6


#Venn diagram
library(eulerr)

fit <- euler(c(
  Acanthogorgia = a_acan,
  Bathypathes = a_bathy,
  Madracis = a_madr,
  Stylaster = a_styl,
  "Acanthogorgia&Bathypathes" = a_acan_bathy,
  "Acanthogorgia&Madracis" = a_acan_madr,
  "Acanthogorgia&Stylaster" = a_acan_styl,
  "Bathypathes&Madracis" = a_bathy_madr,
  "Bathypathes&Stylaster" = a_bathy_styl,
  "Madracis&Stylaster" = a_madr_styl,
  "Acanthogorgia&Bathypathes&Madracis" = a_abm3,
  "Acanthogorgia&Bathypathes&Stylaster" = a_abs3,
  "Acanthogorgia&Madracis&Stylaster" = a_ams3,
  "Bathypathes&Madracis&Stylaster" = a_bms3,
  "Acanthogorgia&Bathypathes&Madracis&Stylaster" = a_all4
), shape="ellipse")

plot(fit, fills = TRUE)

