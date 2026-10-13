## Code authored by Megan K. B. Nolan ##

## Associated to the manuscript in review: Projected climate change impacts on deep-sea corals of the Red Sea (in review) ##

## Boxplots showing background data (whole study area) vs presences of each coral- ultimately replced code from "mkbnolan/EcologicalNicheModels/environmental_boxplots.R" ##

## Supplementary Figure 4 ##


## Creating bar plots of species conditions in highly suitable areas (not only presences)
library(terra) ; library(dplyr) ;library(ggplot2) ; library(tidyr) ; library(cowplot)
setwd("Z:/nolanmc/ENMs/R_BRTs/Present/Plots")

  #Read in coral and environmental data
speciesname <- list("Acanthogorgia","Bathypathes", "Madracis", "Stylaster")
thresholds <- data.frame(Acanthogorgia = 0.237, Bathypathes = 0.195, Madracis = 0.391, Stylaster = 0.232)

slope <- rast("Z:/nolanmc/ENMs/R_BRTs/EnvAll/Cropped/slope.tif.asc")
depth <- rast("Z:/nolanmc/ENMs/R_BRTs/EnvAll/Cropped/depth.tif.asc")
minTemp <- rast("Z:/nolanmc/ENMs/R_BRTs/EnvAll/Cropped/TempMin.tif.asc")
maxTemp <- rast("Z:/nolanmc/ENMs/R_BRTs/EnvAll/Cropped/TempMax.tif.asc")
vrm <- rast("Z:/nolanmc/ENMs/R_BRTs/EnvAll/Cropped/vrm.tif.asc")

depth_df <- data.frame(species=character(), depth=numeric())
slope_df <- data.frame(species=character(), slope=numeric())
minTemp_df <- data.frame(species=character(), minTemp=numeric())
maxTemp_df <- data.frame(species=character(), maxTemp=numeric())
vrm_df <- data.frame(species=character(), vrm=numeric())

  #Cycle through each coral species
for (a in speciesname){
  threshold = thresholds[[a]] 
  
  presenceraster <- rast(paste0(a,"/mean", a,".tif")) #Based on model results
  crs(presenceraster) <- "EPSG:4326" #ensure correct CRS
  binary_raster <- ifel(presenceraster > threshold, 1, 0)
  
  depth_sp <- values(depth)[values(binary_raster) == 1] #Find depth where species presence is 1
  depth_df1 <- data.frame(depth = depth_sp, species=a)%>%
    na.omit(.)
  depth_df <- rbind(depth_df, depth_df1)
  
  slope_sp <- values(slope)[values(binary_raster) == 1]
  slope_df1 <- data.frame(slope = slope_sp, species=a)%>%
    na.omit(.)
  slope_df <- rbind(slope_df, slope_df1)
  
  minTemp_sp <- values(minTemp)[values(binary_raster) == 1]
  minTemp_df1 <- data.frame(minTemp = minTemp_sp, species=a)%>%
    na.omit(.)
  minTemp_df <- rbind(minTemp_df, minTemp_df1)
  
  maxTemp_sp <- values(maxTemp)[values(binary_raster) == 1]
  maxTemp_df1 <- data.frame(maxTemp = maxTemp_sp, species=a)%>%
    na.omit(.)
  maxTemp_df <- rbind(maxTemp_df, maxTemp_df1)
  
  vrm_sp <- values(vrm)[values(binary_raster) == 1]
  vrm_df1 <- data.frame(vrm = vrm_sp, species=a)%>%
    na.omit(.)
  vrm_df <- rbind(vrm_df, vrm_df1)
}

  #Create all boxplots
depth_plot <- ggplot(depth_df, aes(x = species, y = depth)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  ylab("Depth (m)")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

slope_plot <- ggplot(slope_df, aes(x = species, y = slope)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  ylab("Slope (°)")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

minTemp_plot <- ggplot(minTemp_df, aes(x = species, y = minTemp)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  ylab("Min. Temperature (°C)")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

maxTemp_plot <- ggplot(maxTemp_df, aes(x = species, y = maxTemp)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  xlab("Taxa") +
  ylab("Max. Temperature (°C)")+
  scale_x_discrete(labels=c(expression(italic("Acanthogorgia") ~ "spp."),
                            expression(italic("Bathypathes thermophila")), 
                            expression(italic("Madracis interjecta")),
                            expression(italic("Stylaster tritoni"))))

vrm_plot <- ggplot(vrm_df, aes(x = species, y = vrm)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  xlab("Taxa") +
  ylab("Vector Ruggedness Measure")+
  scale_x_discrete(labels=c(expression(italic("Acanthogorgia") ~ "spp."),
                            expression(italic("Bathypathes thermophila")), 
                            expression(italic("Madracis interjecta")),
                            expression(italic("Stylaster tritoni"))))

plot_list <- list(depth_plot, slope_plot, minTemp_plot, maxTemp_plot, vrm_plot)
custom_labels <- c(paste0(letters[1:5], ")"))
combined_plot <- plot_grid(plotlist = plot_list, ncol=2, nrow=3, labels=custom_labels, align="hv")
plot(combined_plot)
ggsave(plot=combined_plot, "species_env_data.pdf", width=12, height=9)

