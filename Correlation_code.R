## Code authored by Megan K. B. Nolan, 2025 ##

## Associated to the manuscript in review: Projected climate change impacts on deep-sea corals of the Red Sea (in review) ##

## Correlation test and matrix, and creation of Principal Components to combined highly correlated variables ##

envdata <- list.files("Environmental", pattern="\\.tif.asc$", full.names = T)
stack_env <- rast(envdata)
new_names <- c("bpi17", "bpi25", "bpi3", "bpi33", "bpi65", "bpi9", #1-6
               "ci3", "ci5", "ci9", "cx3", "cx5", "cx9", "depth", "east", "meancurve", #7-15
               "nitrate", "north", "oxy", "ph", "phosphate", "plan", "primaryprod", #16-22
               "profile", "sal", "silicate", "slope", "temp", "tex3", "tex5", "tex9", #23-30
               "velx", "vely", "vrm3", "vrm5", "vrm9") #31-35
names(stack_env) <- new_names
stack_env2 <- stack_env[[-c(29, 30)]] #Also correct with monotonicity: removed tex5 and tex9

#2A. Correlations: SKIP ----
#plot(stack_env, col=terrain.colors(1000))
#correlation_matrix <- layerCor(stack_env, fun = "pearson", na.rm = TRUE)
#print(corr.matrix)
#Correlations and progress bar
n_layers <- nlyr(stack_env)
correlation_matrix <- matrix(NA, nrow = n_layers, ncol = n_layers)
rownames(correlation_matrix) <- names(stack_env)
colnames(correlation_matrix) <- names(stack_env)

# Progress bar setup
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed: :elapsed eta: :eta",
  total = n_layers * (n_layers - 1) / 2, # Total comparisons
  clear = FALSE
)

# Calculate correlations
for (i in 1:(n_layers - 1)) {
  for (j in (i + 1):n_layers) {
    # Extract pairwise layers
    layer1 <- stack_env[[i]]
    layer2 <- stack_env[[j]]
    
    # Calculate correlation and store it
    correlation_matrix[i, j] <- cor(values(layer1), values(layer2), use = "complete.obs", method = "pearson")
    correlation_matrix[j, i] <- correlation_matrix[i, j] # Symmetric matrix
    
    # Update progress bar
    pb$tick()
  }
}

# Convert to a data frame for easier viewing (optional)
correlation_matrix <- as.data.frame(correlation_matrix)

# Print or save the result
print(correlation_matrix)
write.csv(correlation_matrix, "correlation_matrix.csv")


  #PCA creation for correlated variables ----
rm(list=ls())
setwd("Z:/nolanmc/ENMs/R_BRTs")   #overall folder containing all oceanographic layers
phos <- rast("SSP585/phosphateMean.asc")
nitrate <- rast("SSP585/nitrateMean.asc")
envdata <- rast(list.files("SSP585", pattern="\\.asc$", full.names = T))


values_matrix <- values(envdata, mat = TRUE)
complete_cases <- complete.cases(values_matrix)
values_matrix <- values_matrix[complete_cases, ]

# Standardize the data
scaled_values <- as.data.frame(scale(values_matrix, center=T, scale=T))

# Perform PCA on specific combinations
#Nutrients: phosphate and nitrate
nutrients1 <- scaled_values[,c(3, 7)] #select whichever are the nutrients
colnames(nutrients1)

pca_result <- prcomp(nutrients1, center = TRUE, scale. = TRUE)
summary(pca_result)
pca_scores <- pca_result$x[, 1:2]  #get 1st two PCA scores

# Map PCA scores back to raster format
template <- envdata[[1]] #empty raster template

# Create new rasters for each PC
pc1_raster <- template
values(pc1_raster)[complete_cases] <- pca_scores[, 1]
writeRaster(pc1_raster, filename = "Z:/nolanmc/ENMs/R_BRTs/SSP585/nutrients_pca1.tif", overwrite = TRUE)

pca1 <- raster("Z:/nolanmc/ENMs/R_BRTs/SSP585/nutrients_pca1.tif")
writeRaster(pca1, "Z:/nolanmc/ENMs/R_BRTs/SSP585/nutrients_pca1.asc", format="ascii", overwrite = TRUE)



pc2_raster <- template
values(pc2_raster)[complete_cases] <- pca_scores[, 2]
writeRaster(pca_rasters, filename = "nutrients_pca2.tif", overwrite = TRUE)

#Oxygen and Salinity
nutrients1 <- values_matrix[[c(29, 30)]] #select whichever are the correct variables

pca_result <- prcomp(nutrients1, center = TRUE, scale. = TRUE)
summary(pca_result)
pca_scores <- pca_result$x[, 1:2]  #get 1st two PCA scores

# Map PCA scores back to raster format
template <- raster_stack[[1]] #empty raster template

# Create new rasters for each PC
pc1_raster <- template
values(pc1_raster)[complete_cases] <- pca_scores[, 1]
writeRaster(pca_rasters, filename = "oxsal_pca1.tif", overwrite = TRUE)

pc2_raster <- template
values(pc2_raster)[complete_cases] <- pca_scores[, 2]
writeRaster(pca_rasters, filename = "oxsal_pca2.tif", overwrite = TRUE)


# Combine PCs into a new SpatRaster object
pca_rasters <- c(pc1_raster, pc2_raster)
names(pca_rasters) <- c("PC1", "PC2")


# Visualise correlation matrix to create supplementary figure ----

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)

setwd("Z:/nolanmc/ENMs/R_BRTs") 

# Read the CSV file
cor_matrix <- read.csv("full_correlation_matrix.csv", row.names = 1, check.names = FALSE)
new_names <- c(
  "velyMin" = "Min. Velocity Y",
  "velyMean" = "Mean Velocity Y",
  "velyMax" = "Max. Velocity Y",
  "velxMin" = "Min. Velocity X",
  "velxMax" = "Max. Velocity X",
  "velxMean" = "Mean Velocity X",
  "tempMin" = "Min. Temperature",
  "tempMean" = "Mean Temperature",
  "tempMax" = "Max. Temperature",
  "slope" = "Slope",
  "SilicateMin" = "Min. Silicate",
  "SilicateMean" = "Mean Silicate",
  "silicateMax" = "Max. Silicate",
  "salMin" = "Min. Salinity",
  "salMean" = "Mean Salinity",
  "salMax" = "Max. Salinity",
  "profile" = "Profile Curvature",
  "prodMin" = "Min. Productivity",
  "prodMean" = "Mean Productivity",
  "prodMax" = "Max. Productivity",
  "plan" = "Planar Curvature",
  "phosMin" = "Min. Phosphate",
  "phosMean" = "Mean Phosphate",
  "phosMax" = "Max. Phosphate",
  "phMin" = "Min. pH",
  "phMean" = "Mean pH",
  "phMax" = "Max pH",
  "OxMin" = "Min. Oxygen",
  "OxMean" = "Mean Oxygen",
  "oxMax" = "Max. Oxygen",
  "north" = "Northness",
  "nitrateMin" = "Min. Nitrate",
  "nitrateMean" = "Mean Nitrate",
  "nitrateMax" = "Max. Nitrate",
  "east" = "Eastness",
  "depth" = "Depth",
  "vrm3" = "VRM"
  # Add all your variable name mappings here
)

rownames(cor_matrix) <- new_names[rownames(cor_matrix)]
colnames(cor_matrix) <- new_names[colnames(cor_matrix)]


# Convert to matrix
cor_matrix <- as.matrix(cor_matrix)

# Mask the upper triangle
cor_matrix[upper.tri(cor_matrix)] <- NA

# Melt into long format
cor_melted <- melt(cor_matrix, na.rm = TRUE)

# Plot the heatmap
heatmap <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 1)), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limit = c(-1, 1), space = "Lab",
    name = "Correlation"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=11),
    axis.title = element_blank(), axis.text.y=element_text(size=11)) +
  coord_fixed()
  
ggsave("CorrelationMatrix.pdf", plot=heatmap, width=12, height=12)



