#####Bathypathes ####

rm(list=ls()) #----
setwd(".../R_Maxent")
library(devtools) ; library(plyr) ; library(Rcpp) ; library(ENMeval) ; library(dismo) ; library(sp)
library(plotROC) ; library(ggplot2) ; library(sdmpredictors) ; library(SDMtune) ; library(RColorBrewer)
library(terra) ; library(gridExtra) ; library(progress); library(dplyr)

bathy <- read.csv("SpeciesData/bathypathes.csv")[,c(3,2)] %>% na.omit(.) %>% unique(.)
envdata <- rast(list.files("EnvAll/Cropped", pattern="\\.tif.asc$", full.names = T))
environmentalCond <- envdata[[c(1,12,28,30,32,34:35,37:38)]]
new_names <- c("depth", "phMean", "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

presences_all <- bathy
spname <- "Bathypathes"

cell_ids <- cellFromXY(environmentalCond[[1]], presences_all[, c("Long", "Lat")])

#4. Standardisation and k-means clustering for pseudo absences ----
env_df <- as.data.frame(environmentalCond, na.rm=F)
complete_cells <- complete.cases(env_df)
env_df <- env_df[complete_cells,]

comp_cells <- cell_ids[cell_ids %in% row.names(env_df)]
species_n <- length(unique(comp_cells)) #Check how many presences after spatial thinning
presences <- xyFromCell(environmentalCond[[1]], comp_cells)  #make sure presences all have environmental data
presences <- as.data.frame(unique(presences))
colnames(presences) <- c("Long", "Lat")

env_df_scaled <- scale(env_df, center=T, scale=T) #Zero-mean standardisation

k=10*(species_n)
set.seed(144)
kmeans_result <- kmeans(env_df_scaled, centers=k, nstart=10) #k-means clustering
#10 random starts, no. of clusters = 10*species number, changed algorithm as default gave computational warnings

clusters <- rep(NA, nrow(as.data.frame(environmentalCond, na.rm=F)))
clusters[complete_cells] <- kmeans_result$cluster  #add cluster assignments to the valid cells

cluster_raster <- environmentalCond[[1]]
values(cluster_raster) <- clusters
rm(kmeans_result, env_df_scaled, complete_cells, env_df)

writeRaster(cluster_raster, paste0("Present/Plots/", spname, "/clustered_raster", k, ".tif"), overwrite=T)

pdf(paste0("Present/Plots/", spname, "/clustered_raster", k, ".pdf"), width = 7, height = 7)
plot(cluster_raster, col = rainbow(k))
dev.off()
print("clustering complete")

#5. Generate 10 sets of pseudo-absences ----
pseudoAbsencesClusters <- function(rasters, records, n, cluster_raster, n_runs = 10) {
  # Prepare the base raster (e.g., spatial extent)
  shape <- subset(rasters, 1)
  nonNACells <- which(!is.na(terra::values(shape)))
  sink.points <- terra::xyFromCell(shape, nonNACells)
  cluster_values <- terra::values(cluster_raster)
  
  sink.points <- cbind(sink.points, cluster = cluster_values[nonNACells])
  colnames(sink.points) <- c("Long", "Lat", "cluster")
  
  sink.points <- as.data.frame(sink.points)
  
  # Remove points within the buffer zone of presence records
  sink.points.poly <- as.data.frame(records)
  coordinates(sink.points.poly) <- c("Long", "Lat")
  proj4string(sink.points.poly) <- CRS("+proj=longlat +datum=WGS84")
  
  sink.points.polygon <- buffer(sink.points.poly, width = 400 / 111.699, byid = TRUE) #400 m buffer, converted to degrees
  #plot(sink.points.polygon)
  sink.points.pts <- as.data.frame(sink.points[, 1:2])
  colnames(sink.points.pts) <- c("Long", "Lat")
  coordinates(sink.points.pts) <- c("Long", "Lat")
  proj4string(sink.points.pts) <- CRS("+proj=longlat +datum=WGS84")
  
  to.remove.id <- sp::over(sink.points.pts, sink.points.polygon)
  
  if (class(to.remove.id) == "integer") {
    to.keep <- which(is.na(to.remove.id))
  }
  
  if (class(to.remove.id) == "data.frame") {
    if (ncol(to.remove.id) > 1) {
      to.keep <- which(is.na(to.remove.id[, 1]))
    }
    if (ncol(to.remove.id) == 1) {
      to.keep <- which(is.na(to.remove.id))
    }
  }
  
  absences <- as.data.frame(sink.points[to.keep, ])
  absences_clean <- na.omit(absences)
  absences <- absences_clean
  
  cluster_ids <- unique(absences$cluster) # Extract unique clusters
  results <- list() # Initialize a list to store results for each run
  
  for (run_id in 1:n_runs) {
    run_absences <- data.frame() # Data frame to store absences for the current run
    
    for (cluster_id in cluster_ids) {
      # Get points in the current cluster
      cluster_points <- absences[absences$cluster == cluster_id, 1:2]
      
      if (nrow(cluster_points) > 0) {
        # Randomly select one point from the cluster
        selected_point <- cluster_points[sample(1:nrow(cluster_points), 1), ]
        selected_point <- as.data.frame(selected_point) # Ensure it's a data frame
        colnames(selected_point) <- c("Long", "Lat")
        selected_point$cluster <- cluster_id
        selected_point$run_id <- run_id
        
        # Append to the current run's absences
        run_absences <- rbind(run_absences, selected_point)
      }
    }
    
    # Store the run's results in the list
    results[[run_id]] <- run_absences
  }
  
  # Combine all runs into a single data frame
  final_results <- do.call(rbind, results)
  return(final_results)
}
pseudoAbs <- pseudoAbsencesClusters(environmentalCond, presences, n=1000, cluster_raster, n_runs = 10)
print("pseudo absences generated")

# Plot pseudo-absences and presences
world <- shapefile("globalLandmass/world.shp")
RS <- extent(32, 44, 12, 30)
myRegion <- crop(world, RS)

plot(myRegion,col="gray",border="gray")
points(pseudoAbs,  pch= 19)  #col=pseudoAbs$cluster,
points(presences$Long, presences$Lat, col="darkred" , pch= 19)

#6. Fit a model- loop through pseudo absences ----
#Need to test the model first, then run through 10 runs with different pseudoabsences
#Set up model function and empty dataframes
prepareModelData <-  function(p, a, env) {
  # Extract the first layer of the SpatRaster
  m <- env[[1]]
  
  # Remove points with NA values in the extracted layer
  p.i <- terra::extract(m, p)
  p <- p[!is.na(p.i[,2]),]
  
  # Snap points to the nearest cell center and remove duplicates
  p_cells <- terra::cellFromXY(m, p)
  p <- unique(terra::xyFromCell(m, p_cells))
  
  # Remove absence points with NA values
  a.i <- terra::extract(m, a)
  a <- a[!is.na(a.i[,2]),]
  
  # Snap absence points to the nearest cell center and remove duplicates
  a_cells <- terra::cellFromXY(m, a)
  a <- unique(terra::xyFromCell(m, a_cells))
  
  # Call prepareSWD with the cleaned data
  return(prepareSWD(species = "Model species", p = p, a = a, env = env))
}

model_repeats <- list()
accuracy_metrics <- data.frame(species = character(), run_id = integer(), threshold = numeric(), tss=numeric(),
                               sensitivity=numeric(), specificity=numeric(),auc=numeric(), boyce=numeric(),stringsAsFactors = FALSE)
variableImportance <- data.frame(species = character(), run_id = integer(), Variable = character(), 
                                 Permutation_Importance=numeric(), sd=numeric(), stringsAsFactors = FALSE)

setwd("/Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent")

for (i in 1:10){
  pseudoAbs_subset <- subset(pseudoAbs, run_id == i)
  pseudoAbs_subset <- pseudoAbs_subset[,1:2]
  
  modelData <- prepareModelData(presences, pseudoAbs_subset, environmentalCond)
  # Generate cross validation folds
  folds <- get.block(modelData@coords[modelData@pa == 1,], modelData@coords[modelData@pa == 0,])
  
  #Train Maxent
  trainMaxent <- function(data, folds, reg = 2.5, fc = "h", iter = 500) {
    
    result <- SDMmodelCV(data = data, folds = folds)
    
    df <- data@data
    coords <- data@coords
    pa <- data@pa
    
    occs_fold_ids <- unique(folds$occs.grp)
    bg_fold_ids   <- unique(folds$bg.grp)
    
    # Occurrence and background rows are stored sequentially in the SWD object
    n_occs <- sum(pa == 1)
    
    fold_models <- list()
    oof_df <- data.frame()
    
    for (f in seq_along(occs_fold_ids)) {
      
      occs_train_indices <- which(folds$occs.grp != occs_fold_ids[f])
      occs_test_indices  <- which(folds$occs.grp == occs_fold_ids[f])
      bg_train_indices   <- which(folds$bg.grp   != bg_fold_ids[f]) + n_occs
      bg_test_indices    <- which(folds$bg.grp   == bg_fold_ids[f]) + n_occs
      
      train_indices <- c(occs_train_indices, bg_train_indices)
      test_indices  <- c(occs_test_indices,  bg_test_indices)
      
      train_swd <- new("SWD", species = data@species, coords = coords[train_indices, ],
                       data = df[train_indices, ], pa = pa[train_indices])
      
      test_swd  <- new("SWD", species = data@species, coords = coords[test_indices, ],
                       data = df[test_indices, ], pa = pa[test_indices])
      
      result_sdm <- train(method = "Maxent", data = train_swd, reg = reg,  fc = fc, iter = iter)
      
      test_pred <- predict(result_sdm, data = test_swd, type = "logistic")
      oof_df <- rbind(oof_df, data.frame(obs = test_swd@pa, pred = test_pred, fold = f))
      
      result_sdm@data <- test_swd
      fold_models[[f]] <- result_sdm
    }
    
    result@models <- fold_models
    return(list(model = result, oof = oof_df))
  }
  
  #Initial model
  model <- trainMaxent(modelData, folds, reg = 2.5, fc = "h", iter = 500)$model
  
  plot_list <- list()
  for (c in new_names){
    plot <- plotResponse(model = model, var = c, color = "black", marginal = T, rug = F)
    plot_list[[c]] <- plot
  }
  
  ggsave(
    filename = paste0("Present/Plots/", spname, "/MAXVarResponses", spname, i, ".pdf"),
    plot = gridExtra::arrangeGrob(grobs = plot_list, nrow = 7, ncol = 6),
    width = 15, height = 9
  )
  
  # Variable importance and reduction
  model_reduced2 <- trainMaxent(modelData, folds, reg = 2.5, fc = "h", iter = 500)$model
  
  viModel2 <- varImp(model_reduced2, permut = 5)
  final_vars <- viModel2$Variable[viModel2$Permutation_importance > 5]
  
  model_data <- modelData@data
  fin_model_data <- model_data[, final_vars, drop = FALSE]
  
  finModelData <- modelData
  finModelData@data <- fin_model_data
  
  finalModel_out <- trainMaxent(finModelData, folds, reg = 2.5, fc = "h", iter = 500)
  finalModel <- finalModel_out$model
  finalModel_oof <- finalModel_out$oof
  model_repeats[[i]] <- finalModel
  
  # Final variable importance
  viModel <- varImp(finalModel, permut = 5)
  
  variableImportance <- rbind(variableImportance, data.frame(species = spname,
                                                             run_id = rep(i, length(viModel$Variable)), Variable = viModel$Variable,
                                                             Permutation_Importance = viModel$Permutation_importance, sd = viModel$sd))
  
  pdf(paste0("Present/Plots/", spname, "/MAXvariableImp", spname, i, ".pdf"), width = 7, height = 7)
  plotVarImp(viModel)
  dev.off()
  
  # Response curves
  response_plots <- list()
  for (a in unique(viModel$Variable)) {
    responseplot <- plotResponse(model = finalModel, var = a, marginal = T, rug = F, color = "black")
    response_plots[[a]] <- responseplot
  }
  
  ggsave(
    filename = paste0("Present/Plots/", spname, "/MAXResponsePlots", spname, i, ".pdf"),
    plot = gridExtra::arrangeGrob(grobs = response_plots, nrow = 4, ncol = 5),
    width = 15, height = 9
  )
  
  #Model performance
  # determine the threshold maximizing the sum of sensitivity and specificity
  getAccuracy <- function(oof, threshold = 0.26) {
    
    predMain <- oof$pred
    acc <- data.frame()
    
    for(threshold.i in threshold) {
      
      pred <- predMain
      pred[pred >= threshold.i] <- 1
      pred[pred < threshold.i] <- 0
      
      confusionDF <- data.frame(obs = oof$obs, pred)
      sensitivity <- sum(apply(confusionDF[confusionDF$obs == 1,], 1, function(x) 
      { ifelse(x[1] == 1 & x[2] == 1, 1, 0) })) / sum(confusionDF$obs == 1)
      specificity <- sum(apply(confusionDF[confusionDF$obs == 0,], 1, function(x) 
      { ifelse(x[1] == 0 & x[2] == 0, 1, 0) })) / sum(confusionDF$obs == 0)
      
      acc <- rbind(acc, data.frame(threshold = threshold.i, tss = specificity + sensitivity - 1,
                                   sensitivity = sensitivity, specificity = specificity))
    }
    return(acc)
  }
  
  thresholdMaxTSS <- function(oof) {
    
    r <- getAccuracy(oof, threshold = seq(0, 1, by = 0.01))
    
    plot(r$threshold, r$sensitivity, type = "l", xlab = "Threshold", ylab = "Performance", lty = 1, lwd = 1)
    lines(r$threshold, r$specificity, lty = 1, lwd = 4, col = "gray")
    
    val <- unlist(r$threshold)[which.max(unlist(r$sensitivity) + unlist(r$specificity))]
    
    abline(v = val, lty = 3, lwd = 0.7)
    legend(0.7, 0.9, legend = c("Sensitivity", "Specificity"), col = c("black", "gray"), lty = 1, lwd = c(1, 4), cex = 1)
    
    return(val)
  }
  thresholdMaxTSS <- function(model) {
    
    r <- getAccuracy(model,threshold = seq(0,1,by=0.01))
    
    plot(r$threshold,r$sensitivity,type="l", xlab="Threshold",ylab="Performance" , lty=1, lwd=1 )
    lines(r$threshold,r$specificity , lty=1, lwd=4 , col= "gray")
    
    val <- unlist(r$threshold)[which.max( unlist(r$sensitivity) + unlist(r$specificity) )]
    
    abline(v=val , lty=3, lwd=0.7 )
    legend(0.7, 0.9, legend=c("Sensitivity", "Specificity"),col=c("black", "gray"), lty=1, lwd=c(1,4) , cex=1)
    
    return(val)
    
  }
  pdf(paste0("Present/Plots/", spname, "/MAXSensitivity_Specificity", spname, i, ".pdf"), width = 8, height = 6)
  threshold <- thresholdMaxTSS(finalModel_oof)
  dev.off()
  
  # determine model performance
  # determine model performance
  mannWhitneyAUC <- function(obs, pred) {
    r <- rank(pred)
    n_pos <- sum(obs == 1)
    n_neg <- sum(obs == 0)
    if (n_pos == 0 || n_neg == 0) return(NA)
    (sum(r[obs == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }
  
  calcBoyce <- function(oof) {
    boyce_res <- ecospat::ecospat.boyce(
      fit = oof$pred,
      obs = oof$pred[oof$obs == 1],
      PEplot = FALSE
    )
    return(boyce_res$cor)
  }
  
  accuracy <- getAccuracy(finalModel_oof, threshold = threshold)
  per_fold_auc <- sapply(split(finalModel_oof, finalModel_oof$fold),
                         function(d) mannWhitneyAUC(d$obs, d$pred))
  auc <- mean(per_fold_auc, na.rm = TRUE)
  boyce <- calcBoyce(finalModel_oof)
  
  accuracy_metrics <- rbind(accuracy_metrics, data.frame(species = spname, run_id = i, threshold = accuracy$threshold, 
                                                         tss = accuracy$tss, sensitivity=accuracy$sensitivity, specificity=accuracy$specificity, auc=auc, boyce=boyce))}
#THESE NEED TO BE ADJUSTED FOR GREATER NUMBER OF RUNS
accuracy_metrics[11,] <- c(spname, "Average",(colSums(accuracy_metrics[1:10,3:8])/10))
write.csv(accuracy_metrics, paste0("Present/Plots/", spname, "/MAXaccuracy_", spname, ".csv"))
write.csv(variableImportance, paste0("Present/Plots/", spname, "/MAXvariables_", spname, ".csv"))

variables <- tidyr::pivot_wider(variableImportance[,1:4], names_from = "run_id", values_from= "Permutation_Importance")
variables[is.na(variables)] <- 0
variables$average <- rowSums(variables[3:12])/10

write.csv(variables, paste0("Present/Plots/", spname, "/MAXvariables_summary_", spname, ".csv"))

#7. Predict to produce maps ----
predict_rasters <- list()
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Present/Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Present/Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Present/Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)

#8. Different climate scenarios: SSP119 ----
setwd("/Volumes/WD_Megan/ENMs/R_BRTs")
predict_rasters <- list()
envdata <- rast(list.files("SSP119", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata[[c(1,6, 10, 12, 14:18)]]
new_names <- c("depth", "phMean",  "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names
setwd("/Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent/SSP119")
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)

#9. Different climate scenarios: SSP126 ----
setwd("/Volumes/WD_Megan/ENMs/R_BRTs")
predict_rasters <- list()
envdata <- rast(list.files("SSP126", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata[[c(1,6, 10, 12, 14:18)]]
new_names <- c("depth", "phMean",  "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names
setwd("/Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent/SSP126")
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)

#10. Different climate scenarios: SSP245 ----
setwd("/Volumes/WD_Megan/ENMs/R_BRTs")
predict_rasters <- list()
envdata <- rast(list.files("SSP245", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata[[c(1,6, 10, 12, 14:18)]]
new_names <- c("depth", "phMean",  "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names
setwd("/Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent/SSP245")
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)

#11. Different climate scenarios: SSP370 ----
setwd("/Volumes/WD_Megan/ENMs/R_BRTs")
predict_rasters <- list()
envdata <- rast(list.files("SSP370", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata[[c(1,6, 10, 12, 14:18)]]
new_names <- c("depth", "phMean",  "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names
setwd("/Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent/SSP370")
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)

#12. Different climate scenarios: SSP585 ----
setwd("/Volumes/WD_Megan/ENMs/R_BRTs")
predict_rasters <- list()
envdata <- rast(list.files("SSP585", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata[[c(1,6, 10, 12, 14:18)]]
new_names <- c("depth", "phMean",  "slope", "tempMean", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names
setwd("//Users/megannolan/Desktop/GCB_Review/Results_August/R_Maxent/SSP585")
for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  #predict_map <- terra::plot(mapPresent, col = terrain.colors(100))
  writeRaster(mapPresent, paste0("Plots/", spname, "/MAXmap", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("Plots/", spname, "/MAXmean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("Plots/", spname, "/MAXsd", spname, ".tif"), overwrite = TRUE)
