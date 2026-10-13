## Code authored by Megan K. B. Nolan and Jorge Assis, 2025 ##

## Associated to the manuscript in review: Projected climate change impacts on deep-sea corals of the Red Sea (in review) ##

## Ecological niche models for Red Sea corals ##

#1. Set up and read in data #----
rm(list=ls()) 
setwd(".../R_BRTs")
library(devtools) ; library(plyr) ; library(Rcpp) ; library(ENMeval) ; library(dismo) ; library(sp)
library(plotROC) ; library(ggplot2) ; library(sdmpredictors) ; library(SDMtune) ; library(RColorBrewer)
library(terra) ; library(gridExtra) ; library(progress)

presences_all <- read.csv("SpeciesData/Acanthogorgia.csv")[,c(3,2)] %>% na.omit(.) %>% unique(.) #Read in occurrences for Acanthogorgia spp.
spname <- "Acanthogorgia"

envdata <- rast(list.files("EnvAll/Cropped", pattern="\\.tif.asc$", full.names = T)) #Read in all oceanographic and physical data as rasters
environmentalCond <- envdata[[c(1:2,6:7,10,12,17,21,23,28,29,31,32,34:35,37:38)]]
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope","tempMax", "tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") #Select relevant rasters
names(environmentalCond) <- new_names

cell_ids <- cellFromXY(environmentalCond[[1]], presences_all[, c("Long", "Lat")])
  
#2. Standardisation and k-means clustering for pseudo absences ----
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
  
#3. Generate 10 sets of pseudo-absences ----
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

#4. Fit a model with best hyperparameters- loop through pseudo absences ----
#Need to test the model first, then run through 10 runs with different pseudoabsences
#Set up model function and empty dataframes
prepareModelData <- function(p, a, env) {
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
                               sensitivity=numeric(), specificity=numeric(),auc=numeric(), boyce=numeric(), stringsAsFactors = FALSE)
variableImportance <- data.frame(species = character(), run_id = integer(), Variable = character(), 
                                 Permutation_Importance=numeric(), sd=numeric(), stringsAsFactors = FALSE)
parameters <- data.frame(species = character(), run_id = integer(), interaction.depth = character(), 
                         shrinkage=numeric(),stringsAsFactors = FALSE)
#Function to calculate Boyce's Index
calcBoyce <- function(oof) {
  boyce_res <- ecospat::ecospat.boyce(
    fit = oof$pred,
    obs = oof$pred[oof$obs == 1],
    PEplot = FALSE
  )
  return(boyce_res$cor)
}

for (i in 1:10){
  print(paste("Starting run ", i))

  pseudoAbs_subset <- subset(pseudoAbs, run_id == i)
  pseudoAbs_subset <- pseudoAbs_subset[,1:2]
  
  modelData <- prepareModelData(presences, pseudoAbs_subset, environmentalCond)
  
  # Generate cross validation folds
  folds <- get.block(modelData@coords[modelData@pa == 1,], modelData@coords[modelData@pa == 0,])
  
  # define monotonicity constrains (-1 for negative, +1 for positive, 0 for non-monotonicity)
  monotonicity = data.frame(depth=0,phMean=0,slope=1,tempMax=-1,tempMin=1,velxMax=-1,velxMin=1,
                            velyMax=-1,velyMin=1,vrm=0)
  
  #Train BRT
  trainBRT <- function(data, folds, distribution = "bernoulli", n.trees = 1000, interaction.depth = 1, 
                         shrinkage = 0.1, bag.fraction = 0.5, 
                         n.minobsinnode = 10, var.monotone=monotonicity) {
      
      monotonicity <- var.monotone
      result <- SDMmodelCV(data = data, folds=folds)
      
      df <- data@data
      df <- cbind(pa = data@pa, df)
      # Prepare the custom folds using occs.grp and bg.grp
      occs_fold_ids <- unique(folds$occs.grp)  # Unique fold ids for occurrences
      bg_fold_ids <- unique(folds$bg.grp)      # Unique fold ids for background
      
      # Placeholder for storing models for each fold
      fold_models <- list()
      
      # Train a separate model for each fold
     oof_df <- data.frame() 
     for (i in seq_along(occs_fold_ids)) {
        # Get the train and test indices for occurrences
        occs_train_indices <- which(folds$occs.grp != occs_fold_ids[i])
        occs_test_indices <- which(folds$occs.grp == occs_fold_ids[i])
        
        # Get the train and test indices for background
        bg_train_indices <- which(folds$bg.grp != bg_fold_ids[i])
        bg_test_indices <- which(folds$bg.grp == bg_fold_ids[i])
        
        # Combine the occurrence and background data for training and testing
        train_indices <- c(occs_train_indices, bg_train_indices)
        test_indices <- c(occs_test_indices, bg_test_indices)
        
        train_data <- df[train_indices, ]
        test_data <- df[test_indices, ]
        result_sdm <- SDMmodel(data = data)  
        
        #print(var.monotone)
        
       model <- gbm::gbm(pa ~ ., data = train_data, distribution = distribution,
                          n.trees = n.trees, interaction.depth = interaction.depth, 
                          shrinkage = shrinkage, bag.fraction = bag.fraction, n.minobsinnode = n.minobsinnode,
                          var.monotone=monotonicity)
        #Added for folds
        test_pred <- gbm::predict.gbm(model, newdata = test_data, n.trees = n.trees, type = "response")
        oof_df <- rbind(oof_df, data.frame(obs = test_data$pa, pred = test_pred, fold = i))
        
        
        model_object <- BRT(n.trees = n.trees, distribution = distribution,
                            interaction.depth = interaction.depth,
                            shrinkage = shrinkage, bag.fraction = bag.fraction,
                            model = model)
        
        result_sdm@model <- model_object
        fold_models[[i]] <- result_sdm
      }
      
      result@models <- fold_models
      return(list(model = result, oof = oof_df))
    }
    
  model <- trainBRT(modelData, folds, distribution = "bernoulli", n.trees = 1000, interaction.depth = 1, shrinkage = 0.1,
                      bag.fraction = 0.5, n.minobsinnode = 10, var.monotone=monotonicity)
    
  plot_list <- list()
    for (c in new_names){
     plot <- plotResponse(model=model, var=c, color="black", marginal=T, rug=F)
     plot_list[[c]] = plot
  }
  
  ggsave(
    filename = paste0("Present/Plots/", spname, "/VarResponses", spname, i, ".pdf"), 
    plot = gridExtra::arrangeGrob(grobs=plot_list, nrow=7, ncol=6), 
    width = 15, height = 9
  )
  
  # given a set of possible hyperparameter values for BRT
  h <- list(interaction.depth = c(1,2,3,4) , shrinkage = c(0.1,0.05,0.01,0.005,0.001))
  
  # test all possible combinations of hyperparameter values
  exp1 <- gridSearch(model, hypers = h, metric = "auc")
  plot(exp1)
  results <- exp1@results
  
  best_params <- results[which.max(results$test_AUC), c("interaction.depth", "shrinkage")]
  
  # Extract the best parameters
  best_depth <- as.numeric(best_params["interaction.depth"])
  best_shrinkage <- as.numeric(best_params["shrinkage"])
  parameters <- rbind(parameters, data.frame(speciers=spname, run_id= i, 
                                             interaction.depth = best_depth, shrinkage=best_shrinkage))
  
  # fit a BRT model to the dataset with the best hyperparameter values
  model <- trainBRT(modelData, folds, distribution = "bernoulli", n.trees = 1000, 
                    interaction.depth = best_depth, shrinkage = best_shrinkage,
                    bag.fraction = 0.5, n.minobsinnode = 10, var.monotone=monotonicity)
  
  #Assess for variable contribution and response functions ----
  reduced_monotonicity <- monotonicity
  
  ModelData_reduced <- modelData
  
  model_reduced2 <- trainBRT(ModelData_reduced,folds, distribution = "bernoulli", n.trees = 1000, 
                             interaction.depth = best_depth, shrinkage = best_shrinkage,
                             bag.fraction = 0.5, n.minobsinnode = 10, 
                             var.monotone = reduced_monotonicity)
  
  # determine the final relative variable contribution
  viModel2 <- varImp(model_reduced2, permut = 5)
  
  final_vars <- viModel2$Variable[viModel2$Permutation_importance > 5]
  final_monotonicity <- monotonicity[names(monotonicity) %in% final_vars]
  
  
  model_data <- modelData@data
  fin_model_data <- model_data[,final_vars, drop=F]
  
  finModelData <- modelData
  finModelData@data <- fin_model_data
  
  finalModel_out <- trainBRT(finModelData,folds, distribution = "bernoulli", n.trees = 1000, 
                         interaction.depth = best_depth, shrinkage = best_shrinkage,
                         bag.fraction = 0.5, n.minobsinnode = 10, 
                         var.monotone = final_monotonicity)
  finalModel <- finalModel_out$model
  finalModel_oof <- finalModel_out$oof
  model_repeats[[i]] <- finalModel
  
  # determine the final relative variable contribution
  viModel <- varImp(finalModel, permut = 5)
  
  variableImportance <- rbind(variableImportance, data.frame(species= spname, run_id = rep(i, length(viModel$Variable)), 
                     Variable = viModel$Variable, Permutation_Importance = viModel$Permutation_importance, sd=viModel$sd))
  
  pdf(paste0("Present/Plots/", spname, "/variableImp", spname,i, ".pdf"), width = 7, height = 7)
  plotVarImp(viModel)
  dev.off()
  
  #response curves
  response_plots <- list()
  for(a in unique(viModel$Variable)){
    
    responseplot <- plotResponse(model=finalModel, var = a, marginal = T, rug = F, color="Black")
    response_plots[[a]] = responseplot
    
  }
  
  ggsave(
    filename = paste0("Present/Plots/", spname, "/ResponsePlots", spname,i, ".pdf"), 
    plot = gridExtra::arrangeGrob(grobs=response_plots, nrow=4, ncol=5), 
    width = 15, height = 9
  )
  
  #Model performance
  # determine the threshold maximizing the sum of sensitivity and specificity
   getAccuracy <- function(oof,threshold=0.26) {
    
    predMain <- oof$pred 
    acc <- data.frame()
    
    for(threshold.i in threshold ) {
      
      pred <- predMain
      pred[pred >= threshold.i] <- 1
      pred[pred < threshold.i] <- 0
      
      confusionDF <- data.frame(obs=oof$obs,pred)
      sensitivity <- sum(apply(confusionDF[confusionDF$obs == 1,],1,function(x) 
      { ifelse(x[1] == 1 & x[2] == 1 , 1 , 0 ) })) / sum(confusionDF$obs == 1)
      specificity <- sum(apply(confusionDF[confusionDF$obs == 0,],1,function(x) 
      { ifelse(x[1] == 0 & x[2] == 0 , 1 , 0 ) })) / sum(confusionDF$obs == 0)
      
      acc <- rbind(acc,data.frame(threshold=threshold.i,tss=specificity+sensitivity-1,
                                  sensitivity=sensitivity,specificity=specificity))
      
    }
    return(acc)
    
  }   #Change threshold as required
  thresholdMaxTSS <- function(oof) {
    
    r <- getAccuracy(oof,threshold = seq(0,1,by=0.01))
    
    plot(r$threshold,r$sensitivity,type="l", xlab="Threshold",ylab="Performance" , lty=1, lwd=1 )
    lines(r$threshold,r$specificity , lty=1, lwd=4 , col= "gray")
    
    val <- unlist(r$threshold)[which.max( unlist(r$sensitivity) + unlist(r$specificity) )]
    
    abline(v=val , lty=3, lwd=0.7 )
    legend(0.7, 0.9, legend=c("Sensitivity", "Specificity"),col=c("black", "gray"), lty=1, lwd=c(1,4) , cex=1)
    
    return(val)
    
  }
  pdf(paste0("Present/Plots/", spname, "/Sensitivity_Specificity", spname, i, ".pdf"), width = 8, height = 6) #check this
  threshold <- thresholdMaxTSS(finalModel_oof)
  dev.off()
  
  # determine model performance
  accuracy <- getAccuracy(finalModel_oof,threshold = threshold)
  mannWhitneyAUC <- function(obs, pred) {
    r <- rank(pred)
    n_pos <- sum(obs == 1)
    n_neg <- sum(obs == 0)
    if (n_pos == 0 || n_neg == 0) return(NA)
    (sum(r[obs == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }
  
  per_fold_auc <- sapply(split(finalModel_oof, finalModel_oof$fold),
                         function(d) mannWhitneyAUC(d$obs, d$pred))
  auc <- mean(per_fold_auc, na.rm = TRUE)
  
  boyce <- calcBoyce(finalModel_oof)
  
  accuracy_metrics <- rbind(accuracy_metrics, data.frame(species = spname, run_id = i, threshold = accuracy$threshold, 
               tss = accuracy$tss, sensitivity=accuracy$sensitivity, specificity=accuracy$specificity,auc=auc, boyce=boyce))

}

#Find average and SD of predictions and save
accuracy_metrics[11,] <- c(spname, "Average",(colSums(accuracy_metrics[1:10,3:8])/10))
write.csv(accuracy_metrics, paste0("Present/Plots/", spname, "/accuracy_", spname, ".csv"))
write.csv(parameters, paste0("Present/Plots/", spname, "/parameters_", spname, ".csv"))
write.csv(variableImportance, paste0("Present/Plots/", spname, "/variables_", spname, ".csv"))

variables <- tidyr::pivot_wider(variableImportance[,1:4], names_from = "run_id", values_from= "Permutation_Importance")
variables[is.na(variables)] <- 0
variables$average <- rowSums(variables[3:12])/10

write.csv(variables, paste0("Present/Plots/", spname, "/variables_summary_", spname, ".csv"))

#6. Different climate scenarios: SSP119 ----
predict_rasters <- list()
envdata <- rast(list.files("SSP119", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope", "tempMax","tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  writeRaster(mapPresent, paste0("SSP119/Plots/", spname, "/map", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("SSP119/Plots/", spname, "/mean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("SSP119/Plots/", spname, "/sd", spname, ".tif"), overwrite = TRUE)

#7. Different climate scenarios: SSP126 ----

predict_rasters <- list()
envdata <- rast(list.files("SSP126", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope", "tempMax","tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  writeRaster(mapPresent, paste0("SSP126/Plots/", spname, "/map", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("SSP126/Plots/", spname, "/mean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("SSP126/Plots/", spname, "/sd", spname, ".tif"), overwrite = TRUE)

#8. Different climate scenarios: SSP245 ----
predict_rasters <- list()
envdata <- rast(list.files("SSP245", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope", "tempMax","tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  writeRaster(mapPresent, paste0("SSP245/Plots/", spname, "/map", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("SSP245/Plots/", spname, "/mean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("SSP245/Plots/", spname, "/sd", spname, ".tif"), overwrite = TRUE)

#9. Different climate scenarios: SSP370 ----
predict_rasters <- list()
envdata <- rast(list.files("SSP370", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope", "tempMax","tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  writeRaster(mapPresent, paste0("SSP370/Plots/", spname, "/map", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("SSP370/Plots/", spname, "/mean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("SSP370/Plots/", spname, "/sd", spname, ".tif"), overwrite = TRUE)

#10. Different climate scenarios: SSP585 ----
predict_rasters <- list()
envdata <- rast(list.files("SSP585", pattern="\\.asc$", full.names = T))
environmentalCond <- envdata
new_names <- c("depth", "east", "north", "nutrients", "OxMin",  "phMean", "plan","profile", "salMean",
               "slope", "tempMax","tempMin", "velxMax", "velxMin", "velyMax", "velyMin", "vrm") 
names(environmentalCond) <- new_names

for (b in 1:10){
  model_ob <- model_repeats[[b]]
  mapPresent <- terra::predict(model_ob, environmentalCond, type=c("logistic"))
  predict_rasters[[b]] <- mapPresent
  writeRaster(mapPresent, paste0("SSP585/Plots/", spname, "/map", spname, b, ".tif"), overwrite = TRUE)
}

#Find average and SD of predictions and save
stack_pred <- rast(predict_rasters)
mean_pred <- app(stack_pred, fun=mean, na.rm=T)
writeRaster(mean_pred, paste0("SSP585/Plots/", spname, "/mean", spname, ".tif"), overwrite = TRUE)
sd_pred <- app(stack_pred, fun=sd, na.rm=T)
writeRaster(sd_pred, paste0("SSP585/Plots/", spname, "/sd", spname, ".tif"), overwrite = TRUE)


