## Code authored by Megan K. B. Nolan ##

## Associated to the manuscript in review: Projected climate change impacts on deep-sea corals of the Red Sea (in review) ##

## Correcting extent of oceanographic layers: Includes two methods ##

#Set up
setwd("Z:/nolanmc/ENMs/R_BRTs/EnvAll")   #overall folder containing all oceanographic layers

library(raster) ; library(ENMTools) ; library(maps) ; library(remotes) ; library(dismo) ; library(terra)

#Load in rasters and stack ----
files <- list.files(pattern = "\\.tif$")
geo <- list.files(pattern = "\\.asc$")
summary(geo)

for(x in geo){
  a <- raster(x)
  b <- extent(a)
  print(x)
  print(b)
}

depth <- raster("enm_depth.tif.asc")
extent(depth)

env <- raster("NitrateMax.tif")
extent(env)

extent <- c(32.32229, 44.06008, 12.54287,30.11080)

list_clipped <- list()
for(x in files){
  a <- raster(x)
  b <- terra::crop(a, depth, snap="near", mask=T, touches=T)
  print(x)
  print(extent(b))
  #y <- trim(b)
  list_clipped <- append(list_clipped, list(b))
  writeRaster(b, paste0("ascii/", x), format='ascii', overwrite=TRUE)
}

stack <- raster::stack(list_clipped)
stack_crop <- crop(stack, extent(rugosity))
list <- as.list(stack_crop)
print(list)

for(x in list){
  a <- names(x)
  writeRaster(x, paste0("Z:/nolanmc/ENMs/PredictiveLayers/", a, ".asc"), format='ascii', overwrite=TRUE)
}

setwd("Z:/nolanmc/ENMs/PredictiveLayers/")
files <- list.files(pattern = "\\.asc$")

temp <- raster("ENM_temp.asc")
extent(temp)
sal  <- raster("Oceanographic/sal_fin.asc")
ox   <- raster("Oceanographic/oxNew_Fin.asc")

backscatter <- raster("layers_asciiNew/backscatter.asc")

temp2 <- extractByMask(temp, backscatter, spatial=TRUE)
sal2 <- extractByMask(sal, backscatter, spatial=TRUE)
ox2 <- extractByMask(ox, backscatter, spatial=TRUE)

extent(backscatter)
extent(ox2)

writeRaster(temp2, "layers_asciiNew/Temp.asc", format='ascii', overwrite=TRUE)
writeRaster(sal2, "layers_asciiNew/Salinity.asc", format='ascii', overwrite=TRUE)
writeRaster(ox2, "layers_asciiNew/OxygenNew2.asc", format='ascii', overwrite=TRUE)


#________________________________________
rm(list=ls())
setwd("Z:/nolanmc/ENMs/R_BRTs/SSP585")   #overall folder containing all oceanographic layers

library(raster) ; library(ENMTools) ; library(maps) ; library(remotes) ; library(dismo) ; library(terra)

#Load in rasters and stack ----
all <- list.files(pattern = "\\.tif$")    #\\.asc$")
#geo <- raster::stack(all[c(1:2, 6, 16, 20, 27, 37)])
#env <- raster::stack(all[c(3:5, 7:15, 17:19, 21:26, 28:36)])

env <- raster::stack(all)
geo <- raster("depth.tif.asc")

env_crop <- crop(env,extent(geo))
env_all <- unstack(env_crop)
names(env_all) <- names(env)
names(env_all) <- "tempMin"
env_all_list <- as.list(env_all)

for (i in seq_along(env_all)) {
  writeRaster(env_all[[i]], paste0("Z:/nolanmc/ENMs/R_BRTs/SSP585/",names(env_all)[[i]], ".asc"), format="ascii")
}

list_clipped <- list()
for(x in env){
  crop <- crop(x, ref)
  extend <- extend(crop, ref)
  list_clipped[[i]] <-
  a <- rast(x)
  b <- terra::extend(crop(a, ref, snap="near", mask=T, touches=T), ref)
  print(x)
  #y <- trim(b)
  list_clipped <- append(list_clipped, list(b))
  #writeRaster(b, paste0("ascii/", x), format='ascii', overwrite=TRUE)
}

stack <- rast(list_clipped)
stack_crop <- terra::crop(stack, ext(ref))
list <- as.list(stack_crop)
ext(list[[1]])

library(raster)
setwd("Z:/nolanmc/ENMs/R_BRTs/EnvAll") 
files <- list.files(pattern = "\\.tif$")

#setwd("Z:/nolanmc/ENMs/PredictiveLayers/Extent/ASCII") 
for (i in files) {
  tif <- raster(i)
  writeRaster(tif, paste0(i, ".asc"), format="ascii")
}

#If necessary:
#convert tif to ascii
library(raster)
setwd("Z:/nolanmc/ENMs/PredictiveLayers/ExtentNewData") 
files <- list.files(pattern = "\\.tif$")

#setwd("Z:/nolanmc/ENMs/PredictiveLayers/Extent/ASCII") 
for (i in files) {
  tif <- raster(i)
  writeRaster(tif, paste0("ASCII/", i, ".asc"), format="ascii")
}

