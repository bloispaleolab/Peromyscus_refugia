#### Peromyscus maniculatus species group modeling ##########
options(java.parameters = "-Xmx8000m")
library("raster")
library("maptools")
library("spThin")
library("dismo")
library("rJava")
library("ENMeval")
library("rgeos")
library("rgdal")

#set working directory 
wd <- "."
setwd(wd)

# load in peromyscus maniculatus lineage 1 data 
PMK <- read.csv("data/Pk_coords.csv", header=T)[,c(1,3:4)]
PMS <- read.csv("data/PMS_coords.csv", header=T)
PMG <- read.csv("data/PMG_coords.csv", header=T)

#load evironemntal variables 

env_w <- list.files(path="data/env/0BP/", pattern='tif', full.names=TRUE)
WC <- stack(env_w)
projection(WC) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
WC
plot(WC, 1)
points(PMK$LONG, PMK$LAT, pch=16, col="red", cex=0.5)
points(PMS$LONG, PMS$LAT, pch=16, col="yellow", cex=0.5)
points(PMG$LONG, PMG$LAT, pch=16, col="blue", cex=0.5)


## function to create a minimum convex polygon 
# simpleMCP written by Jamie Kass spring 2014 
# Rob Boria vetted code, spring 2014, during his second chapter analyses by comparing code to MCP created in ArcGIS
simpleMCP <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}


#Clip WC to NA 
data(wrld_simpl)

# plot the data

North_AM <- wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"),
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),]

PM_locs <- list(PMK[3:2], PMS[3:2], PMG[3:2])

# Generate MCP for each species   
PM_mcp <- lapply(PM_locs, simpleMCP)
plot(PM_mcp[[1]], add=TRUE, col= "red")
# Generate a buffered MCP (here 5.0 degrees)
PM_buf <- lapply(PM_mcp, gBuffer, width=5.0)
plot(PM_buf[[1]], add=T)

# Generate a buffered MCP (here 7.5 degrees) for projection 
PM_buf_pro <- lapply(PM_mcp, gBuffer, width=7.5)
plot(PM_buf_pro[[1]], add=T)

#clip layers by full extent 
layers_pk <- crop(WC, PM_buf[[1]])
Pk_pres <- mask(layers_pk, PM_buf[[1]])

plot(Pk_pres[[1]])
#plot(north_am, 1)
points(PMK$LONG, PMK$LAT, pch=16, col="red", cex=0.5)


#clip layers by full extent 
layers_pk_pro <- crop(WC, PM_buf_pro[[1]])
Pk_pres_pro <- mask(layers_pk_pro, PM_buf_pro[[1]])

plot(Pk_pres[[1]])
plot(Pk_pres_pro[[1]])

##########################################################################
##Peromyscus keeni modeling 
#spatial filter in spThin
PM_thin <- thin(loc.data = PMK, 
                lat.col = "LAT", long.col = "LONG", 
                spec.col = "ID", 
                thin.par = 25, reps = 100, 
                locs.thinned.list.return = TRUE, 
                write.files = TRUE, 
                max.files = 5,
                out.dir = "data/thin_PMK/", out.base = "PMK_f", 
                write.log.file = TRUE,
                log.file = "data/thin_PMK/PM_thin_log.txt")

#load in spatial thin dataset
PMK_F <- read.csv("data/thin_PMK/PMK_f_thin1.csv")

#plot 
plot(Pk_pres[[1]])
points(PMK$LONG, PMK$LAT, pch=16, col="black", cex=0.5)
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="red", cex=.5)

#lineage one model testing 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_PMK<- ENMevaluate(PMK_F[2:3], Pk_pres, algorithm = "maxent.jar",
                tune.args = tune.args, partitions = "block", overlap=F, 
                updateProgress = F, numCores = 4)
## look at the results 
Pk_result <- model_PMK@results

write.csv(Pk_result, "output/ENM/PMK/PMK_tune.csv", quote=FALSE, row.names=FALSE)

##Generating final models for P. keeni
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were linear Quadratic Hinge; RM = 5 
args=c("noaddsamplestobackground","noautofeature", "noproduct","nothreshold","betamultiplier=5.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#within SR
PK_model <- maxent(Pk_pres, PMK_F[2:3], args = args)
PK_predict <- predict(Pk_pres, PK_model, args=pred.args) 
plot(PK_predict)
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="black", cex=.35)

##projection
PK_model_pro <- maxent(Pk_pres_pro, PMK_F[2:3], args = args)
PK_predict_pro <- predict(Pk_pres_pro, PK_model_pro, args=pred.args) 
plot(PK_predict_pro)
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="black", cex=.35)



plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PK_predict_pro, add=T) 
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="black", cex=.4)

#load in LGM layers
#load in variables 
env_lgm <- list.files("data/env/21K/", pattern='tif', full.names=TRUE)
env_lgm
predictors_lgm <- stack(env_lgm)
projection(predictors_lgm) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
plot(predictors_lgm[[1]])
#Clipping layers to predict to a larger area 
layers_lgm <-  crop(predictors_lgm, PM_buf_pro[[1]])
env_lgm <- mask(layers_lgm, PM_buf_pro[[1]])
plot(env_lgm[[1]])

pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#load in ice layers
#load("data/GIS/IceSheets.RData")
#crop LGM layer to same extent 
#ice <-  crop(iceList[[43]], extent(lgm_env))

#LGM predictions
PM_lgm <- predict(env_lgm, PK_model , args=pred.args) 
plot(PM_lgm)

par(mfrow=c(1,2))
plot(PK_predict_pro)
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="black", cex=.5)
plot(PM_lgm)
#plot(ice, col="lightblue", add=T)

writeRaster(PK_predict,"output/ENM/PMK/PK_curr.asc", format="ascii", overwrite=T)
writeRaster(PK_predict_pro,"output/ENM/PMK/PK_curr_pro.asc", format="ascii", overwrite=T)
writeRaster(PM_lgm,"output/ENM/PMK/PK_lgm.asc", format="ascii", overwrite=T)

##########################################################################
##Peromyscus maniculatus sonoriensis modeling 
#clip layers by full extent 
layers_pms <- crop(WC, PM_buf[[2]])
Pms_pres <- mask(layers_pms, PM_buf[[2]])

plot(Pms_pres[[1]])
#plot(north_am, 1)
points(PMS$LONG, PMS$LAT, pch=16, col="red", cex=0.5)


#clip layers by full extent 
layers_pms_pro <- crop(WC, PM_buf_pro[[2]])
PMS_pres_pro <- mask(layers_pms_pro, PM_buf_pro[[2]])

plot(Pms_pres[[1]])
plot(PMS_pres_pro[[1]])

#spatial filter in spThin
PMS_thin <- thin(loc.data = PMS, 
                lat.col = "LAT", long.col = "LONG", 
                spec.col = "species", 
                thin.par = 25, reps = 100, 
                locs.thinned.list.return = TRUE, 
                write.files = TRUE, 
                max.files = 5,
                out.dir = "data/thin_PMS/", out.base = "PMS_f", 
                write.log.file = TRUE,
                log.file = "data/thin_PMS/PM_thin_log.txt")

#load in spatial thin dataset
PMS_F <- read.csv("data/thin_PMS/PMS_f_thin1.csv")

#plot 
plot(Pms_pres[[1]])
points(PMS$LONG, PMS$LAT, pch=16, col="black", cex=0.5)
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="red", cex=.5)

#lineage one model testing 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_PMS<- ENMevaluate(PMS_F[2:3], Pms_pres, algorithm = "maxent.jar",
                        tune.args = tune.args, partitions = "block", overlap=F, 
                        updateProgress = F, numCores = 4)
## look at the results 
PMS_result <- model_PMS@results

write.csv(PMS_result, "output/ENM/PMS/PMS_retune.csv", quote=FALSE, row.names=FALSE)

##Generating final models for P. maniculatus sor
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were Hinge; RM = 2
args=c("noaddsamplestobackground","noautofeature", "noproduct","linear", "quadratic", "nothreshold","betamultiplier=2.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#within SR
PMS_model <- maxent(Pms_pres, PMS_F[2:3], args = args)
PMS_predict <- predict(Pms_pres, PMS_model, args=pred.args) 
plot(PMS_predict)
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="black", cex=.35)

##projection
PMS_model_pro <- maxent(PMS_pres_pro, PMS_F[2:3], args = args)
PMS_predict_pro <- predict(PMS_pres_pro, PMS_model_pro, args=pred.args) 
plot(PMS_predict_pro)
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="black", cex=.35)


plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMS_predict_pro, add=T) 
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="black", cex=.4)

#Clipping layers to predict to a larger area 
layers_lgm <-  crop(predictors_lgm, PM_buf_pro[[2]])
env_lgm <- mask(layers_lgm, PM_buf_pro[[2]])
plot(env_lgm[[1]])

pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#load in ice layers
#load("GIS/IceSheets.RData")
#crop LGM layer to same extent 
#ice <-  crop(iceList[[43]], extent(lgm_env))

#LGM predictions
PM_lgm <- predict(env_lgm, PMS_model , args=pred.args) 
plot(PM_lgm)

PMS_lgm <- predict(env_lgm, PMS_model , args=pred.args) 
plot(PMS_lgm)

par(mfrow=c(1,2))
plot(PMS_predict_pro)
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="black", cex=.5)
plot(PMS_lgm)
#plot(ice, col="lightblue", add=T)

writeRaster(PMS_predict,"output/ENM/PMS/PMS_curr.asc", format="ascii", overwrite=TRUE)
writeRaster(PMS_predict_pro,"output/ENM/PMS/PMS_curr_pro.asc", format="ascii", overwrite=TRUE)
writeRaster(PMS_lgm,"output/ENM/PMS/PMS_lgm.asc", format="ascii", overwrite=TRUE)

##########################################################################
##Peromyscus maniculatus gambelii modeling 
#clip layers by study region (5 degree buffer)
layers_pmg <- crop(WC, PM_buf[[3]])
Pmg_pres <- mask(layers_pmg, PM_buf[[3]])

plot(Pmg_pres[[1]])
#plot(north_am, 1)
points(PMG$LONG, PMG$LAT, pch=16, col="red", cex=0.5)


#clip layers to projection extent (7.5 degree buffer)
layers_pmg_pro <- crop(WC, PM_buf_pro[[3]])
PMG_pres_pro <- mask(layers_pmg_pro, PM_buf_pro[[3]])

plot(Pmg_pres[[1]])
plot(PMG_pres_pro[[1]])

#spatial filter in spThin
PMG_thin <- thin(loc.data = PMG, 
                 lat.col = "LAT", long.col = "LONG", 
                 spec.col = "species", 
                 thin.par = 25, reps = 100, 
                 locs.thinned.list.return = TRUE, 
                 write.files = TRUE, 
                 max.files = 5,
                 out.dir = "data/thin_PMG/", out.base = "PMG_f", 
                 write.log.file = TRUE,
                 log.file = "data/thin_PMG/PM_thin_log.txt")

#load in spatial thin dataset
PMG_F <- read.csv("data/thin_PMG/PMG_f_thin1.csv")

#plot 
plot(Pmg_pres[[1]])
points(PMG$LONG, PMG$LAT, pch=16, col="black", cex=0.5)
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="red", cex=.5)

#lineage one model testing 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_PMG<- ENMevaluate(PMG_F[2:3], Pmg_pres, algorithm = "maxent.jar",
                        tune.args = tune.args, partitions = "block", overlap=F, 
                        updateProgress = F, numCores = 4)
## look at the results 
PMG_result <- model_PMG@results

write.csv(PMG_result, "output/ENM/PMG/PMG_retune.csv", quote=FALSE, row.names=FALSE)

##Generating final models for P. maniculatus gambelii
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were Hinge; RM = 6 
args=c("noaddsamplestobackground","noautofeature", "noproduct","linear", "quadratic","nothreshold","betamultiplier=6.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#within SR
PMG_model <- maxent(Pmg_pres, PMG_F[2:3], args = args)
PMG_predict <- predict(Pmg_pres, PMG_model, args=pred.args) 
plot(PMG_predict)
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="black", cex=.35)

##projection
PMG_model_pro <- maxent(PMG_pres_pro, PMG_F[2:3], args = args)
PMG_predict_pro <- predict(PMG_pres_pro, PMG_model_pro, args=pred.args) 
plot(PMG_predict_pro)
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="black", cex=.35)


plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMG_predict_pro, add=T) 
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="black", cex=.4)

#Clipping layers to predict to a larger area 
layers_lgm <-  crop(predictors_lgm, PM_buf_pro[[3]])
env_lgm <- mask(layers_lgm, PM_buf_pro[[3]])
plot(env_lgm[[1]])

pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#LGM predictions
PMG_lgm <- predict(env_lgm, PMG_model , args=pred.args) 
plot(PMG_lgm)

PMG_lgm <- predict(env_lgm, PMG_model , args=pred.args) 
plot(PMG_lgm)

par(mfrow=c(1,2))
plot(PMG_predict_pro)
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="black", cex=.5)
plot(PMG_lgm)
#plot(ice, col="lightblue", add=T)

writeRaster(PMG_predict,"output/ENM/PMG/PMG_curr.asc", format="ascii", overwrite=TRUE)
writeRaster(PMG_predict_pro,"output/ENM/PMG/PMG_curr_pro.asc", format="ascii", overwrite=TRUE)
writeRaster(PMG_lgm,"output/ENM/PMG/PMG_lgm.asc", format="ascii", overwrite=TRUE)


#load in all raster layers 
PK_curr <- raster("output/ENM/PMK/PK_curr_pro.asc")
PK_lgm <- raster("output/ENM/PMK/PK_lgm.asc")
PMS_curr <- raster("output/ENM/PMS/PMS_curr_pro.asc")
PMS_lgm <- raster("output/ENM/PMS/PMS_lgm.asc")
PMG_curr <- raster("output/ENM/PMG/PMG_curr_pro.asc")
PMG_lgm <- raster("output/ENM/PMG/PMG_lgm.asc")

range_expan <- read.csv("output/Range_expansion/PM_range_expansion.csv")
#plot ENMs 
pdf("Figures/PM_ENMs_F.pdf")
par(mfrow=c(3,2))
#Pk
plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PK_curr, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(PMK_F$LONG, PMK_F$LAT, pch=16, col="black", cex=.15)

plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PK_lgm, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(range_expan$longitude[[5]], range_expan$latitude[[5]], pch=16, col="red", cex=.25)

#PMS
plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMS_curr, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(PMS_F$LONG, PMS_F$LAT, pch=16, col="black", cex=.15)

plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMS_lgm, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(range_expan$longitude[[1]], range_expan$latitude[[1]], pch=16, col="yellow", cex=.25)
points(range_expan$longitude[[2]], range_expan$latitude[[2]], pch=16, col="yellow", cex=.25)
points(range_expan$longitude[[3]], range_expan$latitude[[3]], pch=16, col="yellow", cex=.25)

#PMG
plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-150, -50), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMG_curr, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(PMG_F$LONG, PMG_F$LAT, pch=16, col="black", cex=.15)

plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PMG_lgm, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), add=T) 
points(range_expan$longitude[[4]], range_expan$latitude[[4]], pch=16, col="blue", cex=.25)
dev.off()


