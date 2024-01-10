#########################################################
###                   AniMove 2023                    ###
### Script by Kami Safi, Martina Scacco & Anne Scharf ###
#########################################################
###                    Space USe                      ###
#########################################################

library(move2)
library(sf)
library(raster)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

# Load data
bats_mv2 <- mt_read("Parti-colored bat Safi Switzerland.csv")

#_________________________________
## Minimum Convex Polygon MCP ####
#_________________________________

library(adehabitatHR)
bats_mv2$id <- mt_track_id(bats_mv2)
bats_sp <- as_Spatial(bats_mv2[,'id'])

# function mcp is very particular about the input object, it must only
# contain 1 column, and names of the individuals names have to follow
# the validNames() rules
bats_sp <- bats_sp[,(names(bats_sp) %in% "id")]
levels(bats_sp$id) <- validNames(levels(bats_sp$id))

X330 <- bats_sp[bats_sp$id=="X330",]
# calculate mcp:
# by default, this package excludes the 5% locations farthest from the arithmetic mean
mcpX330 <- adehabitatHR::mcp(as(X330,"SpatialPoints"))
plot(X330, type="n", xlab="Longitude", ylab="Latitude")
plot(mcpX330, col="grey90", lty=2, lwd=1.25, add=TRUE)
points(X330, pch=16)
points(X330, pch=1, col="white")
legend("topright", as.character("95% MCP"), fill="grey90", bty="n")
mcpX330
# Note: area value seems strange. That is because our used locations are in the geographic coordinates system (long/lat). adehabitatHR calculated the area according to the units of the projection, in this case decimal degrees


# Therefore we have to project our data into a equidistant projection
#  We take again the median of the coordinates to center our projection
apply(st_coordinates(bats_mv2), 2, median)
AEQ_1p <- "+proj=aeqd +lat_0=47.22550 +lon_0=7.60154 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs"
habibaAEQD <- sf::st_transform(bats_mv, AEQ_1p)

library(rgeos)
#first option: reproject locations, than calculate mcp
bats_sp.proj <- spTransform(bats_sp, AEQ_1p)
mcpData.proj <- mcp(bats_sp.proj)
#second option: calculate mcp, than reproject mcp
mcpData <- mcp(bats_sp)
projection(mcpData) <- CRS("+proj=longlat +datum=WGS84") #declare original latlong proj before reprojecting
mcpData <- spTransform(mcpData, AEQ_1p)
plot(bats_sp.proj[bats_sp.proj$id=="X021",], bty="na", xlab="Longitude", ylab="Latitude")
plot(mcpData.proj[mcpData.proj$id=="X021",], add=TRUE)
plot(mcpData[mcpData$id=="X021",], add=TRUE, lty=2)
legend("bottomleft", c("First reproject then mcp", "First mcp then reproject"), lty=c(1,2), bty="n")
legend("topleft", sprintf("Area = %.2f", c(rgeos::gArea(mcpData.proj, byid=TRUE)["X021"],
                                           rgeos::gArea(mcpData, byid=TRUE)["X021"])/1000^2), lty=c(1,2), bty="n")
# Note to plot:
# - the 2 options result in different area calculations
# - adehabitatHR uses distance between locations to do the calculation
# - therefore always: FIRST project locations, and THEN calculate MCP


## Size of MCP changes with sampling effort or sampling size
library(move)
bats_mv <- to_move(bats_mv2)
hrBootstrap(bats_mv[['X021']], rep=500, levelMax=95)
legend("bottomright", legend=c("MCP size with all data","100% percentil","75% percentil","50% percentil","25% percentil","0% percentil"), 
       lty=c(4,2,3,1,3,2), col=c("black","cyan","red","black","red","cyan"), bty="n")
# Note to plot:
# - if sampling is large enough a saturation between sample size and area is reached

# In conclusion:
# simple and intuitive
# but does not account for intensity of utilization
# assumes independence of locations


#____________
# Kernel ####
#____________

# Kernel density estimations represent highly utilized locations with an increased "utilization value"
# Using a smoothing factor they model a 2 dimensional probability landscape of finding an animal at a certain location
# But has big assumption of observations being independent, which is violated here
# they are aimed at quantifying range distribution, but given this assumption, 
# with highly sampled data they measure occurrence distribution
# the higher the sampling the tighter data are hugged, so we don't see a home range but just the track

library(raster)
library(sf)

# first we just rasterize the data (n. of data point pers cell/pixels)
bats_mv2.proj <- sf::st_transform(bats_mv2, AEQ_1p)
ind191 <- filter_track_data(bats_mv2.proj, .track_id="191")
template <- raster(extent(sf::st_bbox(ind191)))
res(template)<-500
count <- rasterize(as_Spatial(ind191), template,field=1,  fun="count")
plot(count, col=grey(10:0/12))
plot(mcpData.proj[mcpData.proj$id=="X191",], add=TRUE)
points(as_Spatial(ind191), pch=16, cex=0.5)
# Note: the result is highly dependent on the chosen cell size. In this case 500x500m


# kernel implementation by "adehabitatHR" library
library(adehabitatHR)
library(scales)
X021 <- bats_sp.proj[bats_sp.proj$id=='X021',]
kern1 <- kernelUD(as(X021, "SpatialPoints"), h=500)
kern2 <- kernelUD(as(X021, "SpatialPoints"))
kern3 <- kernelUD(as(X021, "SpatialPoints"), h=2000)
kern4 <- kernelUD(as(X021, "SpatialPoints"), h="LSCV")
par(mfrow=c(2,2))
par(mar=c(1,0.5,3,0.5))
kern <- c("kern1", "kern2", "kern3", "kern4")
hName <- c("h=500",
           "h='ad-hoc'",
           "h=2000",
           "h=LSCV")
for(i in 1:4){
  plot(getverticeshr(get(kern[i])))
  points(X021, pch=16, cex=0.75, col=alpha("black", 0.2))
  points(X021, cex=0.75)
  title(hName[i])
}
# Note to plot:
# - h: degree of smoothness or how tightly the data should be hugged by the distribution function
# - h="LSCV": h calculated from the data via least square cross validation
# - h="ad-hoc": h calculated from the data via sample size and spatial spread

# kernel implementation by "ks" library
library(ks)
library(scales)
pos <- coordinates(X021)
H.bcv <- Hbcv(x=pos)
H.pi <- Hpi(x=pos)
H.lscv <- Hlscv(x=pos)
H.scv <- Hscv(x=pos)

par(mfrow=c(2,2))
par(mar=c(1,0.5,3,0.5))
H <- c("H.bcv", "H.pi", "H.lscv", "H.scv")
hT <- c("Biased cross-validation (BCV)",
        "Plug-in",
        "Least-squares cross-validation",
        "Smoothed cross-validation")
for(i in 1:4){
  fhat <- kde(x=pos, H=get(H[i]), compute.cont=TRUE)
  plot(fhat, cont=c(75, 50, 5), bty="n",
       xaxt="n", yaxt="n",
       xlab=NA, ylab=NA, asp=1,display="filled.contour")
  points(X021, pch=16, cex=0.75, col=alpha("black", 0.2))
  title(hT[i])
}

# In conclusion, to think about:
# - Choice of h: the goal is to find a balance between fitting the data with a minimum 
# number of parameters in order to generalise space use, while avoiding over-fitting
# - The underlying kernel function is assumed to be a Gaussian bi-variate distribution (around each location) 
# thought to represent error of the true location and a diffusion process when no position information was obtained
# - assumes independence of locations, with underestimating space use with increasing samp freq
# - assumes regular sampling

# Note: 
# On Friday you will learn about a new kernel estimator that can incorporate 
# temporal autocorrelation pattern and thereby provide unbiased range estimation (Fleming et al. 2015)

#____________
## LoCoH ####
#____________

# They allow to better respect hard natural boundaries (such as rivers, roads, and lakes or transitions in vegetation like forest edges)
# in the calculation of spece use
# The idea behind Local Convex Hulls is to union a series of locally built 
# minimum convex hulls for each location based on:
# - a fixed number of neighbours
# - the neighbours within a radius from each location
# - the neighbours whose sum of distances is within a certain "a" value
# Easier to minimize the inclusion of non-utilized space in a simple manner

# check vignettes of LoCoH.k() LoCoH.r() LoCoH.a() from package adehabitatHR
library(move2)
library(adehabitatHR)
library(sf)
library(imager)
#library(maptools)

# Takes into account the physical location and their correlation in space
# Exclude areas that are not used
# The choice of k, r, or a strongly affects the results

fishers <- mt_read(mt_example())
leroy <- filter_track_data(fishers, .track_id="M1")
leroy <- leroy[!sf::st_is_empty(leroy),] ## the trasformation to "SpatialPoints" does not allow empty locations

# data need to be transformed into a equidistant projection because method relays on distance calculations
leroy_proj <- st_transform(leroy, crs="EPSG:32618")

leroy_sp <- as(as_Spatial(leroy_proj),"SpatialPoints")
# png(filename="locoh_maps.png", width = 21, height = 15, units = "cm", res=300)
par(list(mfrow=c(2,2), mar=c(2,2,2,2)))
leroy.mcp <- mcp(leroy_sp, percent=95)
plot(leroy_sp, col=grey(0.9), lty=2, lwd=2)
points(leroy_sp, col="#00000060", pch=16, cex=0.5)
lines(coordinates(leroy_sp), col="#00000030")
title("Minimum convex polygon")
# include "k" number of closest neighbour locations
kLoc <- LoCoH.k(leroy_sp, k=75)
plot(kLoc, col=grey((0:length(kLoc)/length(kLoc))*0.7), border=NA)
title("k-NNCH LoCoH")
# include location within a radius "r"
rLoc <- LoCoH.r(leroy_sp, r=800)
plot(rLoc, col=grey((0:length(rLoc)/length(rLoc))*0.7), border=NA)
title("r-NNCH LoCoH")
# sum of distances of included neighbour locations to root location is "a"
aLoc <- LoCoH.a(leroy_sp, a=9000)
plot(aLoc, col=grey((0:length(aLoc)/length(aLoc))*0.7), border=NA)
title("a-NNCH LoCoH")
# dev.off()


x11();plot(imager::load.image("locoh_maps.png"),axes=F)

## area changes depending on the choice of k, r, or a
dev.off()# to reset the plotting environment
# png(filename="locoh_Area.png", width = 21, height = 15, units = "cm", res=300)
par(mfrow=c(1,3))
kLocArea <- LoCoH.k.area(leroy_sp,
                         krange=floor(seq(75, 500, length=10)),
                         percent=90)
title("k-NNCH LoCoH")
rLocArea <- LoCoH.r.area(leroy_sp,
                         rrange=seq(500, 1600, 100),
                         percent=90)
title("r-NNCH LoCoH")
aLocArea <- LoCoH.a.area(leroy_sp,
                         arange=seq(5000, 13000, 1000),
                         percent=90)
title("a-NNCH LoCoH")
# dev.off()

x11();plot(load.image("locoh_Area.png"),axes=F)

#________________
# dBBMM & UD ####
#________________

library(move)

leroy_mv <- to_move(leroy)
leroy_prj_ctr <- spTransform(leroy_mv, center=TRUE)

# The smaller the window size, the more subtle changes in var will be considered as significant. 
# The larger the margins within the windows, the more stable are the breaks that are detected, 
# at the cost of not being able to model utilisation distribution at the start and end of the trajectory the number of locations indicated in the margin.size.

# argument ext gives "room" around the extent of the data for the uncertainty. 

# Check the timeLag of your data!! 
# If there are timelags shorter than intended the calculation could take very long,
# you can use the argument "timestep" in the dBBMM function (see ?brownian.bridge.dyn), to make sure calculation does not take forever
# On the opposite if there are some very big timelags, the uncertainty becomes very big, and the raster might not be big enough to contain the data+uncertainty
# In this case we can ignore points with very big uncertainty (for this we need to estimate variance separately as we did yesterday)

# Raster granularity=number of cells used to divide the largest dimension of the track into
# Granularity of raster affects smoothness of the probability surface as well as strongly affect computation time
summary(timeLag(leroy_prj_ctr,"mins"))
BB.leroy <- brownian.bridge.dyn(leroy_prj_ctr, ext=.45, dimSize=150, 
                                location.error=20, margin=11, window.size=31)
plot(BB.leroy)
# Note to plot:
# - all pixels sum up to 1 (total time tracked)
# - the values correspond to the proportion of the time tracked that the animal spend in that area
# - in this case, a pixel with value 0.035 means that leroy spend 3.5% from the total tracking time in that pixel
# - these results are very useful to find out where the animal spend how much time, and e.g. relate it to environmental variables


## extract the utilization distribution (UD)
udleroy <- getVolumeUD(BB.leroy)
plot(udleroy, col=terrain.colors(100))

## from the ud object, also the contours can be extracted
plot(leroy_prj_ctr, col="#00000060", pch=16, cex=0.5, bty="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)
lines(leroy_prj_ctr, col="#00000030")
contour(udleroy, levels=c(0.5, 0.95), add=TRUE, lwd=c(2, 1), lty=c(2,1))
title("Dynamic brownian bridge")

## plotting the UD95 on google map
library(ggmap)
library(ggplot2)
cl95 <- raster2contour(BB.leroy, levels=0.95)
cl95LL <- spTransform(cl95, CRS("+proj=longlat"))
cl95df <- data.frame(do.call("rbind", coordinates(cl95LL)[[1]]))
cl95df$polyN <- rep(1:length(coordinates(cl95LL)[[1]]), lapply(coordinates(cl95LL)[[1]], nrow))
leroyDF <- as.data.frame(spTransform(leroy_prj_ctr,"+proj=longlat"))
m <- get_map(sp::bbox(extent(cl95LL)*1.5), zoom=13, source="stamen", maptype="terrain")
ggmap(m)+geom_path(data=cl95df, aes(x=X1,y=X2,group=polyN),color="red")+
  geom_path(data=leroyDF, aes(x=coords.x1, y=coords.x2),alpha=0.2)+
  geom_point(data=leroyDF, aes(x=coords.x1, y=coords.x2),alpha=0.3, shape=20)+
  labs(x="",y="")+
  theme(axis.text=element_blank(),axis.ticks=element_blank())+
  theme(legend.position="none")


## checking for effect of margin and window size on the results
par(mfrow=c(3,3), mar=c(1,2,2,1))
margins <- c(15, 9, 3)
windows <- c(101, 67, 33)
runs <- expand.grid(margins, windows)
for(i in 1:nrow(runs)){
  BB.leroy <- brownian.bridge.dyn(leroy_prj_ctr, dimSize=150, location.error=20, margin=runs[i,1],
                                  window.size=runs[i,2], time.step=2, ext=2)
  udleroy <- getVolumeUD(BB.leroy)
  udleroy99 <- udleroy
  udleroy99[udleroy99>0.99] <- NA
  udleroy99 <- trim(udleroy99,padding=5)
  contour(udleroy99, levels=c(0.5, 0.95), bty="n", xaxt="n",
          yaxt="n", xlab=NA, ylab=NA, asp=1)
  mtext(paste("Margin = ", runs[i,1], sep=""), 2)
  mtext(paste("Window size = ", runs[i,2]), 3)
}
# Note to plot: no need to worry all to much about the window size and the margin, as they do not have a major impact on the results. Default values work well.


#________________________________________________________
# Effect of sampling frequency and tracking duration ####
#________________________________________________________

### Effect of sampling frequency and tracking duration on area

library(adehabitatLT)

# Here we simulate a track and play with thinning and shortening overall duration 
# and check the effect of on the calculated kernel UD
dev.off()
par(mfrow=c(1,1))
par(list(mar=c(5, 4, 4, 2) + 0.1), bty="o")
set.seed(3628492)
steps <- 100000
prop=seq(0.01,1,0.01)

track <- as(simm.crw(date=1:steps, h=1, r=0.8), "Move")
thin <- lapply(prop, function(x) track[round(seq(1, steps, length.out=steps * x)), ])
short <- lapply(prop, function(x) track[1:round(steps * x),])
ThinAreas <- lapply(lapply(lapply(thin, as, "SpatialPoints"), kernelUD), kernel.area, percent=95)
ShortAreas <- lapply(lapply(lapply(short, as, "SpatialPoints"), kernelUD), kernel.area, percent=95)

plot(I(unlist(ThinAreas)/min(unlist(ThinAreas)))~seq(1:100),
     xlab="Percent of the track",
     ylab="Relative area", ylim=c(0,1.75),
     type="l", lwd=2, lty=1)
lines(I(unlist(ShortAreas)/max(unlist(ShortAreas)))~seq(1:100), lty=2, lwd=2)
abline(h=1)
legend("topright", c("Thinned trajectory", "Shortened trajectory"), lty=c(1,2), lwd=c(2,2), bty="n")
# Note to plot:
# - in this case using the kernelUD, when sampling frequency is lower, i.e. thinned trajectory, than the estimated UD is larger than when the sampling frequency is higher
# - the longer the trajectory, the larger the UD


