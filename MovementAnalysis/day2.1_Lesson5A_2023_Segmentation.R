#########################################################
###                   AniMove 2023                    ###
### Script by Kami Safi, Martina Scacco & Anne Scharf ###
#########################################################
###             Behavioural segmentation              ###
#########################################################

# link to class recording from Animove 2022, Chloe Bracis
# https://animove.org/elearning/ -> Segmentation

library(move)
library(move2)
library(ggplot2)
library(plotly)
library(scales)
library(units)
library(lubridate)
library(sf)
library(rnaturalearth)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

parDefaults <- par()

#______________________________
### Load data and general check

# load Sierit the stork
Sierit <- mt_read("Sierit (DER AN858).csv")
Sierit <- Sierit[!sf::st_is_empty(Sierit),]

# thinning the data just to make the example code run faster!!!
Sierit <- mt_filter_per_interval(Sierit, criterion="first", unit="hour")

# first inspect the data
bb <- st_bbox(Sierit)
ggplot() + theme_bw() +
  geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", 50)) +
  geom_sf(data = mt_track_lines(Sierit), color="red") +
  geom_sf(data = Sierit, color=alpha("red",0.3)) +
  coord_sf(xlim=c(bb[1], bb[3]), ylim=c(bb[2], bb[4]))

# check for gaps and irregularities in data
plot(mt_time(Sierit), set_units(mt_time_lags(Sierit), "hour"), 
     ylab="time lag in hours", xlab="")
plot(mt_time(Sierit), set_units(mt_time_lags(Sierit), "hour"), 
     ylab="time lag in hours", xlab="", ylim=c(0,10))

# check for the daily schedule (in this case the tag is turned off over night)
sub48h <- Sierit[date(mt_time(Sierit))%in%date(c("2013-06-29","2013-06-30")),]
plot(mt_time(sub48h),set_units(mt_time_lags(sub48h), "hour"), ylab="time lag in hours",xlab="")


#_____________
# K-MEANS ####
#_____________

# Clusters data into k clusters of choice based on their euclidean distance
# can be used on any variable or combinations of variable (multivariate)

#______________________
### Spatial clustering

range(mt_time(Sierit))
Sierit_sub <- filter(Sierit, format(mt_time(Sierit), "%Y-%m") %in% c("2013-06","2013-07"))
range(mt_time(Sierit_sub))

bb <- st_bbox(Sierit_sub)
ggplot() + theme_bw() +
  geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", 50)) +
  geom_sf(data = Sierit_sub, color=alpha("red",0.3)) +
  coord_sf(xlim=c(bb[1], bb[3]), ylim=c(bb[2], bb[4]))

# project in metric coordinates and apply kmeans with 3 clusters
apply(st_coordinates(Sierit_sub), 2, median)
AEQ_1p <- "+proj=aeqd +lat_0=47.75148 +lon_0=8.93148 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs"
Sierit_subProj <- st_transform(Sierit_sub, crs = AEQ_1p)

cl <- kmeans(st_coordinates(Sierit_subProj), centers = 3)
str(cl)

Sierit_subProj$SPclust <- as.character(cl$cluster)
ggplot() + theme_bw() +
  geom_sf(data = Sierit_subProj, aes(color=SPclust))


#________________________
### Space-time clustering

# Multiply time (s) by speed (m/s) to get a time in units of metres, compatible with the x-y locations
time_m <- set_units(as.numeric(mt_time(Sierit_subProj)), s) * set_units(mt_speed(Sierit_subProj), m/s)
spacetime_df <- cbind(time_m, st_coordinates(Sierit_subProj))[-length(time_m),]

clST <- kmeans(spacetime_df, centers = 3)

# plot clustering in 2D
ggplot(as.data.frame(spacetime_df), aes(X, Y, color=as.character(clST$cluster))) + 
  geom_point() +
  theme_bw()
# or 3D
plot_ly(x=spacetime_df[,"X"], y=spacetime_df[,"Y"], z=spacetime_df[,"time_m"], 
        type="scatter3d", mode="markers", 
        color=as.character(clST$cluster))


#______________
# LAVIELLE ####
#______________

# The Lavielle method (Lavielle 1999, 2005) uses a penalized contrast paradigm to detect changes / find potential breaks in any time series.
# The method tries to minimize the contrast between a segmented time series and the original time series 
# given a number of pre-defined segments by comparing either mean, variance or both measures of a time series between segments.
# We show two examples of Lavielle applied to two metrics you are now familiar with: FPT and NSD
# adehabitatLT provides the function lavielle()
# it has a lot of defaults, make sure you check the possible arguments
?lavielle

library(adehabitatLT)

#______________________________________________
### Applied to Net Square Displacement NSD ####

## for migrating birds for which we have the nest location, NSD often gives a nice overview of what is happening
# subsetting the data to one position a day
Sierit1day <- mt_filter_per_interval(Sierit, criterion="first",unit="day")
SieritNSD <- set_units((st_distance(Sierit1day,Sierit1day[1,]))^2," km^2")

plot(mt_time(Sierit1day), SieritNSD, type="l",xlab="Time", ylab="Net square distance (Km²)", main="Sierit")
points(mt_time(Sierit1day), SieritNSD, pch=20, cex=0.7)

nsd_l <- lavielle(drop_units(SieritNSD), Lmin=4, Kmax=15, type="var")
(nsegms <- chooseseg(nsd_l, output = "opt", draw=F))


par(parDefaults)
par(mfrow=c(1,2))
(fp <- findpath(nsd_l, nsegms))

sieritCoords <- st_coordinates(Sierit1day)
cols <- rainbow(length(fp))
plot(sieritCoords[,"X"], sieritCoords[,"Y"], 
     xlab="Longitude", ylab="Latitude", type="l", col=alpha("black", 0.8))
lapply(1:length(fp), function(i){
  lines(sieritCoords[fp[[i]][1]:fp[[i]][2], "X"],
        sieritCoords[fp[[i]][1]:fp[[i]][2], "Y"], col=cols[i])})
legend("bottomright", paste("Segment ",c(1:length(fp))), col=cols, lwd=2, bty="n", cex=0.5)


#__________________________________________
###  Applied to First Passage Time FPT ####

library(waddle)
library(adehabitatLT)

# project the data as this and other methods methods assume equidistant projection
# since it is a migration track we use again the two-points equidistant projection we saw yesterday
st_bbox(Sierit)
AEQ_2p <- "+proj=tpeqd +lon_1=-6.833576 +lat_1=33.826017 +lon_2=9.027470 +lat_2=47.761884 +datum=WGS84"
Sierit.prj <- st_transform(Sierit, crs=AEQ_2p)

# create "as.ltraj" class object
Serit.ltraj <- as.ltraj(xy=st_coordinates(Sierit.prj), date=mt_time(Sierit.prj),id=mt_track_id(Sierit.prj), typeII=T)
# either give predefined radius, or check where there may be changes in behavioural state
ftpradii <- 10^seq(2, 6, length.out=50)
fptSieritallR <- fpt(Serit.ltraj, radii=ftpradii, units="days")
# checking the variance of log(ftp)
par(mfrow=c(1,1))
vars <- varlogfpt(fptSieritallR, graph=F)
plot(as.numeric(vars)~ftpradii,
     type="l", lwd=1, lty=2,
     log="x", ylab="Variance of log first passage time",
     xlab="Radius in meters")

# getting radii of the max and min variances
radiiSelec <- which(diff(floor(diff(as.numeric(vars))))==-1)+1
round(ftpradii[radiiSelec])
# calculate FPT, for the selected radii
radii <- round(ftpradii[radiiSelec])
Sierit.fpt <- fpt(Serit.ltraj, radii, units="days")
#look at the results
x11();par(mfrow = c(3,2))
plot.fpt(Sierit.fpt, radii, xlab = "Time", main ="Radius: ")

# calculate the breakpoints using lavielle()
# have to adjust Lmin and Kmax
Sierit.lav  <-  apply(Sierit.fpt[[1]], 2, function(x) lavielle(na.omit(x), Lmin = 50, Kmax = 20))
# plot represent the estimation of the optimal number of segments per radius
x11();par(mfrow = c(3,2))
Sierit.optseg  <-  unlist(lapply(Sierit.lav, chooseseg, output = "opt", draw=F))

# plot the segments
x11();par(mfrow = c(3,2))
temp <-  lapply(1:length(radii), function(x){
  findpath(Sierit.lav[[x]], Sierit.optseg[x])
  title(paste("Radius: ", radii[x])) })



#_______________________________________________________
# EMbC - Expectation Maximization Binary Clustering ####
#_______________________________________________________

# This method splits time series observations into binary categories finding a best split based on machine learning processes. 
# The so-called expectation maximization (EM), is based on two metrics (usually speed and turning angle).
# High speed, low turn angle = directed movement; Low speed high turn angle = foraging
# Other variables can be used, and solar height can be added (see vignette)

## see vignette of EMbC package for details
library(EMbC)
# create data.frame as needed by Embc
dfSierit <- data.frame(timeStamp=mt_time(Sierit),
                       lon=st_coordinates(Sierit)[,1],
                       lat=st_coordinates(Sierit))
# speed/turning angle clustering. Other variables can be used for clustering (see vignette)
bc <- stbc(dfSierit,smth=24)
stts(bc) #get the statistics of the two variables per cluster
# obtain 4 categories, resulting from the combinations of high/low speed/turning angle
par(mfrow=c(1,1))
sctr(bc)

Sierit <- mutate(Sierit, speed=drop_units(mt_speed(Sierit)), turnAngle=drop_units(mt_turnangle(Sierit)))
SieritClass <- filter(Sierit, complete.cases(speed,turnAngle))
bc2 <- embc(as.matrix(cbind(sqrt(SieritClass$speed), abs(SieritClass$turnAngle))))
sctr(bc2)
# the upper line on the plot reflects to the time spend in each category

## lets make a nicer plot
# reclassifing categories into high and low speed
Cat <- bc2@A
Cat[Cat == 1] <- "LL" # low speed + low turn angle
Cat[Cat == 2] <- "LH" # low speed high turn angle
Cat[Cat == 3] <- "HL" 
Cat[Cat == 4] <- "HH" 
SieritClass$embc_cat <- as.factor(Cat)
SieritClass$segments <- mt_segments(SieritClass)

ggplotly(ggplot(SieritClass)+
  theme_void() +
  geom_sf(aes(geometry=segments, color=embc_cat)))


#______________________________________________
# BCPA - Behavioural Change Point Analysis ####
#______________________________________________

# The behavioural change point analysis (BCPA: Gurarie et al. (2009))is a parametric method that uses a sliding window 
# and likelihood estimates to decide whether a metric of the trajectory, usually persistence velocity, has possibly changed. 
# Uses the autocorrelation in velocity as a feature rather than getting rid of it
# Same rationale behind methods applied previously: finding break points in a time series.
# The smaller the window, the shorter the switches in the detected changes can become, at the cost of finding spurious changes. 
# The larger the window, the more reliable the detected changes become, while missing out on the short term changes.

library(bcpa)
# also check out the package "smoove" at https://github.com/EliGurarie/smoove
## see detailed vignette from bcpa R package
# create object needed for the WindowSweep function, input needs to be a object class "track"
SieritT <- MakeTrack(X=st_coordinates(Sierit.prj)[,1], Y=st_coordinates(Sierit.prj)[,2], Time=mt_time(Sierit.prj))
# bcpa::GetVT table computes speeds, step lengths, orientations, turning angles
Sierit.VT <- bcpa::GetVT(SieritT, units = "min")

# set the 3 necessary parameters
?WindowSweep
windowsize <- 40 # number of locations within window (not time!). min 20. Larger windows are more robust but more coarse, smaller windows are more sensitive but more likely to give slightly spurious results, compensated by adjusting the K
windowstep <- 1 # if very high resolution data use value greater then 1
K <- 0.5 # default 2, with a lower value, a simpler model is more likely to be selected

Sierit.ws <- WindowSweep(Sierit.VT, "V*cos(Theta)",
                         time.var = "T.mid",
                         windowsize, windowstep,
                         plotme = FALSE, K = K,
                         tau = TRUE, progress = FALSE)

# plotting the results for smooth and flat BCPA
x11();par(mfrow = c(2,1))
plot.bcpa(Sierit.ws, type="smooth", threshold = 10, legend = FALSE)
plot.bcpa(Sierit.ws, type="flat", clusterwidth = 5, legend = FALSE)
# threshold: indicates how many of the windows that were swept over that data point 
# must have selected it as a significant change point
# clusterwidth (only in flat):temporal range within which different change points are still considered to be within the same cluster

# get dignostic plots
DiagPlot(Sierit.ws, "smooth") #looks at each point, many behavioural states
DiagPlot(Sierit.ws, "flat") #looks at more discrete states

# plot the segmented trajectory
PathPlot(SieritT, Sierit.ws, type="smooth", clusterwidth=10, main="Smooth BCPA")
PathPlot(SieritT, Sierit.ws, type="flat", clusterwidth=5, main="Flat BCPA")


#___________________________________________________
# BPMM - Bayesian Partitioning of Markov Models ####
#___________________________________________________
# Assumptions:
# 1. Movement composed of discrete number of homogeneous models
# 2. Markovian transition between the models (no long-term memory)
# 3. movement variable is independent within state
# Steps
# 1. regularize data
# 2. check autocorrelation
# 3. select movement variable (e.g. step length)
# 4. define set of candidate models, and select optimal number
# 5. partition track

library(bcpa)
# First check for autocorrelation in speed
distSierit <- drop_units(mt_distance(Sierit.prj)) %>% na.exclude()
cor(distSierit[-length(distSierit)], distSierit[-1])
# thin and recalculate
Sierit.prj.thin <- mt_filter_per_interval(Sierit.prj, criterion="first", unit="8 hours")
distSierit <- drop_units(mt_distance(Sierit.prj.thin)) %>% na.exclude()
cor(distSierit[-length(distSierit)], distSierit[-1])
# now transform into a track object
SieritT <- MakeTrack(X=st_coordinates(Sierit.prj.thin)[,1], Y=st_coordinates(Sierit.prj.thin)[,2], Time=mt_time(Sierit.prj.thin))
# and regularize the track (fill gaps if needed)
SieritT.reg <- waddle::InterpolatePoints(SieritT, 8, "hour")$Data
# convert it into a object of class "ltraj"
SieritT.traj <-  as.ltraj(data.frame(SieritT.reg$X, SieritT.reg$Y), SieritT.reg$Time, id = "Sierit")
# plot step length over time
plot(SieritT.traj[[1]]$date, SieritT.traj[[1]]$dist, xlab="", ylab="Distance (m)", type="l", col="darkgrey")
points(SieritT.traj[[1]]$date, SieritT.traj[[1]]$dist, pch=19, cex=0.5)

# check for sd of step length
sd(SieritT.traj[[1]]$dist[500:550],na.rm=T)
sd(SieritT.traj[[1]]$dist,na.rm=T)
# estimation of the optimal number of partitions, see also modpartltraj in adehabitatLT
?Prep.segments
# arguments that have to be set and adjusted: sd, nmodels,Km,log
SieritT.segments <- Prep.segments(SieritT.traj, units = "min", sd = 1800, nmodels=20, Km=50, log=T)

# perform the BPMM partitioning
SieritT.partition <- Partition.segments(SieritT.segments)
plot.segments(SieritT.partition, xlab="", ylab="Step length (m)")

# get diagnostic plots
DiagPlot.segments(SieritT.partition)
# QQ plot looks ok
# Histogram should represent a normal distribution
# autocorrelation is still present


#___________________________________
# MRW - Multistrate Random Walk ####
#___________________________________
# Assumptions:
# 1. Animal transitions between several discrete states
# 2. Each state is associated with unique parameters
# Steps:
# 1. regularize data
# 2. check for autocorrelation
# 3. define distributions of step length and turning angle
# 4. fit the model
# 5. evaluate the model

library(moveHMM)
# regularize track and thin to remove autocorrelation. In this case we use the object in longitude/latitude
Sierit.thin <- mt_filter_per_interval(Sierit, criterion="first", unit="8 hours")
SieritTLL <- MakeTrack(X=st_coordinates(Sierit.thin)[,1], Y=st_coordinates(Sierit.thin)[,2], Time=mt_time(Sierit.thin))
SieritTLL.reg <-  InterpolatePoints(SieritTLL, n=8, units="hour")$Data
# prepare data
Sierit.hmm <- prepData(SieritTLL.reg, type = "LL", coordNames = c("X", "Y"))
# plot the step length and turning angle distributions
plot(Sierit.hmm, compact = T)
summary(Sierit.hmm)

# give distribution of step length and turning angle based on the previous plots. 
# Hint: google weibull,gamma,etc and select shape and scale that best fits the step length distribution. 
# Same for wcauchy mu and rho for turning angle distribution
## Two-state model
weibullShape <- c(3, 0.5)
weibullScale <- c(10,10)
wcauchyMu <- c(0, pi)
wcauchyRho <- c(0.7, 0.2)

doubleStateModel <- fitHMM(Sierit.hmm, nbStates = 2,
                           stepPar0 = c(weibullShape, weibullScale),
                           anglePar0 = c(wcauchyMu, wcauchyRho),
                           stepDist = "weibull", angleDist = "wrpcauchy")
# output
doubleStateModel
# confidence intervals
CI(doubleStateModel)
# plots of fitted distributions
par(mfrow=c(1,2))
plot(doubleStateModel,ask=F)
# plots of segments
plotStates(doubleStateModel,ask=F)

# Three-state model
gammaShape <- c(3, 1, 1)
gammaScale <- c(20, 5, 2)
wcauchyMu <- c(0, pi, 0)
wcauchyRho <- c(0.9, 0.2, 0.6)

tripleStateModel <- fitHMM(Sierit.hmm, nbStates = 3,
                           stepPar0 = c(gammaShape, gammaScale),
                           anglePar0 = c(wcauchyMu, wcauchyRho),
                           stepDist = "gamma", angleDist = "wrpcauchy")

par(mfrow=c(1,2))
plot(tripleStateModel,ask=F)
x11();par(mfrow=c(1,1))
plotStates(tripleStateModel,ask=F)

## check which model has the lowest AIC
AIC(doubleStateModel, tripleStateModel)


#_______________
# CORRIDORS ####
#_______________

# The function corridor identifies pinch points in the movement characterised 
# by parallel movement (LaPoint et al. 2013).
# Finding corridors makes most sense for trajectories where the animal is 
# predominantly territorial/within home range.

library(move)
# identify corridor behaviour, i.e. parallel fast movement. For details see ?corridor
Leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
LeroyCorr <- corridor(Leroy)
x11();plot(LeroyCorr, type="l", xlab="Longitude", ylab="Latitude", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")






