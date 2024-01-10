#########################################################
###                   AniMove 2023                    ###
### Script by Kami Safi, Martina Scacco & Anne Scharf ###
#########################################################
###               Trajectory Analysis                 ###
#########################################################

library(move2)
library(move)
library(adehabitatLT)
library(sf)
library(lubridate)
library(circular)
library(ggplot2)
library(rnaturalearth)
library(units)
library(dplyr)
library(purrr)
library(tidyr)
library(scales)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

# We will look at simple geometric properties of the tracks and 
# how movement processes are affected by small changes in these properties

parDefaults <- par()

#______________________________________
# AZIMUTH, TURNING ANGLE and SPEED ####
#______________________________________

# load and plot Leo
Leo <- mt_read("Leo-65545.csv")

ggplot() + theme_bw() +
  geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", 50)) +
  geom_sf(data = mt_track_lines(Leo), color="firebrick")+
  coord_sf(xlim=c(min(sf::st_coordinates(Leo)[,1]),max(sf::st_coordinates(Leo)[,1])), ylim=c(min(sf::st_coordinates(Leo)[,2]),max(sf::st_coordinates(Leo)[,2])))

# get locations per year and month
tapply(mt_time(Leo), list(year(mt_time(Leo)), month(mt_time(Leo))), length)
table(year(mt_time(Leo)), month(mt_time(Leo)))

# removing years 2012-13 because there is a large gap
Leo <- Leo[year(mt_time(Leo)) < 2012,]

# Annotate speed, azimuth and turning angle to the trajectory.
Leo <- Leo %>% mutate(azimuth = mt_azimuth(Leo), 
                      speed = mt_speed(Leo), 
                      turnangle = mt_turnangle(Leo))
# azimuth and turnangle take values between -pi and pi
# in mt_azimuth positive is east
# in mt_turnangle positive is right turn

###________________________________________________
### Azimuth distribution with speed and season ####

### Categorize seasons. Numbers refer to month.
# assign the categories to a new variable based on the timestamp
# assign season = NA to segments where subsequent locations fall in different seasons
Leo  <- mutate(Leo,
               month = month(mt_time(Leo)),
               season = recode_factor(month,
                                      "1"="Wintering", "2"="Wintering",
                                      "3"="Wintering", "4"="North migration",
                                      "5"="North migration", "6"="Breeding",
                                      "7"="Breeding", "8"="Breeding",
                                      "9"="South migration", "10"="South migration",
                                      "11"="Wintering", "12"="Wintering"),
               season = if_else(season == lead(season, 1) &
                                  mt_track_id(Leo) == lead(mt_track_id(Leo), 1),
                                season, NA)
)

### Look at direction of movement when movement is happening, that is:
# select segments above 2 m/s, we are only interested in segments when Leo is moving, and not the stationary error
# remove missing values in season and group data by season
Leo_moving <- filter(Leo, speed > set_units(2, "m/s") & !is.na(season)) %>% group_by(season)


### Speed ~ azimuth scatter plot
# speed here is ground speed (how fast is the animal in making way)
ggplot(Leo_moving, aes(x = azimuth, y = speed)) +
  geom_point(color=alpha("black", 0.3)) +
  scale_x_units(unit = "degrees", breaks = c(-2:2) * 90, expand = c(0L, 0L)) +
  scale_y_units(unit = "m/s") +
  theme_linedraw()
# Note to plot:
# -high speeds predominantly occur when azimuths are around roughly -180-180∘ and 0∘ (N and S), which are the result of Leo’s migratory behaviour.


### Azimuth density per season
# store the information in a new data frame
DF <- data.frame(D=Leo_moving$azimuth,
                 V=Leo_moving$speed,
                 Season=Leo_moving$season)
# Define the direction as a circular object
DF$Dcirc <- as.circular(drop_units(DF$D), # this function does not "understand" units
                        rotation="clock",
                        units="radians",
                        type="angles",
                        modulo="asis",
                        zero=0,
                        template="geographic")

# remove missing values
DF <- DF[complete.cases(DF),]
# vector with the order of the factor levels, used to set the order of plotting
seasons <- levels(Leo_moving$season)
# store margin defaults before changing them
parDefaults <- par()
par(mar=rep(1,4))
# plot all the azimuths (circular histogram)
plot(DF$Dcirc, stack=T, shrink=1.6, pch=16, sep=0.05, col="grey")
# loop through seasons and plot a line density per season
for(i in 1:length(seasons)){
  # subset the azimuth
  x <- DF[DF$Season==seasons[i],'Dcirc']
  # calculate density and plot as a line
  lines(density(x, bw=180, kernel="vonmises"), lwd=2, lty=i)
  # draw an arrow showing mean and resultant length
  arrows.circular(mean(x), y=rho.circular(x), lwd=2, length=0.1, lty=i)
}
# add a legend
legend("bottomleft", lty=c(1,2,3,4), seasons, bty="n", cex=0.85)
# Note to plot (only contains directional information):
# - grey histogram: number of locations
# - lines: circular density function
# - arrows: mean direction; arrow length: mean resultant length, measure of variance of the data (wintering and breeding can be barely seen)
# - wintering might need better classification

# restore default margins
par(parDefaults)

### Wind rose, similar plot but showing azimuth AND speed, again per season
Leo_gg <- mutate(Leo_moving, speed_categorical = cut(set_units(speed, m/s), breaks = c(2, 5, 10, 15, 35)))
Leo_gg <- filter(Leo_gg, !is.na(season))

ggplot(Leo_gg) +
  coord_polar(start = pi) +
  geom_histogram(
    aes(x = set_units(azimuth, "degrees"),fill = speed_categorical),
    breaks = set_units(seq(-180L, 180L, by = 10L), "degrees"),
    position = position_stack(reverse = TRUE)) +
  scale_x_units(
    name = NULL,
    limits = set_units(c(-180L, 180L), "degrees"),
    breaks = (-2L:1L) * 90L,
    labels = c("S", "W", "N", "E")) +
  scale_y_continuous(name = NULL, limits = c(-1, 120), expand = c(0L, 0L)) +
  facet_wrap(~season) +
  scale_fill_ordinal("Speed [m/s]") +
  theme_linedraw()
# Note to plot:
# - migration has directionality with high speeds
# - breeding has only low speeds
# - wintering seems to contain some migration. But we just took "month: to categorize.


###______________________________________________________
### Turning angle distribution with speed and season ####

summary(Leo$turnangle)

### Speed ~ turning angle scatter plot, per season
ggplot(filter(Leo,!is.na(season)),aes(x = turnangle, y = speed)) +
  geom_point(color=alpha("black", 0.3))+
  scale_x_units("Turning angle", unit = "degrees", breaks = c(-2:2) * 90, expand = c(0L, 0L)) +
  scale_y_units(unit = "m/s")+
  facet_wrap(~season)+
  theme_linedraw()
# Note to plot:
# - during migration movement is directional and speeds are higher
# - in breeding speeds are low and movement is in all directions
# - here again we see that wintering contains some migration


### Wind rose, showing turning angles AND speed, again per season
pi_r <- set_units(pi, "rad")

ggplot(Leo_gg) +
  geom_histogram(
    aes(x = turnangle, fill = speed_categorical),
    position = position_stack(reverse = TRUE)) +
  scale_fill_ordinal("Speed") +
  coord_polar(start = pi) +
  scale_x_units(
    name = NULL,
    limits = c(-pi_r, pi_r),
    breaks = (-2L:1L) * (drop_units(pi_r/2)),
    labels = c("180", "360", "0", "90")) +
  scale_y_continuous(limits = c(-500L, 650L), breaks = c(0L, 250L, 500L)) +
  theme_linedraw()
# Note to plot:
# - probability of turning around is dependent on speed: the faster, the harder to turn
# - high speeds only have low turning angles
# - low speeds have all turning angles


#__________________________________________________
## CORRELATION IN DIRECTION & MOVEMENT PROCESS ####
#__________________________________________________
### Simulate 500 tracks with 1000 steps, 
### with different levels of correlation of azimuth (0.8 <= r <= 0.995)
set.seed(1512)
library(adehabitatLT)
rVals <- sort(rep((0.8 + log10(c(seq(1,100, length.out=10)))/10)[1:9], each=500))
rCRW <- lapply(lapply(rVals, simm.crw, date=1:1000, h=1), as, "Move") #convert to move
# calculate NSD for all tracks, from origin of trajectory
rNSD <- unlist(lapply(lapply(lapply(rCRW, coordinates), spDistsN1, pt=c(0,0)), "^", 2))
# calculate the mean NSD of all 500 trajectories of each value of r
mNSD <- tapply(rNSD, list(sort(rep(rVals,1000)), rep(1:1000, length(rVals))), mean)

### Plot NSD as a function of increasing number of steps for different correlations (r) in directions:
par(mar=c(5, 4, 4, 4) + 0.1)
plot(0,0, type="n", xlim=c(0,1300), ylim=c(0, max(mNSD)),
     bty="n", xlab="Step", ylab="Net square distance", xaxt="n")
axis(1, at=c(0,200,400,600,800,1000))
test <- apply(as.matrix(mNSD), 1, lines, x=1:1000)
text(cbind(rep(c(1250, 1100), length.out=length(row.names(mNSD))), mNSD[,ncol(mNSD)]),
     paste("r=", as.character(round(as.numeric(row.names(mNSD)),3)),sep=""), cex=0.5)
# Note to plot:
# - slope: how fast it gets away from the origin 

# group trajectories by correlation values
rCRW_ls <- split(rCRW, round(rVals, 3))
par(mfrow=c(3,3))
lapply(rCRW_ls[["0.8"]][sample(1:500, 3)], plot, pch=19, main="r = 0.8")
lapply(rCRW_ls[["0.975"]][sample(1:500, 3)], plot, pch=19, main="r = 0.975")
lapply(rCRW_ls[["0.995"]][sample(1:500, 3)], plot, pch=19, main="r = 0.995")
# Note to plot:
# - high r = very directional movement, steep increase in NSD per step made, gets away fast from origin
# - low r = movement in all direction around the origin, wiggly (brownian motion), takes long time to get away

# We saw from simulated trajectories how the correlation in direction directly affects 
# how long an animal takes to leave a certain area (its movement process/diffusive process from an area).
# To directly quantify this same parameter from an empirical trajectory we can use:
# - "reference" point -> NSD Net Squared Displacement
# - NO "reference" point -> FPT First Passage time

#___________________________________
# NET SQUARE DISPLACEMENT (NSD) ####
#___________________________________

# NSD = distance of each point in the trajectory from a given reference point (nest, colony..)
# NSD ~ time = how fast an animal moves away from a point of reference
# In this case the point of reference is the first point of the track, that is Leo's nest

# restore plotting device
dev.off()
par(parDefault)

# calculate and plot NSD
LeoNSD <- (st_distance(Leo,Leo[1,]))^2
LeoNSD <- units::set_units(LeoNSD, "km^2")
layout(matrix(c(1,1,2,3), ncol=2, byrow=T))
plot(mt_time(Leo), LeoNSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)", main="All data")

leoBreed08 <- Leo[which(Leo$season=="Breeding" & year(mt_time(Leo))==2008),]
leoBreed08NSD <- (st_distance(leoBreed08,leoBreed08[1,]))^2
units(leoBreed08NSD) <- units::make_units(km^2)

plot(mt_time(leoBreed08),leoBreed08NSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)", main="Breeding 2008")

leoWinter <- Leo[which(Leo$season=="Wintering"),]
leoWinter <- leoWinter[c(which(year(mt_time(leoWinter))==2008 &month(mt_time(leoWinter))%in%c(11,12)),
                         which(year(mt_time(leoWinter))==2009 &month(mt_time(leoWinter))%in%c(1,2,3))),]
leoWinterNSD <- (st_distance(leoWinter,leoWinter[1,]))^2
units(leoWinterNSD) <- units::make_units(km^2)
plot(mt_time(leoWinter),leoWinterNSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)",main="Winter 2008/2009")
layout(matrix(c(1), ncol=1, byrow=T))
# Note to plot:
# - "All data": shows nicely migration and winter/summer range, high nest fidelity and winter site fidelity
# - "Breeding": Leo moves very locally (it's squared, so < 20 m from nest)
# - "Winter08/09": shows that data is not classified correctly, data should be cut off at both ends


#__________________________________________________
# FIRST PASSAGE TIME (FPT) AND VARIANCE IN FPT ####
#__________________________________________________

# If you don't have a point of reference, every location of your trajectory can represent a "reference point"
# FPT calculates how long it takes to leave a circle with a radius r, centered on each location
# what is the mean FPT and the variation of FPT along the track?

library(adehabitatLT)
# functions in the adehabitatLT family assume planar coordinate systems where distances and angles are preserved
# therefore it cannot deal with spherical (geographic) coordinate system and we need to reproject
# the best is to use a local projection, centered in the middle of the track, to minimize the effects of spherical distortion on distances. 
# but this is a problem in large scale datasets, long distance migration in this case
# here we take the start and end point of its migration for a 2points equidistant projection
# you will learn a much better/automated way to do this on Friday!

Leo <- mt_read("Leo-65545.csv")

AEQ_2p <- "+proj=tpeqd +lon_1=-110.17999 +lat_1=55.64575 +lon_2=-78.10711 +lat_2=18.81255 +datum=WGS84"
Leoprj <- st_transform(Leo, crs=AEQ_2p)

# create "as.ltraj" class object
Leo_ltraj <- adehabitatLT::as.ltraj(xy=st_coordinates(Leoprj), 
                                    date=mt_time(Leoprj), 
                                    id=mt_track_id(Leoprj), typeII=T)
# calculate FPT, for different radii, in this case between 1 km and 1000 km in 150 exponential steps (samples more smaller radii than bigger)
fptLeo <- fpt(Leo_ltraj, radii=10^seq(3, 6, length.out=150), units="days")
# fpt object: n. columns = n.radii, n.rows = n.locations, values = fpt in days
str(fptLeo)
dim(fptLeo[[1]])

# calculate mean nb of days to leave each radii
meanFPT <- colMeans(fptLeo[[1]], na.rm=T)
radiiFPT <- attributes(fptLeo)$radii
plot(meanFPT ~ radiiFPT,
     type="l", lwd=2, xlab="Radii in meters",
     ylab="First passage time in days", log="xy")
# Note to plot:
# - with increasing radii size, on average it takes the animal longer to leave the circle
# - we are interested in the changes of slope


### variance of the log(FPT)
vars <- varlogfpt(fptLeo, graph=F)
plot(as.numeric(vars)~radiiFPT,
     type="l", lwd=1, lty=2,
     log="x", ylab="Variance of log first passage time",
     xlab="Radius in meters")
# Note to plot:
# - minima and maxima can indicate change in the movement process


### fitting LM to min/max peaks of variance of log(fpt)
plot(log10(meanFPT)~log10(radiiFPT),
     type="l", lwd=2, xlab="Log radii in meters",
     ylab="Log first passage time in days")
# fit a model to the largest valley, and largest peak of variance
lm1 <- lm(log10(meanFPT[1:which.min(vars[1:which.max(vars)])])~
            log10(radiiFPT[1:which.min(vars[1:which.max(vars)])]))
lm2 <- lm(log10(meanFPT[which.min(vars[1:which.max(vars)]):which.max(vars)])~
            log10(radiiFPT[which.min(vars[1:which.max(vars)]):which.max(vars)]))
abline(lm1, lty=2)
abline(lm2, lty=3)
text(4, 0.1, paste(signif(summary(lm1)$coefficients[2,1], 2),
                   "±",
                   signif(summary(lm1)$coefficients[2,2], 2)), pos=4, cex=0.75)
text(4, 1, paste(signif(summary(lm2)$coefficients[2,1], 2),
                 "±",
                 signif(summary(lm2)$coefficients[2,2], 2)), pos=4, cex=0.75)
# Note to plot:
# - flat slope (< 2): directional movement
# - steep slope (around 2): brownian movement


### breaks in the trend of the variance of log(fpt)
plot(as.numeric(vars)~radiiFPT,
     type="l", lwd=1, lty=2,
     ylab="Variance of log first passage time",
     xlab="Radius in meters", log="x")
breaks <- which(diff(floor(diff(as.numeric(vars))))==-1)+1
abline(v=radiiFPT[breaks])
# Note to plot:
# - besides the local minimum and the global maximum, there is a post maximum hump in the variance


### fitting LM to all segments identified by the breaks (changes in slope of variance of log(fpt))
plot(log10(meanFPT)~log10(radiiFPT),
     type="n", lwd=4, xlab="Log radii in meters",
     ylab="Log first passage time in days")

lm1 <- lm(log10(meanFPT[1:breaks[1]])~log10(radiiFPT[1:breaks[1]]))
lm2 <- lm(log10(meanFPT[breaks[1]:breaks[2]])~log10(radiiFPT[breaks[1]:breaks[2]]))
lm3 <- lm(log10(meanFPT[breaks[2]:breaks[3]])~log10(radiiFPT[breaks[2]:breaks[3]]))
lm4 <- lm(log10(meanFPT[breaks[3]:breaks[4]])~log10(radiiFPT[breaks[3]:breaks[4]]))
lm5 <- lm(log10(meanFPT[breaks[4]:length(as.numeric(vars))])~log10(radiiFPT[breaks[4]:length(as.numeric(vars))]))

abline(lm1, lty=2, lwd=1 + summary(lm1)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm2, lty=2, lwd=1 + summary(lm2)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm3, lty=2, lwd=1 + summary(lm3)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm4, lty=2, lwd=1 + summary(lm4)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm5, lty=2, lwd=1 + summary(lm5)$coefficient[2,1], col=alpha("black", 0.8))

lines(log10(meanFPT)~log10(radiiFPT),type="l", lwd=4, col=alpha("grey40", 0.8))
legend("bottomright",title="Radii (m)", lty=c(2,2,2,2,2),
       lwd=signif(c(1+summary(lm1)$coefficient[2,1],
                    1+summary(lm2)$coefficient[2,1],
                    1+summary(lm3)$coefficient[2,1],
                    1+summary(lm4)$coefficient[2,1],
                    1+summary(lm5)$coefficient[2,1]),2),
       c(paste(c(1000, round(radiiFPT[breaks],0))[1:2], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[2:3], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[3:4], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[4:5], collapse=" - "),
         paste(c(round(radiiFPT[breaks],0)[4], 100000), collapse=" - ")),
       bty="n", cex=0.75)
# Note to plot:
# - the different radii represent the different scales at which Leo is operating


### We can better understand these scales by looking at FPT in relation to time
par(mfrow=c(2,2))
for(i in c(1,2,3)){
  plot(fptLeo[[1]][,breaks[i]]~ mt_time(Leo), type="n",
       xlab="Time", ylab="FPT (days)",
       main=paste("Radius ", round(radiiFPT[breaks[i]],0), "meters"),
       bty="n")
  points(fptLeo[[1]][,breaks[i]]~ mt_time(Leo), pch=16, col=alpha("grey", 0.1))
  lines(fptLeo[[1]][,breaks[i]]~ mt_time(Leo))
}
par(mfrow=c(1,1))
# Note to plot:
# - 1st: ~80Km seems to be highlighting the migration times, when Leo migrates it moves in steps of about 80 km every couple of days before stopping.
# - 4th: ~3Km seems to be Leo's day range size, it never takes her more than 5-10 days to leave a 3Km circle.
# - 2nd, 3rd: ~20Km is probably about the home range. The longest amount of days could represent breeding (up to 100 days), the shortest spikes wintering, when Leo is a bit more mobile (up to 25 days).
# - The points at the bottom of plot 1,2, and 3 (shortest fpt) should represent migration.


#________________________________________________________________
# VARIANCE OF dBBMM (dynamic Brownian Bridge Movement Model) ####
#________________________________________________________________

# dBB allows to assigns a probability of where the animal could have been when we did not observe it (between two known locations).
# It assumes brownian motion between them and gives the probability of being anywhere between point 1 and 2 at any given time.
# Probability is based on the movement variance of the trajectory: the faster and more direct, the less the variance between the two points
# Variance calculated leaving out the 2nd of three locations and calculating its distance from the straight line connecting location 1 and 3

Leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
# we need to reproject the data to metres. Here we can use center=T: the center of the coordinate system is the center of the track (only possible with move objects)
Leroy <- spTransform(Leroy, center=T) 

# we choose a window size of odd n. of locs; within a window, 
# one global variance is estimated using all locations and local variances are estimated based on the margin.
# the window slides through the entire trajectory and uses the margin to find breakpoints and estimate variances. 
# Size of the window and margins define the expected scale of change in behaviour.
LeroyVar <- brownian.motion.variance.dyn(Leroy, location.error=25, window.size=71, margin=21)

# This method allows the variance to change along the track, that's why the brownian bridge is "dynamic" 
# so we can use boxplots to show the variance of the movement variance at different times of the day:
VarDat <- data.frame(var=getMotionVariance(LeroyVar), 
                     hour=hour(LeroyVar$study.local.timestamp)) # important to use local time for interpretation
boxplot(VarDat$var~VarDat$hour, xlab="Hour of the day", ylab="mean Brownian variance", pch="*")
# Note to plot:
# - Leroy moves around during the night, and sleeps during the day, it probably avoids people 
# - Low variance in the motion variance means the location I left out is easy to predict (either stationary, or moving with constant speed)

#__________________________________________________________________
# VARIANCE of dBGB (dynamic bi-Gaussian Bridge Movement Model) ####
#__________________________________________________________________

# The variance in dBGB decomposes the variance in two:
# forward-backward (parallel to movement) and side-ways component (perpendicular to movement)
# helps distinguishing stationary from movement (which resulted in similar dBB variance above)

LeroyBGB <- dynBGBvariance(Leroy, locErr=25, windowSize=31, margin=15)
VarDat <- data.frame(var=getMotionVariance(LeroyBGB), hour=hour(LeroyVar$study.local.timestamp))
VarDat$I_d <- ((VarDat$var.para-VarDat$var.orth)/(VarDat$var.para+VarDat$var.orth))
boxplot(I_d~hour, xlab="Hour of the day", ylab="dBGB variance index", data=VarDat, pch="*")
abline(h=0, lty=2)
# Note to plot:
# - zero: true brownian motion, diffusive in all directions as likely to go sideways than front or back.
# (almost nothing in nature is brownian, but error is, if data fall on 0 it could be error e.g. in the den)
# - positive values: directional movement, more front or back than side, less changes in direction and more in velocity
# - negative values: more likely to go sideways, there are more changes in direction and less in velocity


#__________________________
# SIMULATION OF TRACKS ####
#__________________________

### Advantages of simulating tracks ####
# underlying movement process is known and explicit. 
# great for developing and testing new analytical methods as they allow to test for performance in a rigorous way.

### Setting seeds ####
# Computers generally do not really generate random numbers, but they use a pseudo random number generator. 
# If we take the same starting point, the same series of numbers is generated
# which allows to make our simulations reproducible and allows debugging 
# this is done by setting a "seed"
# set.seed(number)
# remember this any time you are generating random numbers (e.g. sample(), rnorm() etc)

#______________________________________________________________
### 1. simulate uncorrelated random walk (Brownian motion) ####

set.seed(123) 

### simulate pure brownian walk
# uncorrelated and unbiased
# random distribution of steps with mean = 0, 0 or positive/negative steps in x or y equally possible
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmBrownDF <- data.frame(x=cumsum(rnorm(steps, 0, 1)),
                          y=cumsum(rnorm(steps, 0, 1)),
                          time=seq(start.time, start.time+duration, length.out=steps),
                          id="SimmBrown")
simmBrown <- mt_as_move2(simmBrownDF, coords = c("x", "y"),
                         time_column = "time", track_id_column = "id")

### simulate biased brownian walk 
# uncorrelated, biased
# mean different from 0, tendency towards a step of a certain value, here 0.1 
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmBiasDF <- data.frame(x=cumsum(rnorm(steps, 0.1, 1)),
                         y=cumsum(rnorm(steps, 0.1, 1)),
                         time=seq(start.time, start.time+duration, length.out=steps),
                         id="SimmBiased")
simmBias <- mt_as_move2(simmBiasDF, coords = c("x", "y"),
                        time_column = "time", track_id_column = "id")
simsB <- mt_stack(list(simmBrown, simmBias))

ggplot()+
  geom_sf(data = mt_track_lines(simsB), aes(color = id))+
  theme_linedraw()


#_______________________________________________
### 2. simulate uniform random distribution ####

set.seed(456) 

# Here steps are taken from an uniform distribution, any value is equally likely between min and max (not a normal distribution as above)
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmUdf <- data.frame(x=cumsum(runif(steps, min=-1.96, max=1.96)),
                      y=cumsum(runif(steps, min=-1.96, max=1.96)),
                      time=seq(start.time, start.time+duration, length.out=steps),
                      id="Simm")
simmU <- mt_as_move2(simmUdf, coords = c("x", "y"),
                     time_column = "time", track_id_column = "id")
ggplot()+
  geom_sf(data = mt_track_lines(simmU), aes(color = id))+
  theme_linedraw()
# Note to plot:
# - looks like a brownian random walk, but it isn't. 
# - unrealistic. This type of walk is never observed, as normally there is not the same probability to do all types of movement

### Comparison of distribution of speeds
hist(drop_units(mt_speed(simmU)), col="grey", xlim=c(0,1.1), main=NA, xlab="Speed",breaks="FD")
hist(drop_units(mt_speed(simmBrown)), breaks="FD", ylim=c(0,250), add=T, col=alpha("white", 0.5))
legend("topright",fill=c("white","grey"), legend=c("simm. Brownian", "simm. Uniform"), bty="n")
# Note to plot:
# - uniform: speeds are skewed to the right, suggesting that the animal is travelling mostly at high speeds => this is highly unnatural
# - brownian: more realistic speeds, normally speeds are left skewed

#____________________________________________________
### 3. simulate multivariate normal distribution ####

# This method allows to simulate tracks taking into account 
# a covariate structure between multiple variables
# For instance we said that there is an intrinsic relationship between speed and turning angle, so we can integrate it
# you can also have a look at the momentum package
set.seed(1512) 

# First, create multivariate normally distributed values with defined variance/covariance structure, taken from Leo
library(mvtnorm)
v <- set_units(Leo$speed, m/s)
turn <- data.frame(Vm=rowMeans(cbind(as.numeric(v)[-1], as.numeric(v)[-length(as.numeric(v))])), 
                   angle=drop_units(Leo$turnangle[-nrow(Leo)]))
turn <- turn[turn$Vm>2,]
turn <- turn[!turn$angle==0,]
sig <- cov(cbind(log(turn$Vm), log(abs(rad(turn$angle)))), use="complete.obs")
mV <- mean(log(turn$Vm), na.rm=T)
mT <- mean(log(abs(turn$angle)), na.rm=T)
rNum <- rmvnorm(10000, mean=c(mV, mT), sigma=sig)
dists <- exp(rNum[,1])
turnVec <- deg(exp(rNum[,2]))%%180
rNeg <- sample(1:length(turnVec), length(turnVec)*0.5)
turnVec[rNeg] <- turnVec[rNeg] * -1
plot(turnVec,dists, pch=16, col=alpha("grey", 0.5))
# Note to plot:
# - function provides realistic looking speeds and turning angles

## Second, use the simulated speeds and turning angles to simulate a trajectory
library(geosphere)
AZ <- (cumsum(c(176, turnVec))%%360)[-1]
start <- matrix(c(0,0), ncol=2)
track <- NULL
for(i in 1:length(AZ)){
  tmp <- destPoint(start, AZ[i], dists[i])
  track <- rbind(track, tmp)
  start <- tmp
}
plot(track, type="l", ylab="Latitude", xlab="Longitude", asp=1)
# Note to plot:
# This track is much less directional than Leo's original, because 
# although relation between speeds and turning angles is more realistic, the values lack autocorrelation
# in nature there is usually persistence in maintaining speed 


#__________________________________________
### 4. simulate correlated random walk ####

library(adehabitatLT)

set.seed(5323) 
# requires regular timesteps
crw <- simm.crw(1:100, r=.99)

plot(crw)
text(100,30,paste("Correlation coefficient of azimuth: ",
                  round(cor(crw[[1]]$abs.angle[-99:-100],
                            crw[[1]]$abs.angle[c(-1,-100)]), 2),
                  sep=""), pos=2, cex=0.75)
text(100,20,paste("Correlation coefficient of turning angle: ",
                  round(cor(crw[[1]]$rel.angle[c(-1, -99, -100)],
                            crw[[1]]$rel.angle[c(-1, -2 ,-100)]), 2),
                  sep=""), pos=2, cex=0.75)
# Note to plot:
# - direction is correlated (r=.99), there is a high consistency in direction
# - turning angles are uncorrelated

# adehabitatLT has different functions for simulated walks
# simm.brown simulates a brownian walk (same we did manually in point 1)
# also you can have a look at ?simm.levy ?simm.mou and ?simm.mba


#___________________________
### 5. randomized track ####

# this function creates random tracks that have the same start and end point as original track, 
# but with a shuffled sequence of the original segments
# therefore maintains step lengths and azimuths, but not turning angles
# make the random tracks realistic by maintaining some of the geometrical properties of the species' movement
# useful to create a null expectation based on an empirical trajectory

rand.seg <- function(x, repeats=1){
  coordinateDiffs <- cbind(diff(st_coordinates(x)[,1]-st_coordinates(x)[1,1]) , diff(st_coordinates(x)[,2]-st_coordinates(x)[1,2]))
  t <- replicate(repeats, apply(rbind(st_coordinates(x)[1,], coordinateDiffs[sample(1:(nrow(x)-1)),]), 2, cumsum), simplify=F)
  options(warn=-1)
  MO <- mt_stack(lapply(t, function(s){
    df <- data.frame(x=s[,1], y=s[,2], time=mt_time(x), as.data.frame(s),animal="RandTrack")
    mt_as_move2(df, coords=c("x","y"), time_column="time", track_id_column= "animal")
  }
  ),.track_combine="rename")
  options(warn=0)
  return(MO)
}


set.seed(43597)
rw_traj <- simm.crw(1:50, h=1, r=0.95, c(10,10))
spdf <- adehabitatLT::ltraj2spdf(rw_traj)
simm.df <- data.frame(spdf, track = paste0(attr(rw_traj[[1]], "id"), "_", attr(rw_traj[[1]], "id")))
rw <- mt_as_move2(simm.df, coords = c("x", "y"), time_column = "date", track_id_column = "track")
repl <- rand.seg(rw, 99)

ggplot()+
  geom_sf(data = mt_track_lines(repl),aes(color="RandomTracks"), alpha=0.5)+
  geom_sf(data = mt_track_lines(rw), aes(color="EmpiricalTrack"))+
  geom_point(aes(x=10, y=10, color="Start"))+
  geom_point(aes(x=tail(st_coordinates(rw),1)[1],y=tail(st_coordinates(rw),1)[2], color="End"))+
  scale_color_manual("",breaks = c("EmpiricalTrack", "RandomTracks", "Start","End"),
                     values = c(EmpiricalTrack = "blue", RandomTracks = "darkgrey", Start = "green", End="red"))+
  guides(color = guide_legend(order = 1, override.aes = list(shape = c(NA, NA,16,16), linetype = c("solid", "solid", "blank","blank"))))+
  theme_void()


# Final notes on simulated trajectories:
# Understanding movement as a continuous time process will eventually
# allow to predict and simulate trajectories more realistically.
# More complex ways of creating random trajectories are already provided 
# in the libraries ctmm and crawl.
# The capability to generate random trajectories can be used to define 
# null models for a wide range of inferential statistics, 
# but the characteristics of the random tracks 
# starkly influences how good or specific the inferences are.


#________________________________
# AUTO-CORRELATION STRUCTURE ####
#________________________________

# Trajectories, empirical or simulated, are often generated based on a correlation of their geometric properties.
# In almost all natural movements we can find some degree of persistence in maintaining the direction and often also speed. 
# (in fact correlated random walks, look more "natural" than the Brownian random walks 
# The correlation structure of a trajectory can not only inform about the underlying movement processes 
# but also be used to investigate changes in these processes, 
# which can be used to find breaks in the behaviour of the animal from the movement, defined as segmentation

# Autocorrelation is calculated as the correlation of values in a vector with a certain lag

### Correlation between any 2 subsequent speeds (autocorrelation of lag 1)
Leroy_mv2 <- mt_as_move2(Leroy)
speed_leroy <- drop_units(set_units(mt_speed(Leroy_mv2),m/s))
speed_leroy <- speed_leroy[!is.na(speed_leroy)]
cor(speed_leroy[-length(speed_leroy)], speed_leroy[-1])


### Correlation of all speeds separated by one segment (autocorrelation of lag 2)
lag <- 2
cor(speed_leroy[seq(1, length(speed_leroy)-lag, 1)],
    speed_leroy[seq(1+lag, length(speed_leroy), 1)])
# autocorrelation decreases with increasing lag, let's see for 100

### Autocorrelation of speeds between lag 1 to lag 100
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.test(speed_leroy[seq(1, length(speed_leroy)-lag, 1)],
                     speed_leroy[seq(1+lag, length(speed_leroy), 1)], method="spearman")$estimate)
  p <- c(p, cor.test(speed_leroy[seq(1, length(speed_leroy)-lag, 1)],
                     speed_leroy[seq(1+lag, length(speed_leroy), 1)], method="spearman")$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)
# Note to plot:
# - significant autocorrelation up to lag ~7-10, which in leroy corresponds to ~6h
# - oscillation correspond to periodicity in activity, Leroy is active for a while, and than inactive for a while

# different plot, same information
acf(speed_leroy)

### Autocorrelation of azimuth and turning angle
par(mfcol=c(1,2))
# coerce to circular object
Leroy_ll <- st_transform(Leroy_mv2, crs="EPSG:4326")
az0 <- as.circular(drop_units(mt_azimuth(Leroy_ll)), type="direction",
                   units="radians", template="geographic",
                   zero="0", rotation="clock", modulo="asis")
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.circular(az0[seq(1, length(az0)-lag, 1)],
                         az0[seq(1+lag, length(az0), 1)]))
  p <- c(p, cor.circular(az0[seq(1, length(az0)-lag, 1)],
                         az0[seq(1+lag, length(az0), 1)], test=T)$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient", ylim=c(-0.1, 1), main="Azimuth")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)


turn0 <- as.circular(drop_units(mt_turnangle(Leroy_ll)), type="angle",
                     units="radians", template="geographic",
                     zero="0", rotation="clock", modulo="asis")
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.circular(turn0[seq(1, length(turn0)-lag, 1)],
                         turn0[seq(1+lag, length(turn0), 1)]))
  p <- c(p, cor.circular(turn0[seq(1, length(turn0)-lag, 1)],
                         turn0[seq(1+lag, length(turn0), 1)], test=T)$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient",
     ylim=c(-0.1, 1), main="Turning angles")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)
# Note to plot:
# - there is little autocorrelation, Leroy moves randomly in all directions


### Correlation between lag in time and space (step length)
par(mfrow=c(1,1))
leroy <- filter(leroy,!st_is_empty(leroy))
plot(drop_units(mt_distance(leroy)[mt_time_lags(leroy)>set_units(17, "min")]),
     drop_units(mt_time_lags(leroy)[mt_time_lags(leroy)>set_units(17, "min")]),
     xlab="Distance in meters",
     ylab="Time lag in minutes", log="xy")
# Note to plot:
# - usually for large timelags we expect large distances, unless you sample at 1 Hz with 25 m gps sampling error
# - long distances (error) with very short time lags, which would result in very high speeds
# - this plot is a good sanity check


#________________________________________________
# IRREGULAR DATA: INTERPOLATION AND THINNING ####
#________________________________________________

# Most empirical trajectories have gaps/irregularity in the sampling, 
# which can be a problem with some analytical methods assuming regular sampling.
# To regularize movement data there are technically two solutions: 
# 1. reduce the data to the largest gap in time or 
# 2. interpolate positions in the gappy sections. 
# Both come with their problems, but especially the second one makes strong assumptions with the following consequences:

#__________________________________________
### Consequence of interpolating point ####

library(scales)
fishers <- mt_read(mt_example())
leroy <- filter_track_data(fishers, .track_id="M1")
# in move2 empty locations are keps with their timestamps, so interpolation is as easy as:
leroy <- mutate(leroy, locationType=ifelse(st_is_empty(leroy), "interpolated","true"))
leroyInt <- mt_interpolate(leroy)

ggplot()+
  geom_sf(data = mt_track_lines(leroyInt), color="grey50")+
  geom_sf(data = leroyInt,  aes(color=locationType))+
  scale_color_manual("",values = c(true = alpha("blue", 0.3), interpolated = alpha("firebrick", 0.5)))+
  theme_void()
# Note to plots:
# - gaps are biologically meaningful: Leroy is denning. By interpolating we are "forcing Leroy to move"

# What happens to the autocorrelation structure?
leroyInt <- filter(leroyInt, !st_is_empty(leroyInt)) # first and last locations are empty and cannot be interpolated
acf(drop_units(mt_speed(leroyInt)[-nrow(leroyInt)]), lag.max=100, main="speeds")
# Note to plots:
# - pattern in autocorrelation plot is result from the introduction of these linearly interpolated positions

# Final notes on interpolation:
# rather work with gaps and errors then introducing points (unless we are sure they were in a static location (e.g. den)
# or for smooth visualizations


#___________________________________________
### Consequence of subsampling/thinning ####

# Difference in speed due to differences in sampling
set.seed(7478)
# creating one track, and thin the same track to 1 every third observation
r.track <- as(simm.crw(1:1000, 1, 0.99), "Move")
r.trackThin <- r.track[seq(1,nrow(coordinates(r.track)),3),]
# compare the distribution of their speeds (non-parametric test to compare distributions)
ks.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))

hist(sqrt(speed(r.trackThin)), freq=F, xlim=c(0,2), breaks="FD", col=alpha("grey", 0.5), xlab="Square root transformed speed", main=NA)
hist(sqrt(speed(r.track)), freq=F, add=T, breaks="FD", col=alpha("white", 0.5))
legend("topleft", c("all positions", "every 3rd position"), fill=c("white", "grey"), bty="n")
# high segment values are only possible with high sampl freq, in the thinned track we underestimate distance covered in the same time lag
# ground speed should always come together with the rate at which it was sampled

# Their distribution is different but their means quite similar
# non-parametric test to compare ranks
wilcox.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))
# parametric test to compare means
t.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))

# Final notes on thinning, especially when comparing groups:
# Always important to compare the "operational sampling rate" how many fixes I get per unit time, 
# irrespective of the sampling schedule I choose
# (e.g. if we are comparing sexes, and one is always between leaves battery is lower and they get less fixes)


#____________________
# TRACK ANALYSIS ####
#____________________

# Compare movement metrics between individuals/groups
# Ways to incorporate autocorrelation structure

library(move2)
library(move)
library(units)
library(dplyr)
library(sf)
library(MASS)
library(scales)
library(adehabitatLT)
library(mgcv)

movebank_store_credentials("RBook", "Obstberg1")
bats <-  movebank_download_study(study_id = "Parti-colored bat Safi Switzerland",
                                 sensor_type_id = "radio-transmitter")

# Annotate speed, azimuth and turning angle to the trajectory
bats <- mutate(bats, speed = set_units(mt_speed(bats),m/s), 
               timelag=set_units(mt_time_lags(bats),day), 
               distance=set_units(mt_distance(bats),m))
# add column sex to event attribute
bats <- mt_as_event_attribute(bats, sex)
#check id column and summarize values per id
mt_track_id_column(bats)
#bats <- group_by_track_data(bats, individual_local_identifier) # this works in the github version, or 
bats <- group_by(bats, individual_local_identifier)
batsPerTrack <- summarize(bats,
                          sex=unique(sex),
                          medianSpeed=median(speed,na.rm=T),
                          timeTracked=sum(timelag,na.rm=T),
                          distanceTracked=sum(distance,na.rm=T)) #cumulative distance tracked
indData <- st_drop_geometry(batsPerTrack)
head(indData, 4)


### Test for differences in distance traveled per day between sexes
indData <- droplevels(indData)
boxplot(log10(I(distanceTracked/timeTracked)) ~ sex, data=indData, names=c("Males", "Females"),
        ylab=expression(paste(Log_10, " of cumulative distance in m per day", sep="")))
# Note to plot: there may be a small difference
indDataNoUnits <- drop_units(indData)
t.test(log(I((distanceTracked/1000))/timeTracked) ~ sex, data=indDataNoUnits)
# but this difference does not seem to be significant
# with a glm we can add co-variates
mod <- glm(sqrt(distanceTracked) ~ as.factor(sex) + timeTracked, data=indDataNoUnits)
# we check that assumptions are met
par(mfrow=c(2,2))
plot(mod, ask=F)
# and we see that the longer we track them, the higher the cumulative distance (of course)
summary(mod)


### Test for differences in speed between sexes
par(mfrow=c(1,1))
boxplot(log(medianSpeed) ~ sex, data=indDataNoUnits, names=c("Males", "Females"), ylab="Median speed in m/s")
# Note to plot: males seem to travel faster than females.. one obs with interesting negative speed, but we leave it for now
wilcox.test(medianSpeed ~ sex, data=indDataNoUnits)
# this time the difference seems to be significant

# Speeds and distances have naturally very long tails, and to run linear models we often have to transform our data
# Usually log or sqrt, but if we are not sure, boxcox is a nice function to find a right transformation of the data
# If we minimise the bias in the residuals the estimated slope is correct
bc <- boxcox(medianSpeed ~ as.factor(sex) + timeTracked, data=indDataNoUnits[-15,])
modII <- glm(I(medianSpeed^bc$x[which.max(bc$y)]) ~ as.factor(sex)+timeTracked, data=indDataNoUnits[c(-15),])
# plot(modII)
summary(modII)
# there is a difference between sexes in speed
# time tracking is not significant anymore, because is integrated in the speed calculation
# given what we said about the effect of sampling frequency we could check that too
bats$timelag_min <- set_units(bats$timelag, min) 
group_by(bats, individual_local_identifier, sex) %>% 
  summarize(medianTimelag=median(timelag_min, na.rm=T))


### Model using lme, glmm and gamm
# - can include complex correlation structures
# - can include random factors
# - can include different families (glm and gam)
# - can include smooth terms (gam)
# Interesting link: https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html
 
library(mgcv)
library(MASS)
library(nlme)
testDat <- data.frame(v=drop_units(bats$speed),
                      id=mt_track_id(bats),
                      time=mt_time(bats),
                      sex=bats$sex)
testDat <- testDat[which(testDat$v < 20 & testDat$v > 2),]
testDat <- arrange(testDat, id, time)
acf(testDat$v)
cs1 <- corAR1(0.2, form = ~ time|id) # correlation structure, value between -1 and 1
# we start with every subsequent observation being correlated with each other)

# using linear models
mod1 <- nlme::lme(log(v) ~ as.factor(sex), 
                  random=~1|id, #random factor
                  correlation=cs1, #corr structure
                  data=testDat)
summary(mod1)
acf(residuals(mod1))

# using generalized linear models
mod2 <- MASS::glmmPQL(log(v) ~ as.factor(sex), 
                      random=~1|id, #random factor
                      correlation=cs1, #corr structure
                      family=gaussian, #could choose diff family
                      data=testDat, verbose=FALSE)
summary(mod2)
acf(residuals(mod2))

# using generalized additive models
mod3 <- mgcv::gamm(log(v) ~ as.factor(sex), 
                   random=list(id=~1), #random factor
                   correlation=cs1, #corr structure
                   data=testDat, family=gaussian)
summary(mod3$gam)
summary(mod3$lme)
# significant but very low Rsq
acf(residuals(mod3$lme))
par(mfrow=c(2,2))
gam.check(mod3$gam)

# We tried to do the best we could do with only looking at the geometry of the track, 
# We used space and time, and comparisons between sex, all without looking at the environmental context yet, we'll do that later


