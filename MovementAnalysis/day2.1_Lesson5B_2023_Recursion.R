##############################################################
###                     AniMove 2023                       ###    
### Script by Chloe Bracis, Thomas Mueller, Martina Scacco ###
##############################################################
###                  Recursion Exercise                    ###
##############################################################
# chloe.bracis@gmail.com

library(recurse)
library(move2)
library(mapview)
library(sf)
library(ggplot2)
library(ggforce)
library(ggmap)
library(RgoogleMaps)
library(fields)
library(lubridate)
library(terra)
library(tidyterra)
library(scales)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

#_______________________________
### Functions for color palettes

# get continuous palette from blue (low values) to red (high values)
getContinuousPalette = function(n, alpha = 1)
{
  cols = alpha(brewer_pal(palette = "RdYlBu", direction = -1)(9), alpha)
  return( gradient_n_pal(cols)(seq(0, 1, length = n)) )
}

# discrete sequential palette, e.g. for one color per day 
getDiscreteSequentialPalette = function(n, alpha = 1)
{
  return(alpha(rainbow(n), alpha))
}

# store default graphic parameters
parDefault <- par()

#___________________________
### Initial exploration ####

### Read movebank elephant "Habiba" from Wall et al. 2014

elephants <- mt_read("Elliptical Time-Density Model (Wall et al. 2014) African Elephant Dataset (Source-Save the Elephants).csv")
habiba <- filter_track_data(elephants, .track_id = "Habiba")
mapview::mapView(mt_track_lines(habiba)$geometry) #as points

### Summarise data
summary(habiba)
unique(habiba$`study-timezone`)
range(mt_time(habiba))
diff(range(habiba$`utm-easting`))
diff(range(habiba$`utm-northing`))

### Color code data day by day
uniqueDays <- unique(day(mt_time(habiba)))

ggplot() + theme_void() +
  geom_sf(data = mt_segments(habiba), aes(color = as.character(day(mt_time(habiba))))) +
  scale_color_manual(values = getDiscreteSequentialPalette(length(uniqueDays))) +
  guides(color = "none")
# -> it appears that the looping behavior within one day

#______________________________________
### Examine the number of revisits ####

### project to azimuthal equidistant projection (important to project your data in metres before applying these methods)
# and we center the projection to the median of the data coordinates
apply(st_coordinates(habiba), 2, median)
AEQ_1p <- "+proj=aeqd +lat_0=0.5714825 +lon_0=37.5964800 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs"
habibaAEQD <- sf::st_transform(habiba, AEQ_1p)
st_crs(habibaAEQD)

### The recursion package does not deal with move2 objects yet, so we convert it to a move object (also possible with a data.frame)
habibaAEQD_mv <- to_move(habibaAEQD)

### Calculate recursions for a 500m Radius (you only need coords, time and radius)
r <- 500 # 500m 
habibavisits <- getRecursions(habibaAEQD_mv, radius = r, timeunits = "hours") #radius in metres

# examine output object
class(habibavisits)
names(habibavisits)
str(habibavisits)
head(habibavisits$revisits)

### Plot revisits in graphics
# plot their distribution
hist(habibavisits$revisits, breaks = 10, col = "darkgrey", main = "", xlab = "Revisits (radius 500m)")
# plot revisits on trajectory
plot(habibavisits, habibaAEQD_mv, alpha = 1, legendPos = c(4177000, 72200))
drawCircle(4180000, 72100, radius = r) # draw a 500 m circle as a reference for scale
text(4181500, 72100, "500m")

### Plot revisits on a background map
# note: these maps expects latlong coordinates, so we use habiba (not habibaAEQD)
bb <- sf::st_bbox(habiba)
names(bb) <- c("left","bottom","right","top")

m <- ggmap::get_map(location = bb, zoom = 13, source = "stamen", maptype = 'terrain', legend="topright")

habibaDF <- as(to_move(habiba), 'data.frame')
ggmap(m) + 
  geom_path(data = habibaDF, aes(x = coords.x1, y = coords.x2), color = "black", linewidth = 1) + 
  geom_point(data = habibaDF, aes(x = coords.x1, y = coords.x2), 
             color = getContinuousPalette(max(habibavisits$revisits), 0.5)[habibavisits$revisits])

### There are many more information in the recursion object:
# accessing the data frame for the recursion statistics
habibavisits$revisitStats[1:15,]

#______________________________________
### Examine the first passage time ####

# the first visits passes through the center of the circle, therefore
# first visit at each location = first passage time
# (this only works in the recurse package if you use single individuals, not available for multiple individuals)
habibavisits$firstPassageTime <- habibavisits$revisitStats$timeInside[habibavisits$revisitStats$visitIdx == 1]

hist(as.numeric(habibavisits$firstPassageTime), breaks = 20, col = "darkgrey", main = "FPT in 500 m", 
     xlab = "First passage time (hrs)")
# the histogram indicates a bimodal distribution potentially indicating two different behaviors
# natural split at 6 h

# Plot on map locations with first passage larger or smaller than 6 hrs
cutOff = 6
ggmap(m) + 
  geom_path(data = habibaDF, aes(x = coords.x1, y = coords.x2), color = "black", linewidth = 1) + 
  geom_point(data = habibaDF, aes(x = coords.x1, y = coords.x2), 
             color = alpha(ifelse(habibavisits$firstPassageTime > cutOff, "blue", "grey"),.5))


# It seems there is a different behavior at sections of the loop where they spend more time at night
boxplot(as.numeric(habibavisits$firstPassageTime) ~ hour(habiba$timestamp), 
		outline = FALSE, col = "grey", xlab = "Daytime (hrs)", ylab = "First passage (hrs)")


#__________________________________________
### Examine residence/utilization time ####

# Examine residence/utilization time, the sum of all visits around a focal point
head(habibavisits$residenceTime)
hist(as.numeric(habibavisits$residenceTime), breaks = 20, col = "darkgrey", main = "", xlab= "Residence time (hrs)")

# there seems to be a bimodal distribution, separated at about 20 hrs total visit time
colorsResTime <- as.character(cut(habibavisits$residenceTime, breaks = c(0, 10, 20, 30, 40), labels = c("grey", "gold", "darkorange", "firebrick3")))
ggmap(m) + 
  geom_path(data = habibaDF, aes(x = coords.x1, y = coords.x2), color = "black", linewidth = 1) + 
  geom_point(data = habibaDF, aes(x = coords.x1, y = coords.x2), 
             color = alpha(colorsResTime, 0.5))

# this looks different than the previous map
# in fact, even though the previous plot highlighted the long FPT of some locations at night
# the sum of the time spent in other locations is higher if these get rivisited many many times
# this map looks similar to the first map we plotted, which was showing the number of revisits

boxplot(as.numeric(habibavisits$residenceTime) ~ hour(habiba$timestamp), 
        outline = FALSE, col = "grey", xlab = "Daytime (hrs)", ylab = "Residence time (hrs)")
# in fact the relationship to time is opposite to that of fpt
# here we have longer cumulative residence time during the day

#_________________________________________
### Relate to environmental resources ####

### Look at NDVI at these areas to explore relation to resources
# note: the NDVI raster is in UTM coordinates
ndvi <- rast("MOD13Q1.A2014033.250m_16_days_NDVI.tif") / 10000

# in this case we reproject to the same crs of the ndvi
habibaUTM <- sf::st_transform(habiba, crs = sf::st_crs(ndvi)) 

# sanity check to see if raster and data overlap
ggplot() +
  geom_spatraster(data = ndvi) + 
  geom_sf(data = mt_track_lines(habibaUTM), linewidth = 1, color = "white")

# or more simply, but then we need to restore the default graphic parameters saved at the beginning of the session
# terra::plot(ndvi); points(to_move(habibaUTM)) 
# par(parDefault)


### Extract NDVI values at the location of the move2 object
habibaUTM$ndvi <- terra::extract(ndvi, st_coordinates(habibaUTM), method = "bilinear")[,1]

# Plot NDVI with data colored by revisits
coords <- st_coordinates(habibaUTM)
ggplot() + theme_bw() +
  geom_spatraster(data = ndvi) + 
  scale_fill_gradientn(colors = gray.colors(50), name = "NDVI") +
  geom_sf(data = mt_track_lines(habibaUTM), linewidth = 1, color = "white") +
  geom_sf(data = habibaUTM, color = alpha(colorsResTime, 0.5)) +
  geom_circle(data = data.frame(x=min(coords[,1]) + r, y=max(coords[,2]) - r), aes(x0=x, y0=y, r=r)) +
  annotate("text", label = "500 m", x=min(coords[,1])+r, y=max(coords[,2])+200) +
  coord_sf(datum = st_crs(ndvi)) #to change the graticule and show axes in correct projection
# Note: we can already see higher revisits in lighter pixels (higher ndvi)

### Explore NDVI for different residence times (> 20 h)
cutOff = 20
graphics::boxplot(habibaUTM$ndvi ~ habibavisits$residenceTime > cutOff, 
        outline = FALSE, col = "grey", notch = FALSE, 
        xlab = paste0("Utilization > ", cutOff," hrs"), ylab = "NDVI")

### Explore NDVI at different times of the day
boxplot(habibaUTM$ndvi ~ hour(habiba$timestamp), 
        outline = FALSE, col = "grey", xlab = "Time of day (hrs)", ylab = "NDVI")


#_____________________________
### Examine revisitations ####

### Time since last visit in days
hist(as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24), freq = TRUE, 
     xlab = "Time since last visit (days)", col = "darkgrey", main = "")

### Revisitation time after 1 week
returnsAfterOneWeek <- as.vector(na.omit(habibavisits$revisitStats$coordIdx[as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24) > 7]))

ggmap(m) + 
  geom_path(data = habibaDF, aes(x = coords.x1, y = coords.x2), color = "white", linewidth = 1) + 
  geom_point(data = habibaDF[returnsAfterOneWeek, ], aes(x = coords.x1, y = coords.x2), 
             color = alpha("blue", 0.2))
# Note: places revisted after a week

### Shorter revisitation times
hist(as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24), freq = TRUE, 
   xlab = "Time since last visit (days)", col = "darkgrey", main = "", 
   xlim = c(0, 4), ylim = c(0, 400), breaks = 70)

# there seems to be a hint for periodivity - possible suggests periodogram analyses


#_________________________________
### Play with different radii ####

### What if we pick a smaller radius?
r = 100 # 100m 

habibavisits <- getRecursions(habibaAEQD_mv, radius = r, timeunits = "hours") 

hist(habibavisits$revisits, breaks = 10, col = "darkgrey", main = "", xlab = "Revisits (100m radius)")

plot(habibavisits, habibaAEQD_mv, alpha = 1, legendPos = c(-8000, 4000))
drawCircle(-5000, 4000, radius = r)
text(-4000, 4000, "100 m", cex=0.6)
# check out the package vignette for an example of testing different radii


















