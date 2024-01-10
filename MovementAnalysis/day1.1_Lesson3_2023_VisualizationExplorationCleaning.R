########################################################
###                  AniMove 2023                    ###    
### Script by Kami Safi, Anne Scharf, Martina Scacco ###
########################################################
###       Visualization, exploration, cleaning       ###
########################################################

library(move2)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth) #needs rnaturalearthhires also installed
library(ggmap)
library(ggsn)
library(plotly)
library(mapview)
library(units)
library(lubridate)
library(classInt)
library(moveVis)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

#___________________________
# MAPPING MOVEMENT DATA ####
#___________________________

bats <- mt_read("Parti-colored bat Safi Switzerland.csv")

###_____________________
### Basic plotting ####

#### basic plots colored by individual (with graphics)
bats_ls <- split(bats, mt_track_id(bats))

cols <- rainbow(length(bats_ls))
plot(bats$geometry, xlab="Longitude", ylab="Latitude", type="n")
lapply(1:length(bats_ls), function(i){
  plot(bats_ls[[i]]$geometry, col=cols[i], pch=16, cex=0.5, add=T)})

#### basic plots colored by individual (with ggplot2)
ggplot() + theme_void() +
  geom_sf(data = bats) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`))

###________________________
### Plotting on a map ####

#### plot on the world boundaries
worldMap <- rnaturalearth::ne_countries(returnclass = "sf", scale = "large")

ggplot() + theme_void() +
  geom_sf(data = worldMap) +
  geom_sf(data = bats) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`))
  
# zoom in by subsetting the country
ggplot() + theme_void() +
  geom_sf(data = worldMap[worldMap$name == 'Switzerland',]) +
  geom_sf(data = bats) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`)) 

# or by cropping the map with a bounding box
bb <- st_bbox(bats)
exp <- 2

ggplot() + theme_void() +
  geom_sf(data = worldMap) +
  geom_sf(data = bats) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`)) +
  coord_sf(xlim = c(bb[1]-exp, bb[3]+exp), ylim = c(bb[2]-exp, bb[4]+exp), expand = F)


#### plot on different background maps
# request map data and then plot
names(bb) <- c("left","bottom","right","top") #the bounding box needs renaming
m <- ggmap::get_map(location = bb, 
                    zoom=9, source="stamen", maptype = "terrain")

ggmap(m) +
  coord_sf(crs = st_crs(bats)) + #set the coordinate system to that of the bats data
  geom_sf(data = bats, inherit.aes = FALSE) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`), inherit.aes = FALSE) +
  guides(color = "none")

# other possibilities
m_wc <- ggmap::get_map(location = bb, 
                       zoom=9, source="stamen", maptype = "watercolor")
m_bw <- ggmap::get_map(location = bb, 
                       zoom=9, source="stamen", maptype = "toner")

gg2 <- 
  ggmap(m_wc) +
  coord_sf(crs = st_crs(bats)) + #set the coordinate system to that of the bats data
  geom_sf(data = bats, inherit.aes = FALSE) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`), inherit.aes = FALSE) +
  guides(color = "none")
gg3 <- 
  ggmap(m_bw) +
  coord_sf(crs = st_crs(bats)) + #set the coordinate system to that of the bats data
  geom_sf(data = bats, inherit.aes = FALSE) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`), inherit.aes = FALSE) +
  guides(color = "none")

library(patchwork)
(gg2 + gg3) #+ plot_layout(guides = "collect") & theme(legend.position = "bottom") #to use if one of the plots has legend


# we can also add a scalebar, and choose a different map
xylim <- as.numeric(attributes(m)$bb)

myNiceMap <- 
  ggmap(m_wc) +
  coord_sf(crs = st_crs(bats)) + #set the coordinate system to that of the bats data
  geom_sf(data = bats, inherit.aes = FALSE) +
  geom_sf(data = mt_track_lines(bats), aes(color = `individual-local-identifier`), inherit.aes = FALSE) +
  guides(color = "none") +
  ggsn::scalebar(x.min = xylim[2], x.max = xylim[4], y.min = xylim[1], y.max = xylim[3],
           dist = 10, dist_unit="km", transform=T, model = 'WGS84', 
           st.size=2, anchor=c(x=8.65,y=46.93))

myNiceMap

#_______________________
# INTERACTIVE PLOTS ####
#_______________________

###____________
### Plotly ####

# by wrapping a ggplot in plotly
plotly::ggplotly(myNiceMap,
                 tooltip = c("individual-local-identifier"))


###____________
### Mapview ####

# by using mapview, after changing the class so that it is recognised as an SF object
batsSF <- bats
class(batsSF) <- class(bats) %>% setdiff("move2") # remove class "move2" from object
mapview::mapView(batsSF, zcol="individual-local-identifier", legend=F) #as points
mapview::mapView(mt_track_lines(bats), zcol="individual-local-identifier", legend=F) #as lines

# see this book "Making Maps with R" for more details on making maps in R: 
# https://bookdown.org/nicohahn/making_maps_with_r5/docs/introduction.html


#_________________________________
# ANIMATE TRACKS WITH moveVis ####
#_________________________________

# Not available on CRAN at the moment (maintenance needed) but from github repository: 
# devtools::install_github("16EAGLE/moveVis")

library(moveVis)
library(move)

# we will use the example dataset in the package, of storks in germany
data("move_data", package = "moveVis") # move class object
# if your tracks are present as data.frames, see df2move() for conversion
# if they are in move2 use to_move() to convert them

# align move_data to a uniform time scale
m <- align_move(move_data, res = 4, unit = "mins")
# reproject m to web mercator projection (default for base maps)
m <- sp::spTransform(m, "EPSG:3857") #spTransform is a sp based function that was used with move

# create spatial frames with a OpenStreetMap watercolour map
frames <- frames_spatial(m, path_colours = c("red", "green", "blue"),
                         map_service = "carto", map_type = "dark", alpha = 0.5) %>%
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>%
  add_scalebar() %>%
  add_timestamps(type = "label") %>%
  add_progress()

frames[[100]] # preview one of the frames, e.g. the 100th frame

# animate frames and save a gif in your home directory
animate_frames(frames, out_file = "moveVis.gif")

# Link to the Animove 2022 recording from Jakob Schwalb-Willmann:
# https://animove.org/elearning/ -> moveVis
# Link to the website:
# https://movevis.org/

# If you like animations, you could also check out the library "gganimate"


#________________________________________
# REMOVING OUTLIERS BASED ON MAPPING ####
#________________________________________

buff <- mt_read("Kruger African Buffalo, GPS tracking, South Africa.csv.gz")

# quick map to detect outlier
mapview::mapView(buff$geometry)

# get the position of the coordinate that has the max longitude
which.max(st_coordinates(buff)[,1])
# drop the point with the largest coordinate values
buff_noOut <- buff[-which.max(st_coordinates(buff)[,1]), ]

# The outlier has been removed and we save the cleaned object
mapview::mapView(buff_noOut$geometry)

save(buff_noOut, file="buffalo_noOut_less3.Rdata")


#_______________________________________________
# TEMPORAL ORGANIZATION OF THE TRAJECTORIES ####
#_______________________________________________

### Number of locations
nrow(bats)

### Number of tracks
mt_n_tracks(bats)

#__________________________________
### Time lag between locations ####

timeLags <- mt_time_lags(bats)
#by default it uses the most convenient unit but we can change it
head(timeLags)
timeLags <- units::set_units(timeLags, hours)
head(timeLags)

### Timelag, distance, speed, and direction are properties of the segment (2 locations, not one)
# therefore the number of valid values will be nrow(track) - 1 , and the last value will be NA
# the vector will contain as many NAs as the number of tracks, one at the end of each track
tail(timeLags)
table(is.na(timeLags))

### Distribution of time lags
summary(timeLags)
timeLags_h <- units::drop_units(timeLags) #hist is not compatible with time "units", we drop them
hist(timeLags_h, breaks=50, main=NA, xlab="Time lag in hours") 
arrows(24.5,587.5,20.7,189.7, length=0.1)
arrows(49.5,587.5,45.7,189.7, length=0.1)

### Distribution of timelags shorter than 1h
hist(timeLags_h[timeLags_h < 1], breaks="FD", main=NA, xlab="Time lag in hours")

#___________________________________
### Number of locations in time ####

ts <- mt_time(bats)
# transform timestamps into local time of the study for better interpretation
tsLocal <- lubridate::with_tz(ts, tzone="Europe/Zurich")
# N. of location per hour (geometry gets inherited)
bats %>% group_by(hour(tsLocal)) %>% 
  summarize(n())
# N. of locations per month and hour, rename columns and drop geometry
bats %>% group_by(Month = month(tsLocal), Hour = hour(tsLocal)) %>% 
  summarize(N.locations = n()) %>% 
  sf::st_drop_geometry()

#______________________________________________
# SPATIAL ORGANIZATION OF THE TRAJECTORIES ####
#______________________________________________

#__________________________________
### Distance between locations ####

dist <- set_units(mt_distance(bats), m)
summary(dist)
hist(dist)

#_______________________________
### Speed between locations ####

speeds <- set_units(mt_speed(bats), m/s)
summary(speeds)
hist(drop_units(speeds), breaks="FD")

### Have a look at the realistic speeds (e.g. < 20m/s)
speedsRealistic <- speeds[speeds < set_units(20, m/s)] # remember to set the unit before the operation
hist(drop_units(speedsRealistic), xlab = "Speed [m/s]", breaks = "FD") # this is the common shape for speeds

### Speed vs timelags
speedVsTimeLag <- data.frame(timeLag = timeLags, speeds = speeds)
speedVsTimeLag <- speedVsTimeLag[speedVsTimeLag$timeLag < set_units(10, hour) & speedVsTimeLag$speeds < set_units(20, m/s),]
# with longer timelags the speeds seem lower... you will learn a lot more about this in the following days.
plot(speedVsTimeLag$timeLag, speedVsTimeLag$speeds, xlab='Time lag', ylab='Speed', pch=19) # when "units" are recognised, they are automatically added to the plot

### Plot highlighting segments with high and low speeds

bats <- mt_read("Parti-colored bat Safi Switzerland.csv")
# select bat 191
unique(mt_track_id(bats))
bat191 <- filter_track_data(bats, .track_id = "191")
# store speed
v <- set_units(mt_speed(bat191), m/s)
# find 9 colors for 10 breaks in the speed vector, the myPal vector has the same length as the number of segments
myBreaks <- classIntervals(v[!is.na(v)], n=10, style="equal")$brks
#myPal <- as.character(cut(v, breaks = myBreaks, labels = grey.colors(10)))
ggplot() +
  #geom_sf(data = mt_segments(bat191), linewidth = 1.5, color = myPal) +
  #two alternative ways to assign colors, either with myPal and the line above, or with the line below and scale_color_gradient2
  geom_sf(data = mt_segments(bat191), aes(color = drop_units(v)), linewidth = 1.5) +
  xlab("Longitude") + ylab("Latitude") +
  scale_color_gradient2(low = "black", mid = "#434343", high = "white", 
                        midpoint = drop_units(median(v, na.rm=T)),
                        breaks = drop_units(myBreaks),
                        name = "Speed [m/s]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # get rid of panel grids
        panel.background = element_rect(fill = '#434343')) # Change panel background

#_____________________________________________________________
### Direction of movement / azimuth / heading of movement ####

# NOTE: the words heading or bearing are mostly referred to the direction of body axis
# When analysing tracking data we are not observing the animal but only its movement, so we can only now the direction of movement (not body orientation).
# BUT: many devices record "heading" which is the orientation of the tag (not the direction of movement). SO if the tag does not shift its position it is an indication of body orientation.
direction <- mt_azimuth(bats) # Angles in radians relative to the N
head(direction)
summary(direction)
hist(direction,  breaks = 18, xlab="Direction of movement", main = NA)
# they seem to go in all directions

#______________________
### Turning angles ####

### Turning angles are a property of 2 segments (3 locations)
turnAngles <- mt_turnangle(bats) 
# angles in radians relative to the previous step, with 2 NAs per track (first and last value of each track)
turnAngles[bats$`individual-local-identifier` == 191]
summary(turnAngles)

hist(turnAngles, breaks = 18, xlab="Turning Angle", main = NA)
# this shape can indicate movement along a linear structure

turnAnglesBuf <- mt_turnangle(buff)
hist(turnAnglesBuf)
# this is the common shape for turning angles


#__________________
# MISSED FIXES ####
#__________________

# load data from Leroy the fisher
fishers <- mt_read(mt_example())
leroy <- filter_track_data(fishers, .track_id = "M1")
#leroy <- mt_read(list.files(system.file("extdata", package = "move2"), full.names = T))
# pick one year
leroy <- leroy[year(mt_time(leroy))==2009,]
# get the time stamps of entries without location
missedFixes <- data.frame(time = mt_time(leroy[sf::st_is_empty(leroy),]), status = "Not Successful")
# add the time stamps of positions obtained
obtainedFixes <- data.frame(time = mt_time(leroy[! sf::st_is_empty(leroy),]), status = "Successful")

# rbind them and change time to local time zone of the study
df <- rbind(missedFixes, obtainedFixes)
df$time <- lubridate::with_tz(df$time, tz="America/New_York")

# Plot histogram that is filled out with proportions and is binned per hours
ggplot(df, aes(x=hour(time), fill=status)) +
  geom_histogram(binwidth=1, position='fill') +
  xlab("Hour of the day") + ylab("Proportion") +
  scale_fill_grey()

## histogram filled out with proportions binned per day
ggplot(df, aes(x=time, fill=status)) +
  geom_histogram(binwidth=24*60*60, position='fill') +
  xlab("Hour of the day") + ylab("Proportion") +
  scale_fill_grey()



                                     