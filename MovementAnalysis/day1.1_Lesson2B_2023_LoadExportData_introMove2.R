########################################################
###                  AniMove 2023                    ###    
### Script by Kami Safi, Anne Scharf, Martina Scacco ###
########################################################
###  Intro to Move2, load-export data, manipulation  ###
########################################################

# install.packages("remotes")
# remotes::install_github("AniMoveCourse/animove_R_package")

library(move2)
library(dplyr)
library(sf)
library(readr)

# set wd to the data folder in your computer
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/Teaching/Animove/Animove2023_Canada/LecturesMaterial/MovementAnalysis/move2updated/data")

#_________________________________
# STRUCTURE OF A MOVE2 OBJECT ####
#_________________________________

### Event and track attribute tables
mv2 <- mt_read(mt_example())
class(mv2)
mv2

### Units
str(mv2)

### Variables definitions
movebank_get_vocabulary(mv2)

### Split and bind a list of individual move2 objects
mv2L <- split(mv2, mt_track_id(mv2))
mv2_2 <- mt_stack(mv2L)

#__________________________
# CHECKING PROPERTIES ####
#__________________________

?mt_is_track_id_cleaved

mt_is_track_id_cleaved(mv2) #do data from the same track occur consecutively
mt_is_time_ordered(mv2, non_zero = TRUE) #data are not ordered by default
mt_has_unique_location_time_records(mv2) #duplicates are allowed, checks if they exist
mt_has_no_empty_points(mv2) #empty locations are allowed, check if they exist
mt_is_time_ordered_non_empty_points(mv2, non_zero = FALSE) #each location with time and coords are ordered
mt_is_move2(mv2)


#________________________
# BASIC MANIPULATION ####
#________________________

### Retrieve names of recognised columns
mt_time_column(mv2)
mt_track_id_column(mv2)

### Retrieve information
mt_track_data(mv2)
mt_n_tracks(mv2)
mt_time(mv2) %>% head()
mt_track_id(mv2) %>% unique()
mt_n_tracks(mv2)

# geometry informations are retrieved using sf functions (they start with st_..)
sf::st_coordinates(mv2) %>% head()
sf::st_crs(mv2)
sf::st_transform(mv2, crs = "EPSG:32622") %>% sf::st_crs()

### Move attributes between track table and event table
?mt_as_event_attribute()

mv2

mv2_2 <- mt_as_event_attribute(mv2, c("individual-taxon-canonical-name","tag-local-identifier","study-name"))
names(mv2_2)
names(mt_track_data(mv2_2))

### set a different track identifier (it has to be a variable present in both event and track tables)
mv2_2 <- mt_set_track_id(mv2_2, value="tag-local-identifier")
mt_track_id_column(mv2_2)

### Manipulations on the TRACK TABLE
# 1. Manipulating rows and columns in the TRACK attribute table
?filter_track_data

filter_track_data(mv2, .track_id = c("F1","F2"))
filter_track_data(mv2, `tag-local-identifier` == "1072")

select_track_data(mv2, `individual-local-identifier`) # keep one column

mv2 <- mutate_track_data(mv2, sex = substr(`individual-local-identifier`,1,1)) # add column sex to track data
mt_track_data(mv2)

# 2. Summarise variables by TRACK attribute
mv2 %>%
  group_by_track_data(sex) %>%
  summarize(min_time=min(mt_time()), max_time=max(mt_time()))

mv2 %>%
  group_by_track_data(sex) %>% summarize(n = n())

### Manipulations on the EVENT TABLE
# 1. Manipulating rows and columns in the event attribute table
dplyr::filter(mv2, visible == TRUE)

dplyr::select(mv2, `event-id`, visible, timestamp)

# 2. Summarise variables by EVENT attribute
mv2 %>% 
  group_by(`behavioural-classification`) %>% 
  summarize(n = n())


### Data thinning/subsampling
mt_filter_per_interval(mv2, criterion = "first", unit = "2 hours")


#__________________________________
# CONVERT BETWEEN MOVE & MOVE2 ####
#__________________________________

### This is still a transition phase and some methods we will present still require the "old" move object
mv <- to_move(mv2)
class(mv)

mv2_2 <- mt_as_move2(mv)
class(mv2_2)

#_________________________________
# GETTING TRACKING DATA INTO R ####
#_________________________________

###________________________________________________
### 1. Directly downloading data from Movebank ####

# For the actual download we will look at three functions:
# movebank_download_study_info() movebank_download_deployment() movebank_download_study()

### store the movebank credentials for the most used account.
# you will be promped to set a keyring password
# by default this will be stored in the key list as service="movebank"
movebank_store_credentials("RBook", "Obstberg1")

### store credentials for another movebank account with a key name 
movebank_store_credentials("martina.scacco", key_name = "myOtherAccount")

## The "RBook" account is in this case used by default.
# If for a given session you want to use a different account you can specify it:
options("move2_movebank_key_name" = "myOtherAccount")

# To check which accounts are stored in keyring:
keyring::key_list()

### Remove accounts:
# the default one:
movebank_remove_credentials() #in this case martina.scacco became the default
keyring::key_list()
# a specific account with a key name:
movebank_remove_credentials(key_name = "movebank")
keyring::key_list()

# We restore it
movebank_store_credentials("RBook", "Obstberg1")

### Browse the Movebank database
# get the metadata of the studies visible from this account
(allStudies <- movebank_download_study_info())
rlang::global_entrace()
rlang::last_warnings()[[1]]$problems  #check the error, non-existent timestamp, problem of Movebank

# list studies for which we have download access
movebank_download_study_info(i_have_download_access=T)

# select only columns of interest and studies about bats
movebank_download_study_info() %>% 
  select(id, name, number_of_deployed_locations) %>% 
  filter(grepl("Parti-colored bat", name))

# retrieve all sensor ids recognised by Movebank
movebank_retrieve(
  entity_type = "tag_type",
  attributes = c("external_id", "id")
)

### Download information about a specific Movebank study
### get the id of the study (also part of the name but with unique match)
movebank_get_study_id("Parti-colored bat")

# get the metadata of the study
movebank_download_study_info(study_id = 1918503) %>% 
  print(width = Inf) #all cols

# check for reference data of animals, deployments and tags
# unique character string or study id, both work
movebank_download_deployment(study_id = 1918503)
movebank_download_deployment(study_id = "Parti-colored bat") #same result

# what sensors available in the study?
unique(movebank_download_deployment(study_id = 1918503)$sensor_type_ids)


##### a. Download LOCATION data as a move2 object ----
# download all data (all sensors, all animals)
movebank_download_study(study_id = 1918503) 

# download only data from specific sensor (name or id)
# it is recommendable to ALWAYS SPECIFY THE SENSOR to ensure that all attributes associated to the sensor get downloaded
movebank_download_study(study_id = 1918503, 
                        sensor_type_id = "radio-transmitter")

# download only a minimum amount of columns in event table (deployment, timestamp, geometry)
movebank_download_study(study_id = 1918503,
                        sensor_type_id = "radio-transmitter",
                        attributes = NULL) #this only acts on the event table, the track table is always downloaded in full

# specify one or more animals, either specifying the individual_id or the individual_local_identifier
movebank_download_study(study_id = 1918503,
                        attributes = NULL,
                        sensor_type_id = "radio-transmitter",
                        #individual_id = c(1918727, 1918739))
                        individual_local_identifier = c(239,360))

# and a specific time range
(sub <- movebank_download_study(study_id = 1918503,
                                attributes = NULL,
                                sensor_type_id = "radio-transmitter",
                                timestamp_start = as.POSIXct("2002-06-02 23:06:15", tz="UTC"),
                                timestamp_end = as.POSIXct("2002-06-11 22:18:25", tz="UTC")))
range(mt_time(sub))


##### b. Download NON-LOCATION data ----

# Note the empty geometry
(acc <- movebank_download_study(74496970,
                                sensor_type_id = "acceleration",
                                individual_local_identifier = "DER AR439"))
names(acc)

# Visualize and get basic stats of acceleration data (currently only for eObs tags)
# moveACC : https://gitlab.com/anneks/moveACC
# More with Hannah tomorrow.. stay tuned!


###________________________________________________________
### 2. Reading in a .csv file downloaded from Movebank ####

### We read in this dataset downloaded from Movebank directly as a move2 object
bats <- mt_read("Parti-colored bat Safi Switzerland.csv")
str(bats)
# Note: pay attention to the change in column name!
# When data are downloaded as csv from Movebank the column names change their separator
# compared to when they are directly downloaded from the API
# This unfortunately depends on Movebank so we cannot change it for now

###_________________________________________________
### 3. Creating a move object from any data set ####

### Sierit the stork, we read the data as a data.frame
df <- readr::read_csv("Sierit_DER AN858(eobs2561)_fromJan2020.csv",
                      col_types = list(timestamp = "c")) # to import the time column as character
str(df)

# first make sure the date/time is in POSIXct format
df$timestamp <- as.POSIXct(df$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
# also ensure that timestamps are ordered by individual (the object would be created anyway but it is good practice)
df <- df[order(df$`individual-local-identifier`, df$timestamp),]

# now create a move2 object from a data.frame
srk <- mt_as_move2(df, coords = c("location-long","location-lat"),
            crs = "EPSG:4326",
            time_column = "timestamp",
            track_id_column = "individual-local-identifier",
            na.fail = F) # allows or not empty coordinates

# alternatively, the re-ordering can also be done after creating the object
srk <- dplyr::arrange(srk, mt_track_id(srk), mt_time(srk))

srk

#________________________________________
# OUTLIERS AND DUPLICATED TIMESTAMPS ####
#________________________________________

### Same dataset, directly imported as move2 to have more columns automatically recognised
srk <- mt_read("Sierit_DER AN858(eobs2561)_fromJan2020.csv")

### OUTLIERS
# When downloading data using movebank_download_study(), the argument "remove_movebank_outliers" is TRUE by default
mt_movebank_visible(srk) %>% table()
# but if false, we could then remove them by:
srk <- mt_filter_movebank_visible(srk)

### EMPTY COORDINATES (in case of missed fixes or other sensors as ACC)
# Empty coordinates are allowed in move2, if we want to remove them we can:
sf::st_is_empty(srk) %>% table()
srk <- dplyr::filter(srk, !sf::st_is_empty(srk)) # Omit empty locations

### DUPLICATES
### By default, duplicated timestamps are allowed and kept in the dataset when downloading from Movebank
# But there are different ways to remove them

table(duplicated(mt_time(srk)))
# remove duplicated times without control of which locations are being removed (not recommended usually)
srk_noDups1 <- mt_filter_unique(srk, criterion = "first") #or sample
table(duplicated(mt_time(srk_noDups1)))
# often duplicates are caused by the same location transferred twice via different sources, once with more information than the other 
# add criterion="subsets" to remove duplicates that are subset of others (recommended)
srk_noDups2 <- mt_filter_unique(srk, criterion = "subsets") 
table(duplicated(mt_time(srk_noDups2)))
# BUT! some columns have different values (e.g. event-id), in fact:
# rows containing columns that are neither exact duplicates nor subsets don't get automatically excluded,
# How do we chack this?

# select the timestamps corresponding to duplicated events
tdup <- mt_time(srk)[duplicated(mt_time(srk))]
# and inspect the first duplicate
dups <- filter(srk, mt_time(srk)==tdup[1])
print(dups, width = Inf)

# some have different values, so we need to exclude them (note the use of mt_unique instead of mt_filter_unique):
srk_noDups <- srk[mt_unique(dplyr::select(srk, 
                            -c("event-id","eobs:battery-voltage","eobs:key-bin-checksum","import-marked-outlier"))), ]
table(duplicated(mt_time(srk_noDups)))

# hopefully soon added: select the duplicate with most information, i.e. with least columns containing NA


#______________________
# OUTPUTTING DATA ####
#______________________

### Save the move2 object as RData file
save(srk_noDups, file="SieritStork_cleaned.Rdata")
load("SieritStork_cleaned.Rdata") #cannot be assigned a new name, but can contain multiple objects

### Or as an rds file
saveRDS(srk_noDups, file="SieritStork_cleaned.rds")
sierit <- readRDS("SieritStork_cleaned.rds") # can only contain one object, but you can give it a name when loading


### Save as csv file
# at the moment with an extra step to deal with the geometry, but this will be solved soon
srk_noGeom <- st_drop_geometry(srk_noDups) %>% as.data.frame()
coords <- as.data.frame(st_coordinates(srk_noDups)) %>% as.data.frame() %>% 
  rename("location-long" = "X", "location-lat" = "Y")

srkDF <- cbind(srk_noGeom, coords)
head(srkDF)

write.table(srkDF, file="SieritStork_cleaned.csv", sep=",", row.names = FALSE)

mt_read("SieritStork_cleaned.csv")


## Save as a Shapefile
srk_sub <- dplyr::select(srk_noDups, timestamp, geometry, `individual-local-identifier`)
sf::st_write(srk_sub, getwd(), layer="stork", driver="ESRI Shapefile", delete_layer = TRUE) # for overwriting
sf::st_read("stork.shp")

### Save as a GeoPackage
sf::st_write(srk_sub, dsn="stork.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting
sf::st_read("stork.gpkg")

### Save as kml
sf::st_write(srk_noDups, dsn = "stork.kml", driver="kml", delete_layer = TRUE)
kml <- sf::st_read("stork.kml")

### NOTE: for more details see the vignettes https://cran.r-project.org/web/packages/move2/index.html
# and the articles in the website: https://bartk.gitlab.io/move2/

### To report issues with the package visit https://gitlab.com/bartk/move2/-/issues

### For users already familiar with move, we recommend the article 
# https://bartk.gitlab.io/move2/articles/convert.html 
# (which shows how to achieve a certain task in both packages)
