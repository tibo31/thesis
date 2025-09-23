##########################
##### Load Libraries #####
##########################

library(sf)

################################
##### Load GeoSpatial Data #####
################################

##### COORDINATE REFERENCE SYSTEM (CRS) #####

#WORLD GEODETIC SYSTEM (1984) [PROJ4]
EPSG_4326 <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

##### COUNTRIEs #####
GeoDATA <- read_sf("data/Geo/COUNTRIEs/CCMN_GeoDATA.geojson", crs = EPSG_4326) #WORLD GEODETIC SYSTEM (1984)

##### Continental REGIONs #####
GeoRDATA_ContinentalREGIONs <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/ContinentalREGIONs/CCMN_GeoRDATA_ContinentalREGIONs.geojson", crs = EPSG_4326) 

##### Continental SIREGIONs #####
GeoRDATA_ContinentalSIREGIONs <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/ContinentalSIREGIONs/CCMN_GeoRDATA_ContinentalSIREGIONs.geojson", crs = EPSG_4326)
#NAME => FACTOR-LEVEL(s) => World Map
GeoRDATA_ContinentalSIREGIONs$NAME <- factor(
  x = GeoRDATA_ContinentalSIREGIONs$NAME,
  levels = c(
    "Eastern Africa", "Middle Africa", "Northern Africa", "Southern Africa", "Western Africa",
    "Central Asia", "Eastern Asia", "South-Eastern Asia", "Southern Asia", "Western Asia",
    "Eastern Europe", "Northern Europe", "Southern Europe", "Western Europe",
    "Caribbean", "Central America", "South America",
    "Northern America",
    "Australia and New Zealand", "Melanesia", "Micronesia", "Polynesia",
    "Other Countries/Territories"))

##### Geo. REGIONs #####
GeoRDATA_GeoREGIONs <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/GeoREGIONs/CCMN_GeoRDATA_GeoREGIONs.geojson", crs = EPSG_4326) 
#NAME => FACTOR-LEVEL(s) => World Map
GeoRDATA_GeoREGIONs$NAME <- factor(
  x = GeoRDATA_GeoREGIONs$NAME,
  levels = c("Australia and New Zealand", "Central and Southern Asia", "Eastern and South-Eastern Asia", 
             "Europe and Northern America", "Latin America and the Caribbean", "Northern Africa and Western Asia", 
             "Oceania*", "Sub-Saharan Africa", "Other Countries/Territories"))

##### More/Less Developed COUNTRIEs #####
GeoRDATA_MoreLess <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/DevLevels/MoreLess/CCMN_GeoRDATA_MoreLess.geojson", crs = EPSG_4326) 
#NAME => FACTOR-LEVEL(s) => World Map
GeoRDATA_MoreLess$NAME <- factor(
  x = GeoRDATA_MoreLess$NAME, 
  levels = c("More Developed Countries", "Less Developed Countries", "Other Countries/Territories"))

##### More/Less/Least Developed COUNTRIEs #####
GeoRDATA_MoreLessLeast <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/DevLevels/MoreLessLeast/CCMN_GeoRDATA_MoreLessLeast.geojson", crs = EPSG_4326) 
#NAME => FACTOR-LEVEL(s) => World Map
GeoRDATA_MoreLessLeast$NAME <- factor(
  x = GeoRDATA_MoreLessLeast$NAME, 
  levels = c(
    "More Developed Countries", "Less Developed Countries*", "Least Developed Countries", "Other Countries/Territories"))

##### Dev. Level(s) #####
GeoRDATA_LDC_LLDC_SIDS <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/DevLevels/LDC_LLDC_SIDS/CCMN_GeoRDATA_LDC_LLDC_SIDS.geojson", crs = EPSG_4326) 
#NAME => FACTOR-LEVEL(s) => World Map
GeoRDATA_LDC_LLDC_SIDS$NAME <- factor(
  x = GeoRDATA_LDC_LLDC_SIDS$NAME, 
  levels = c("Least Developed Countries* (LDC*)", 
             "Land Locked Developing Countries* (LLDC*)", 
             "Small Island Developing States* (SIDS*)", 
             "LDC | LLDC", "LDC | SIDS", "Other Countries/Territories"))

##### Income Level(s) #####
GeoRDATA_IncomeLevels <- read_sf( #WORLD GEODETIC SYSTEM (1984)
  "data/Geo/REGIONs/IncomeLevels/CCMN_GeoRDATA_IncomeLevels.geojson", crs = EPSG_4326)
#NAME => FACTOR-LEVEL(s) => World Map 
GeoRDATA_IncomeLevels$NAME <- factor(
  x = GeoRDATA_IncomeLevels$NAME, 
  levels = c("High-Income Countries", 
             "Upper-Middle-Income Countries", "Lower-Middle-Income Countries", 
             "Low Income Countries", "Other Countries/Territories"))