## Process landscape level variable for Voronoi polygons

# update and load needed packages
source("./Scripts/package_load.R")

# Variables to calculate --------------------------------------------------

# population density - census
# proportion ag - land cover (Ches Bay)
# proportion habitat - land cover (Ches Bay)
# impervious surface - land cover (Ches Bay)
# road density - TIGER/Lines

# average temperature - can we have tie varying?
# NDVI - will we have seasonal?


# Load spatial data -------------------------------------------------------

# polygons
hex_polies <- st_read(dsn = "./Data/visit_counts_201801_202001",
                      layer = "visit_counts_201801_202001") |>
  filter(POI_count > 0)

# extract block level 2020 census data for each state(s)
# deer were collared in 2018 and moved for 1-2 years so I think 2020 is ok
# state_block_pop <- tidycensus::get_decennial(
#   geography = "block",
#   variables = "P1_001N",
#   state = "MD",
#   geometry = TRUE,
#   year = 2020,
#   show_call = TRUE
# )

# st_write(state_block_pop, "./Data/GIS/MD_2020_pop", 
#          layer="2024-04-10_MD_2020_block_pop",
#          driver = "ESRI Shapefile")

# Read in saved file from commented out snippet above
state_block_pop <- st_read("Data/GIS/MD_2020_pop/", 
                           layer = "2024-04-10_MD_2020_block_pop")

# load landcover
# data for howard county available here:
# https://www.chesapeakeconservancy.org/projects/cbp-land-use-land-cover-data-project
lc <- rast("Data/GIS/howa_24027_lulc_2018/howa_24027_lu_2018.tif")

# load roads shapefile
# shapefile available here: https://www2.census.gov/geo/tiger/TIGER2023/ROADS/
# make valid helps with some inconsistencies with geometries
howard_roads <- st_read("./Data/GIS/tl_2023_24027_roads/", 
                        layer = "tl_2023_24027_roads") %>%
  st_transform(crs = st_crs(state_block_pop)) %>%
  st_make_valid(.)

# load trails
# data available here:https://data.howardcountymd.gov/InteractiveMap.html
howard_trails <- st_read("./Data/GIS/Trails.kml")

# change crs (wouldn't let me pipe this with kml)
howard_trails <- st_transform(howard_trails, crs = st_crs(state_block_pop))


# Building area -------------------------------------------------------

# building footprint area was already calculated by Szandra

# Road density ------------------------------------------------------------

## In parallel

polies_rp <- st_transform(hex_polies, crs = st_crs(howard_roads))

# set cluster
cl <- parallel::makeCluster(detectCores() - 2)
# register cluster for foreach loop
doParallel::registerDoParallel(cl = cl)

# run the loop in parallel
roads <-foreach(i = 1:nrow(polies_rp), .combine=rbind,
                .export = c("howard_roads", "polies_rp"),
                .packages = "sf") %dopar% {
                  
                  road_length <- st_intersection(howard_roads, polies_rp[i,]) %>%
                    st_length(.) %>%
                    sum(.) %>%
                    as.numeric(.)
                  
                  road_length_km <- road_length/1000
                  
                  if(road_length == 0){
                    road_density <- 0
                  } else {
                    road_density <- road_length_km/
                      (as.numeric(st_area(polies_rp[i,])) / 1e6)
                  }
                  
                  polies_roads <- c(road_length, road_density)
                  
                  return(polies_roads)
                  
                }

stopCluster(cl)

colnames(roads) <- c("road_length", "road_dens")
  
  
# Trail Density -----------------------------------------------------------

polies_rp <- st_transform(hex_polies, crs = st_crs(howard_trails))

# set cluster
cl <- parallel::makeCluster(detectCores() - 2)
# register cluster for foreach loop
doParallel::registerDoParallel(cl = cl)

# run the loop in parallel
trails <-foreach(i = 1:nrow(polies_rp), .combine=rbind,
                 .export = c("howard_trails", "polies_rp"),
                 .packages = "sf") %dopar% {
                   
                   trail_length <- st_intersection(howard_trails, 
                                                   polies_rp[i,]) |>
                     st_length() |>
                     sum() |>
                     as.numeric()
                   
                   trail_length_km <- trail_length/1000
                   
                   if(trail_length == 0){
                     trail_density <- 0
                   } else {
                     trail_density <- trail_length_km/
                       (as.numeric(st_area(polies_rp[i,])) / 1e6)
                   }
                   
                   polies_trails <- c(trail_length, trail_density)
                   
                   return(polies_trails)
                   
                 }

stopCluster(cl)

colnames(trails) <- c("trail_length", "trail_dens")
  
# Land cover --------------------------------------------------------------

polies_albers <- st_transform(hex_polies, crs = st_crs(lc))

# Extract info from each cell from within each v_poly
lc_poly <- terra::extract(lc, polies_albers)

# Calculate the proportion of the total cells
# The margins option says to do this for each row and not the whole table
lc_proportions <- as.data.frame.matrix(prop.table(table(lc_poly), margin = 1))

# see Appendix C: 
# https://cicwebresources.blob.core.windows.net/docs/LU_Classification_Methods_2017_2018.pdf

# Change column names just to clean it up and work with it better
colnames(lc_proportions) <- c("Water", "Impervious_Roads", "Impervious_Structures",
                              "Impervious_Other", "Tree_over_Impervious",
                              "Turf_Grass", "Pervious_Developed", "Tree_over_Turf",
                              "Forest", "Tree_Other", "Natural_Succession",
                              "Cropland", "Pasture_Hay", "Extractive", 
                              "Wetland_Riverine", "Wetland_Terrene")

# Make some combined land cover types

# Habitat - forest, tree other, and natural succession
lc_proportions$habitat <- lc_proportions[,"Forest"] + lc_proportions[,"Tree_Other"] + 
  lc_proportions[,"Natural_Succession"] 

# Wetlands - two wetlands combined
lc_proportions$wetlands <- lc_proportions[,"Wetland_Riverine"] + 
  lc_proportions[,"Wetland_Terrene"]

# Impervious - roads, structures, other, tree over impervious
lc_proportions$impervious <- lc_proportions[,"Impervious_Roads"] + 
  lc_proportions[,"Impervious_Structures"] + 
  lc_proportions[,"Impervious_Other"] + 
  lc_proportions[,"Tree_over_Impervious"]

# Agriculture - crops and hay
lc_proportions$agriculture <- lc_proportions[,"Cropland"] + 
  lc_proportions[,"Pasture_Hay"]

  
# Landscape Metric --------------------------------------------------------

poly_rp <- st_transform(hex_polies, st_crs(lc))

# reclassify habitat into a single habitat raster
# 41; 65; 75; 95 Forest (FORE)
# 42; 64; 74; 94 Tree Canopy, Other (TCOT)
# 16; 54-56 Natural Succession (NATS) 

# change the raster to TRUE or FALSE if habitat
habitat <- lc %in% 
  c(41,65,75,95,42,64,74,94,16,54:56)

# empty list to store results
contig <- list()

for(i in 1:nrow(poly_rp)){
  # crop to the individual polygon
  habitat_crop <- mask(crop(habitat, poly_rp[i,]), poly_rp[i,])
  
  # calculate split index for polygon
  # calculates a category for FALSE so we remove that one
  contig[[i]] <- lsm_c_contig_mn(habitat_crop, directions = 4) |>
    filter(class == 1) |>
    pull(value)
}

# gives empty lists for where there was no patches
# make empty lists equal 0
contig_fix <- lapply(contig, function(x){
  if(length(x)==0){
    x <- 0
  } else {
    x
  }
})

# crunch the list down to a vector
contig_final <- data.frame(contig = do.call("rbind", contig_fix))


# Surface Temperature -----------------------------------------------------

# refer to GEE script to acquire Landsat panels
st <- rast("./Data/GIS/Day_monthlymean/April_2018_Day_Mean_LST.tif")

# just for transforming layer outside loop
poly_rp <- st_transform(hex_polies, crs = st_crs(st))

# get list of raster files
st_files <- list.files("./Data/GIS/Day_monthlymean/")

# split the filename apart to get the month and year
file_split <- strsplit(st_files, "_")

# keep only the month and year and then turn it into a date object
# we will use this to order out loop in the same order as the poly sf columns
convert_date <- lapply(file_split, function(x){
  as.Date(zoo::as.yearmon(paste(x[1:2], collapse = " "), format = "%B %Y"))})

# now we have a vector of dates that are the dates of the raster files
dates <- do.call("c", convert_date)

# reorder files to match our sampling order
files2 <-st_files[order(dates)]

# loop through and extract the mean temp per poly
# putting `order` in loop will call the raster file in chronological order

st_mean <- matrix(NA, ncol = length(dates), nrow = nrow(hex_polies))

for(i in 1:length(files2)){
  
  r <- rast(paste0("./Data/GIS/Day_monthlymean/",files2[i]))
  
  st_mean[,i] <- terra::extract(r, poly_rp, fun = mean)[,2]
  
}

colnames(st_mean) <- paste0("st_mo_", 1:ncol(st_mean))

# Polygon Area ------------------------------------------------------------
poly_area <- st_area(hex_polies) |>
  drop_units()

poly_area_km <- data.frame(poly_area = poly_area / 1e6)


# Combine results ---------------------------------------------------------

# remove columns that are already named this way
# basically replace the old columns with this new data

hex_polies_reduced <- hex_polies |>
  dplyr::select(GRID_ID_p:pop_total)

# put it all together
results <- dplyr:::bind_cols(hex_polies_reduced,
                             lc_proportions,
                             contig_final,
                             roads,
                             trails,
                             st_mean,
                             poly_area_km)


saveRDS(results, "./Data/2025-03-24_full_data.rds")

