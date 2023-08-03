########################
# Territory size analysis
# Author: Ruth Simberloff (rsimberloff@gmail.com)
# Date: October 2021

# This script generates a utilization distribution for each bird using observed locations (spatial points), and then 
# finds the area of the UD at 75% kernel density.

#### setup ####

# attach packages
library(tidyverse)
library(adehabitatHR)

library(sp) # creates Spatial Points Data Frame, a special class of object needed for adehabitat
library(rgdal) # converts coordinates from long/lat to UTM
library(proj4) # A simple interface to lat/long projection and datum transformation of the PROJ.4 cartographic projections library

#### import data ####

raw_data <- read_csv("WCS2021_territories.csv")
glimpse(raw_data)


# CALCULATE TERRITORY SIZE #######################################################################

#### create sp objects ####

# Get coordinates

### transform coordinates to UTM
longlat_coords <- raw_data[3:4]
utm_coords <- project(as.matrix(longlat_coords), proj = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs")

### add UTM coordinates to data frame
raw_data <- bind_cols(raw_data, as_tibble(utm_coords), .name_repair = "universal") %>%
  rename(easting = V1, northing = V2)

# make a Spatial object

spdf <- SpatialPointsDataFrame(data = raw_data[1:2],
                               coords = utm_coords,
                               proj4string = CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs"))

#### estimate UDs ####

# estimation of kernel home range
# (first column is bird ID)
kuds <- kernelUD(spdf[,1], h = "href")


# estimation of home range areas at 75% density
areas <- kernel.area(kuds, percent = 75, unout = "m2")


### transpose areas
areas <- t(areas)
### convert area data from matrix to data frame (& keep the bird IDs as the first column)
areas <- as.data.frame(cbind(rownames(areas), as.data.frame(areas))) %>%
  rename(bird_id = "rownames(areas)")

# store home range contour as SpatialPolygonsDataFrame
vertices75 <- getverticeshr(kuds, percent = 75)

### export home range contours as ESRI shapefiles
### Used for visualization in Supplemental Figure 2
#writeOGR(vertices25, layer = "vertices25", dsn = "territory_verts_25_href", driver = "ESRI Shapefile")
#writeOGR(vertices50, layer = "vertices50", dsn = "territory_verts_50_href", driver = "ESRI Shapefile")
#writeOGR(vertices75, layer = "vertices75", dsn = "territory_verts_75_href", driver = "ESRI Shapefile")
#writeOGR(vertices95, layer = "vertices95", dsn = "territory_verts_95_href", driver = "ESRI Shapefile")

### export all coordinates as shapefile
#writeOGR(spdf, layer = "all_points", dsn = "all_points", driver = "ESRI Shapefile")


# Validating sample size ####################################################################

# this script draws random samples from the full dataset of observed locations for each bird. The sample size of the randomly drawn subset of locations
# increases by increments of 10 until the full sample (maximum 90 locations per bird) is reached. It draws 100 random samples at each sample size, and calculates the territory area
# from each randomly sampled subset of locations. If the number of observed locations is sufficient to estimate territory area accurately for that bird, the territory
# area is expected to plateau as the sample size increases


# extract the smoothing parameters from the above area estimations to use in sample size loop
## write a function to extract the value "h" from the "individual" slot in the SpatialObjectsDataFrame
get.kud.element <- function(individual) {
  individual@h$h
}

## map this function onto each individual in the object "kuds" and output as named list
library(purrr)
smooth_params <- map_dbl(kuds, get.kud.element)
print(smooth_params)

test_output <- tibble(bird_id = character(), ## make empty data frame to store loop results in
                      sample_size = integer(), 
                      area = numeric(),
                      h_value = numeric())

rep.counter <- 0 # a variable to keep track of the number of times it loops through the code

progress.bar <- txtProgressBar(min = 1, max = 100) # make a UI in the console to display loop progress


# WARNING: This loop may take up to 20 minutes to run

while(rep.counter <= 100) { ## The loop repeats until the rep counter reaches 100
  for(n in c(10, 20, 30, 40, 50, 60, 70, 80, 90)) { ## sample sizes to draw (10 locations, 20 locations, etc)
    for(id in names(smooth_params)) { ## 
      test_sample <- raw_data %>%
        filter(bird_id == id) %>% ## filters for each bird consecutively, using bird IDs from the object smooth_params
        dplyr::select(bird_id, easting, northing) %>%
        slice_sample(n = n, replace = FALSE) ## randomly sample n locations from each bird
      test_spdf <- SpatialPointsDataFrame(data = test_sample,
                                          coords = test_sample[2:3], 
                                          proj4string = CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs"))
      
      test_kuds <- kernelUD(test_spdf[,1], h = smooth_params[id]) ## generate kernel UD (using "href" method to select smoothing factor)
      ## pulls from smooth_params the h value calculated for each bird above (in the area estimation with the full data set)
      
      
      test_areas <- kernel.area(test_kuds, percent = 95, unout = "m2") ## calculate territory area from the kernel UD, selecting the 95% contour
      
      test_areas <- t(test_areas) ## the kernel area function spits out a wide df -- transpose to make long (and tidy)
      
      test_areas <- as.data.frame(cbind(rownames(test_areas), as.data.frame(test_areas))) %>%
        rename(bird_id = "rownames(test_areas)") ## Turn test_areas into normal data frame, keep row names as bird ID variable
      
      ## make output for each iteration
      iteration_output <- tibble(bird_id = test_areas$bird_id, sample_size = n, area = test_areas$"95", h_value = smooth_params[id])
      
      test_output <- bind_rows(iteration_output, test_output) ## 
      
      
      
    }}
  
  rep.counter <- rep.counter + 1 ## add 1 to the counter each time the code loops}
  setTxtProgressBar(progress.bar, rep.counter) ## set the progress bar to track how many times the code has looped
  
}


#### visualizing results of sample size calculations (Supplemental Figure 1)

se <- function(x){ sqrt(var(x)/length(x))}
summary_test_output <- test_output %>%
  group_by(bird_id, sample_size) %>%
  summarize(mean = mean(area), se = se(area)) 


p2 <- ggplot() +
  geom_point(data = filter(summary_test_output, bird_id == "WS.KO"), mapping = aes(x = sample_size, y = mean), size = 2) +
  geom_errorbar(data = filter(summary_test_output, bird_id == "WS.KO"), mapping = aes(x = sample_size, y = mean, ymin = (mean - se), ymax = (mean + se)), width = 0) +
  labs(x = "Number of locations", y = "Territory area      ") +
  theme_linedraw() +
  theme(text = element_text(family = "Roboto Condensed", size = 24))
p2
