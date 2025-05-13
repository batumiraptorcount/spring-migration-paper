library(tidyverse)
library(terra)
library(leastcostpath)
library(parallel)
library(ggdark)
library(ggfx)
library(sf)

dem <- rast("~/Dropbox/BRC/Research/Spring Paper/Map/Georgia_DEM_1200x800.tif")
dem[dem < 1] <- NA

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


gradient_autumn <- seq(from = max(as.matrix(dem), na.rm = TRUE), to = min(as.matrix(dem), na.rm = TRUE), length.out = 800)
m_autumn <- rep.col(gradient_autumn, 1200)
m_autumn <- as.matrix(imrotate(as.cimg(m_autumn), -15, interpolation = 2, boundary = 1))[1:800, 1:1200]
dem_autumn_pressure <- rast(m_autumn, extent = ext(dem))

gradient_spring <- rev(gradient_autumn)
m_spring <- rep.col(gradient_spring, 1200)
m_spring <- as.matrix(imrotate(as.cimg(m_spring), -15, interpolation = 2, boundary = 1))[1:800, 1:1200]
dem_spring_pressure <- rast(m_spring, extent = ext(dem))

spring_min_y <- 4445581
spring_x <- seq(from = -10947.54, to = 497155, len = 100)
spring_locs <- expand.grid(x = spring_x, y = spring_min_y)

autumn_max_y <- 4866765
autumn_x <- seq(from = 48693.24, to = 497155, len = 100)
autumn_locs <- expand.grid(x = autumn_x, y = autumn_max_y)

coordinates_spring <- spring_locs[, c("x", "y")]
points_spring <- lapply(1:nrow(coordinates_spring), function(i) {
  sf::st_point(as.numeric(coordinates_spring[i, ]))
})
locs_spring <- sf::st_sf(geometry = sf::st_sfc(points_spring, crs = terra::crs(dem_utm)))

coordinates_autumn <- autumn_locs[, c("x", "y")]
points_autumn <- lapply(1:nrow(coordinates_autumn), function(i) {
  sf::st_point(as.numeric(coordinates_autumn[i, ]))
})
locs_autumn <- sf::st_sf(geometry = sf::st_sfc(points_autumn, crs = terra::crs(dem_utm)))


calculate_lcps <- function(season) {
  if (season == "autumn") {
    dem_used <- dem + dem_autumn_pressure
  } else if (season == "spring") {
    dem_used <- dem + dem_spring_pressure
  } else {
    dem_used <- dem
  }
  plot(dem_used)
  dem_utm <- project(dem_used, "EPSG:32638")

  slope_cs <- create_slope_cs(x = dem_utm, cost_function = "tobler", neighbours = 32, exaggeration = TRUE)
  slope_cs <- add_global_stochasticity(slope_cs, percent_quantile = 0.2)

  num_cores <- detectCores() - 1

  if (season == "autumn") {
    create_lcp_for_index <- function(i) {
      create_lcp(x = slope_cs, origin = locs_autumn[i, ], destination = locs_spring, cost_distance = TRUE)
    }
    lcps_list <- mclapply(1:nrow(locs_autumn), create_lcp_for_index, mc.cores = num_cores)
  } else if (season == "spring") {
    create_lcp_for_index <- function(i) {
      create_lcp(x = slope_cs, origin = locs_spring[i, ], destination = locs_autumn, cost_distance = TRUE)
    }
    lcps_list <- mclapply(1:nrow(locs_spring), create_lcp_for_index, mc.cores = num_cores)
  }
  lcps <- do.call(rbind, lcps_list)
  return(lcps)
}

lcps_autumn <- calculate_lcps("autumn")
lcps_spring <- calculate_lcps("spring")

plot(lcps)

st_write(lcps_spring, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_rotated_15.shp", append = FALSE)
st_write(lcps_autumn, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_rotated_15.shp", append = FALSE)

lcps_autumn_short <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

lcps_autumn_short_rotated <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_rotated.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

lcps_autumn_short_rotated <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_rotated_15.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

lcps_spring_short <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

lcps_spring_short_rotated <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_rotated.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

lcps_spring_short_rotated <- st_read("~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_rotated_15.shp") %>%
  group_by(fromCll) %>%
  slice_min(n = 5, order_by = cost) %>%
  st_jitter()

st_write(lcps_spring_short, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_shortest.shp", append = FALSE)
st_write(lcps_spring_short, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_shortest.kml", append = FALSE)
st_write(lcps_autumn_short, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_shortest.shp", append = FALSE)
st_write(lcps_autumn_short, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_shortest.kml", append = FALSE)
st_write(lcps_spring_short_rotated, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Spring_shortest_rotated.kml", append = FALSE)
st_write(lcps_autumn_short_rotated, "~/Dropbox/BRC/Research/Spring Paper/Map/LCPaths_Autumn_shortest_rotated.kml", append = FALSE)

ggplot() +
  # geom_sf(data = st_jitter(lcps_autumn_short) %>% slice_sample(n = 100), alpha = 1/25, color = "yellow") +
  with_outer_glow(
    geom_sf(data = st_jitter(lcps_autumn_short), alpha = 1/25, color = "yellow"),
    colour = "yellow",
    sigma = 5) +
  # geom_sf(data = st_jitter(lcps_spring_short) %>% slice_sample(n = 100), alpha = 1/25, color = "white") +
  with_outer_glow(
    geom_sf(data = st_jitter(lcps_spring_short), alpha = 1/25, color = "white"),
    colour = "white",
    sigma = 5) +
  dark_theme_gray()

ggplot() +
  # geom_sf(data = st_jitter(lcps_autumn_short) %>% slice_sample(n = 100), alpha = 1/25, color = "yellow") +
  with_outer_glow(
    geom_sf(data = st_jitter(lcps_autumn_short_rotated), alpha = 1/25, color = "yellow"),
    colour = "yellow",
    sigma = 5) +
  # geom_sf(data = st_jitter(lcps_spring_short) %>% slice_sample(n = 100), alpha = 1/25, color = "white") +
  with_outer_glow(
    geom_sf(data = st_jitter(lcps_spring_short_rotated), alpha = 1/25, color = "white"),
    colour = "white",
    sigma = 5) +
  dark_theme_gray()

# lcps <- do.call(rbind, lcps_list)


# dem
# dem[is.na(dem)] <- 0
dem[dem < 1] <- NA

writeRaster(dem, "~/Dropbox/BRC/Research/Spring Paper/Map/Georgia_DEM_NonNeg.tif")

plot(dem)

dim(dem)

gradient <- seq(from = max(as.matrix(dem), na.rm = TRUE), to = min(as.matrix(dem), na.rm = TRUE), length.out = 800)

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

m <- rep.col(gradient, 1200)

d2 <- rast(m, extent = ext(dem))

writeRaster(dem + d2, "~/Dropbox/BRC/Research/Spring Paper/Map/DEM_Autumn.tif", overwrite = TRUE)


gradient_spring <- seq(to = max(as.matrix(dem), na.rm = TRUE), from = min(as.matrix(dem), na.rm = TRUE), length.out = 800)

m_spring <- rep.col(gradient_spring, 1200)

d2_spring <- rast(m_spring, extent = ext(dem))

writeRaster(dem + d2_spring, "~/Dropbox/BRC/Research/Spring Paper/Map/DEM_Spring.tif", overwrite = TRUE)
