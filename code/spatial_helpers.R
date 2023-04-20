
# Haversine distance between different geographic points in m (spherical earth)
haversine <- function(lon1, lat1, lon2, lat2, r = 6378137) { # r = radius of the earth in m
  lat1 <- lat1 * pi/180
  lon1 <- lon1 * pi/180
  lat2 <- lat2 * pi/180
  lon2 <- lon2 * pi/180
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- atan2(sqrt(a), sqrt(1 - a))
  return(c * (2 * r))
}

# Check: geosphere::distHaversine(c(1, 2), c(3, 4)) against haversine(1,2,3,4)

# Optimized version of the same function
haversine_opt <- function(lon1, lat1, lon2, lat2, r = 6378137) {
  lon1 <- lon1 * pi/180
  lat1 <- lat1 * pi/180
  lon2 <- lon2 * pi/180
  lat2 <- lat2 * pi/180
  dlat <- (lat2 - lat1) %/=% 2
  sin_dlat2 <- sin(dlat)^2
  lon2 %-=% lon1 %/=% 2
  sin_dlon2 <- sin(lon2)^2
  a <- sin_dlat2 %+=% (cos(lat1) %*=% cos(lat2) %*=% sin_dlon2)
  a %/=% (1 - a) 
  res <- atan(sqrt(a)) %*=% (2 * r)
  return(res)
}

# Gives the ratio of the distance between two adjacent longitudes to their distance at the equator.
haversine_lon1lat <- function(lat = 0) {
  a0 <- sin(pi/360)^2
  a <- cos(lat * pi/180)^2 * a0
  atan2(sqrt(a), sqrt(1 - a)) / 
  atan2(sqrt(a0), sqrt(1 - a0))    
  # Or: 
  # atan(sqrt(a/(1 - a))) / 
  # atan(sqrt(a0/(1 - a0)))
}

# Note that this is approximately equal to cos(lat * pi/180), but differs very slightly, especially towards the poles..
# See: ts.plot(cos(0:90 * pi/180)/haversine_lon1lat(0:90))

# Optimized version of the same function
haversine_lon1lat_opt <- function(lat = 0) {
  a0 <- sin(pi/360)^2
  denom <- atan2(sqrt(a0), sqrt(1 - a0))
  tmp <- cos(lat * pi/180)
  tmp %*=% tmp %*=% a0
  tmp %/=% (1 - tmp) # a/(1-a) = 1 / (1/a - 1) 
  res <- atan(sqrt(tmp)) %/=% denom 
  res
}

# Transforms longitudes so that they measure distance
trans_lon <- function(lon, lat) {
  lon * haversine_lon1lat_opt(lat)
}

# This function efficiently transforms all points to the centroid of a spatial grid with 
# km kilometers extent, e.g. km = 1 produces a 1km^2 spatial grid. The function has 
# an option to round the degrees-kms to any number of digits to facilitate matching geospatial
# datasets using the grid easier. A precision of 6 digits should be available on every system. 
round_to_kms_exact <- function(lon, lat, km, round = TRUE, digits = 6) {
  degrees <- km / (40075.017 / 360)               # Gets the degree-distance of the kms at the equator
  if(round) div <- round(degrees, digits)         # Round to precision
  res_lat <- TRA(lat, div, "-%%") %+=% (div/2)    # This transforms the data to the grid centroid... 
  scale_lat <- haversine_lon1lat_opt(res_lat)     # Scale factor based on grid centroid
  res_lon <- setTRA(lon * scale_lat, div, "-%%") %+=% (div/2) %/=% scale_lat  
  return(list(lon = res_lon, lat = res_lat))
}

# Same function but using the cos(lat * pi/180) approximation, and thus faster...
round_to_kms_fast <- function(lon, lat, km, round = TRUE, digits = 6) {
  degrees <- km / (40075.017 / 360)               # Gets the degree-distance of the kms at the equator
  if(round) div <- round(degrees, digits)         # Round to precision
  res_lat <- TRA(lat, div, "-%%") %+=% (div/2)    # This transforms the data to the grid centroid... 
  scale_lat <- cos(res_lat * pi/180)              # Approx. scale factor based on grid centroid
  res_lon <- setTRA(lon * scale_lat, div, "-%%") %+=% (div/2) %/=% scale_lat  
  return(list(lon = res_lon, lat = res_lat))
}

# This linearly transforms a variable to a 0-1 scale
ltr <- function(x) {
  r <- frange(x)
  (x - as.double(r[1L])) %/=% (r[2L] - r[1L]) 
}


