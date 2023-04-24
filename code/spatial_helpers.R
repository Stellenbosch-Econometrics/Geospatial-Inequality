
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


# This function finds all cells of the grid in the vicinity of vectors or lon and lat coordinates
# using a radius of km * times/2 kilometers. It is necessary because measures like the IWI and RWI are sparse and 
# only available in populated places. Interpolated inequality estimates only makes sense in areas where we have multiple wealth estimates. 
# The function determines the degree-equivalent of a step-size in km, and then traces out cells with a radius of times/2 * km around a set of lat lon points (wealth estimates)
# E.g. if I have a 1km grid, and I choose a step-size of 0.5km (should be about half of the grid size to be safe), and times = 30, then 
# the function will find all 1km grid cells around my set of points with a radius of 30/2 * 0.5 = 7.5km (sampling the grid in steps of 0.5km which should get all relevant cells).
# Note that the sampling here uses a "square" rather than "round" search radius, but this does not matter in light of the gridded_distance_interpol() function 
# which will sort out cells that are too far from any given set of observed points. The important thing is that we fetch all potentially eligible cells with this function 
upsample_cells <- function(dggrid, lon, lat, km, times) {
  degrees = km / (40075.017 / 360)  # Gets the degree-distance of the kms at the equator
  cells = dgGEO_to_SEQNUM(dggrid, lon, lat)$seqnum
  init = -degrees * times / 2
  lat2 = lat + init; lon2 = lon + init
  for (j in seq_len(times)) {
    for (k in seq_len(times)) {
      lon2 %+=% degrees
      cells = funique(c(cells, dgGEO_to_SEQNUM(dggrid, lon2, lat2)$seqnum))
    }
    lat2 %+=% degrees
    lon2 = lon + init
    cells = funique(c(cells, dgGEO_to_SEQNUM(dggrid, lon2, lat2)$seqnum))
  }
  return(cells)
}

# This function performs spatial interpolation with arbitrary functions that allow inverse-distance-weights 
# (such as a function to calculate weighted GINI coefficeints)
# It has arguments: 
# - degree_grid: lat lon degree grid to limit the sampling region for finding points within a threshold distance (used to increase computational efficiency)
# - y: the data of interest (= wealth estimate) 
# - points_s2: cordinates of those measurements (y) as a apoint geometry created with s2::s2_lnglat(lon, lat)
# - cells_s2: centroids of cells to calculate interpolated measures for
# - FUN: statistic of interest (e.g. w_gini()) that expects a weights vector as the second argument: used to apply inverse-distance weights.
# - thresh: threshold in m under which surrounding points will be considered for a distance-weighted computation
# - step: stepsize in degrees, corresponding to the stepsize chosen to create degree_grid
# - tol: tolerance in degrees around the sampling region given by degree_grid. It should be greater than thresh / (40075017 / 360)
# - w: an optional additional weight vector (such as population), that is multiplied with the inverse-distance weights
# - ...: further arguments to FUN
# - verbose: TRUE prints lon an lat coordinates of the degree_grid as the function loops through it
#
# The function returns a vector of interpolated estimates for each cell corresponding to cells_s2. 
gridded_distance_interpol <- function(degree_grid, y, points_s2, cells_s2, FUN, thresh, 
                                      step = 1L, tol = 0.1, w = NULL, ..., verbose = FALSE) {
  if(length(y) != length(points_s2)) stop("length(y) must match length(points_s2)")
  degrees = qM(degree_grid)
  if(!identical(colnames(degrees), c("lon", "lat"))) stop("degree_grid must be a data frame with columns named 'lon' and 'lat'")
  wnull = is.null(w)
  if(!wnull && length(w) != length(y)) stop("length(w) must match length(y)")
  
  px = s2_x(points_s2); py = s2_y(points_s2); cx = s2_x(cells_s2); cy = s2_y(cells_s2)
  res = alloc(NA_real_, length(cx))
  
  for(d in mrtl(degrees)) {
    lon_tol = tol / cos(d[2L] * pi/180) # Spherical correction of tolerance
    if(verbose) cat("lon =", d[1L], " lat =", d[2L], fill = TRUE)
    which_cells = which(between(cx, d[1L], d[1L]+step) & between(cy, d[2L], d[2L]+step))
    which_points = which(between(px, d[1L]-lon_tol, d[1L]+(lon_tol+step)) & between(py, d[2L]-tol, d[2L]+(tol+step)))
    if(length(which_cells) == 0L || length(which_points) < 2L) next
    points_s2_subset = points_s2[which_points]
    y_subset = y[which_points]
    if(wnull) {
      for(i in which_cells) {
        disti = s2_distance(points_s2_subset, cells_s2[i])
        ind = which(disti < thresh)
        if(length(ind) > 1L) res[i] = FUN(y_subset[ind], thresh - disti[ind], ...)
      }
    } else {
      w_subset = w[which_points]
      for(i in which_cells) {
        disti = s2_distance(points_s2_subset, cells_s2[i])
        ind = which(disti < thresh)
        if(length(ind) > 1L) res[i] = FUN(y_subset[ind], w_subset[ind] * (thresh - disti[ind])/thresh, ...) 
      }
    }
  }
  return(res)
}
