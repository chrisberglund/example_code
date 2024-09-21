library(dplyr)
library(aniMotum)
library(pathroutr)

regularize_tracks <- function(positions, time_interval, vmax, what = "predicted") {
  movement_data <- positions %>%
    arrange(id, datetime) %>%
    mutate(date = datetime) %>%
    group_by(id) %>%
    nest()
  
  future::plan(future::multisession, workers = 4)
  fits <- bind_rows(movement_data %>%
                      furrr::future_pmap(\(id, data) aniMotum::fit_ssm(
                        x = cbind(id, data),
                        model = "crw",
                        vmax = vmax,
                        time.step = time_interval
                      ), .options = furrr::furrr_options(seed = 123)))
  future::plan(future::sequential)
  
  predicted <- aniMotum::grab(fits, what = what) %>%
    dplyr::select(id,lon, lat, date)
  
  return(predicted)
}

reroute_paths <- function(positions) {
  coastline <- st_read("data/antarctica_coastline/add_coastline_high_res_polygon_v7_8.gpkg",
                       layer = "add_coastline_high_res_polygon_v7_8")
  
  # Subset coastline polygons by a convex hull around a 50km buffer of all points
  positions_sf <- positions %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326), remove = FALSE) %>%
    st_transform(3031)
  
  
  boundary_region <- positions_sf %>%
    st_buffer(dist = 70000) %>%
    st_union() %>%
    st_convex_hull() %>%
    st_intersection(coastline) %>%
    st_collection_extract("POLYGON") %>%
    st_sf()
  
  vis_graph <- pathroutr::prt_visgraph(boundary_region)
  
  ids <- distinct(positions_sf, id) %>% pull()
  if (exists("new_tracks")) rm(new_tracks)
  pb <- txtProgressBar(1, length(ids), style = 3)
  for (i in 1:length(ids)) {
    track <- filter(positions_sf, id == ids[i])
    track <- prt_trim(track, boundary_region)
    rerouted_pts <- prt_reroute(track, boundary_region, vis_graph, blend = TRUE)
    pts <- prt_update_points(rerouted_pts, track)
    
    if (exists("new_tracks")) {
      new_tracks <- bind_rows(new_tracks, pts)
    } else {
      new_tracks <- pts
    }
    setTxtProgressBar(pb, value = i)
  }
  close(pb)
  new_tracks <- new_tracks %>% 
    rename(datetime = date) %>%
    dplyr::select(-fid) %>% 
    st_drop_geometry()
  
  return(new_tracks)
}

prepare_tracks <- function(positions, time_interval, vmax, what = "predicted") {
  tracks <- regularize_tracks(positions, time_interval, vmax, what)
  tracks <- reroute_paths(tracks)
  return(tracks)
} 

