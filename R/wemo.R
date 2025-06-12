WEMO <- function(site_point,
                 bathy,
                 shoreline,
                 wind_data = data.frame(time = NA, year = NA, month = NA, wind_direction = NA, wind_speed = NA),
                 max_fetch = 10000,
                 wave_iteration_distance = 100,
                 wind_percentile = 0.95,
                 water_level = 0,
                 depths_or_elev = "elev") {

  #### check input data ####
  # check formats
  if(!inherits(site_point, 'sf')){
    stop("Point isn't an sf object")
  }
  # check names of variables
  if(!all(c("site_name", "Id") %in% names(site_point))){
    stop("site_point needs variables named Id and site_name")
  }

  if(!inherits(shoreline, 'sf')){
    stop("Shoreline isn't an sf object")
  }
  #
  # if(!inherits(bathy, 'RasterLayer')){
  #   stop("bathy isn't a RasterLayer object")
  # }

  # check matching CRS
  if((sf::st_crs(site_point) != sf::st_crs(shoreline))){
    cat("The CRS for shoreline and site_point differ; transforming shoreline CRS to match site_point... ")
    shoreline <- sf::st_transform(shoreline, sf::st_crs(site_point))
    if(sf::st_crs(site_point) == sf::st_crs(shoreline)){
      cat("SUCCESS!   \n")
    } else{
      cat("... unable to transform shoreline to CRS of site_point.   \n")
      stop("CRS of site_point doesn't match CRS of shoreline")
    }
  }

  if(sf::st_crs(site_point) != sf::st_crs(bathy)){
    cat("The CRS for bathy and site_point differ; transforming bathy CRS to match site_point... ")
    terra::crs(bathy) <- terra::crs(site_point)
    if(sf::st_crs(site_point) == sf::st_crs(bathy)){
      cat("SUCCESS!   \n")
    } else{
      cat(" unable to transform bathy to CRS of site_point.   \n")
      stop("CRS of bathy doesn't match CRS of site_point")
    }
  }

  # crop bathy and shoreline to the bounding box
  crop_bbox <- create_bounding_box(site_point, max_fetch * 1.1)
  shoreline <- sf::st_crop(shoreline, crop_bbox)
  # bathy <- crop(bathy, crop_bbox)

  #### WIND ####
  # summarize wind based on what the user inputs
  wind_data <- wind_data %>%
    dplyr::mutate(wind_direction = ifelse(wind_direction == 360, 0, wind_direction)) %>%
    dplyr::filter(wind_speed != 0)

  if(wind_percentile == "mean"){
    wind_data_summary <- wind_data %>%
      dplyr::mutate(directions = ifelse(wind_direction == 360, 0, wind_direction)) %>%
      dplyr::group_by(directions) %>%
      dplyr::summarize(
        n = dplyr::n(),
        proportion = 100 * n/nrow(.),
        speed = mean(wind_speed, na.rm = T)) %>%
      dplyr::ungroup()
  }else{
    wind_data_summary <- wind_data %>%
      dplyr::mutate(directions = ifelse(wind_direction == 360, 0, wind_direction)) %>%
      dplyr::filter(wind_speed != 0) %>%
      dplyr::group_by(directions) %>%
      dplyr::summarize(
        n = dplyr::n(),
        proportion = 100 * n/nrow(.),
        speed = stats::quantile(wind_speed, wind_percentile, na.rm = T)) %>%
      dplyr::ungroup()
  }

  #### FETCH ####
  ## calculate fetch
  n_directions <- wind_data_summary %>% dplyr::filter(!is.na(directions)) %>% dplyr::pull(directions) %>% unique %>% length

  site_point <- site_point %>%
    dplyr::arrange(site_name)

  site_names <- data.frame(site_name = site_point$site_name, key = NA)
  site_names$key <- as.character(1:nrow(site_names))

  cat("Calculating Fetch...    \n")
  fetch <- windfetch::windfetch(shoreline, site_point, progress_bar = T, max_dist = max_fetch / 1000, n_directions = n_directions/4, quiet = T) %>%
    windfetch::as_sf() %>%
    dplyr::mutate(key = stringr::str_remove(site_name, "Site " )) %>%
    dplyr::select(-site_name)

  fetch <- fetch %>%
    dplyr::left_join(site_names, by = 'key') %>%
    dplyr::select(site_name, directions, quadrant, fetch, geometry)

  ## calculate effective fetch
  cat("Calculating Effective Fetch...    \n")
  efetch <- effective_fetch(fetch, n_directions)

  fetch_all <- dplyr::left_join(efetch,
                                fetch %>% as.data.frame() %>% dplyr::select(-geometry),
                                by = c('site_name', 'directions', 'quadrant')) %>%
    dplyr::mutate(depths = NA)

  #### BATHY ####
  ## grab depths along each fetch ray at sampling distance
  cat("Grabbing Bathy Data...    \n")
  fetch_all$depths <- lapply(1:nrow(fetch_all), function(i){

    if(as.numeric(fetch_all[[i, "fetch"]]) > (wave_iteration_distance)){
      d <- get_bathy_depth(bathy_raster = bathy, fetch_ray =  fetch_all[i, ], depth_sample_distance = wave_iteration_distance) %>%
        .[,2]
    } else {
      d = 0
    }

    # add the tide height to extracted bathy values
    if(depths_or_elev == 'depths'){
      d <- d + water_level
    }else if(depths_or_elev == 'elev'){
      d <- -1*d + water_level
    }
    d <- rev(d)
    return(d)
  })
  # join the wind data to the bathy data
  fetch_all <- dplyr::left_join(fetch_all, wind_data_summary, by = "directions")

  #### Build Wave ####
  # Calculate wind wave height along the fetch rays
  cat("Building Wind Waves...    \n")
  wemo_details <- lapply(1:nrow(fetch_all), function(i){
    height_df <- wind_wave_height(fetch = as.numeric(fetch_all$efetch[i]),
                                  depths = as.numeric(fetch_all$depths[[i]]),
                                  sample_distance = wave_iteration_distance,
                                  wind_speed = as.numeric(fetch_all$speed[i]),
                                  wind_proportion = as.numeric(fetch_all$proportion[i]))

    height_df <- height_df %>%
      dplyr::mutate(site_name = fetch_all$site_name[i], directions = fetch_all$directions[i], .before = 1)

    return(height_df)
  })

  wemo_details <- dplyr::bind_rows(wemo_details)

  # heights_and_depths <- data.frame(
  #   wemo_details$depths,
  #   wemo_details$heights
  # )

  # wemo_details <- wemo_details %>%
  #   mutate(heights = NULL,
  #          depths = NULL)
  #

  wemo_output <- wemo_details %>%
    dplyr::mutate(RWE = WEI * wind_proportion/100,
           site_name = as.character(site_name)) %>%
    dplyr::group_by(site_name) %>%
    dplyr::reframe(
      RWE = sum(RWE, na.rm = T),
      max_wave_height = max(wave_height_final, na.rm = T),
      direction_of_max_wave = paste(list(directions[which(wave_height_final == max_wave_height)])),
      avg_wave_height = mean(wave_height_final, na.rm = T),
      avg_wave_period = mean(wave_period, na.rm = T),
      max_wave_period = max(wave_period, na.rm = T)
    )

  # if(velocity_mode == T){
  #   if(ov_method == T){
  #
  #   }
  # }

  return(
    list(wemo_output = wemo_output,
         wemo_details = wemo_details,
         # heights_and_depths,
         fetch = fetch,
         efetch = efetch,
         wind_data_summary = wind_data_summary)
  )
}
