wind_wave_height_new <- function(fetch, depths = list(), distances = list(), wind_speed, wind_proportion){
  # Check if the input list lengths are consistent: depths should have one more element than distances
  if(length(depths) != length(distances) + 1) {
    stop("Distance vector must be shorter than depth vectors by length of one")
  }

  # If there is only one or zero depths, return default wave characteristics of zero
  if(length(depths) <= 1) {
    df <- data.frame(
      wave_height_final = 0,
      WEI = 0,
      wave_period = 0,
      wave_number = 0,
      celerity_final = 0,
      nnumber_final = 0
    )
    return(df)
  }

  # ---- Constants ----
  G <- 9.80665 # acceleration due to gravity (m/s^2)
  A_goda <- 0.17 # empirical constant for Goda breaking height
  Con <- 1.56 * 1025 / 8 # (1025 kg/m^3 seawater density, 1.56 = G/(2*pi))
  alpha <- 1.909 # beach slope angle in degrees for breaking height
  avg_k <- 0 # initialize cumulative wave number


  # ---- Initialize vectors to store intermediate results ----
  crit_ratio <- double(length(depths))     # wave height to wind speed squared ratio
  wave_height <- double(length(depths))    # wave height at each segment
  breaking_flag <- double(length(depths))  # flag for wave breaking
  celerity <- double(length(depths))       # wave speed
  wave_number <- double(length(depths))    # wave number (k)
  group_velocity <- double(length(depths)) # group velocity
  wave_length <- double(length(depths))    # wave length (L)

  # ---- Calculate initial deep-water wave characteristics ----
  gf <- (0.077 * (G * fetch / wind_speed^2)^0.25)
  deep_wave_period <- 1.2 * tanh(gf) * 2 * pi * wind_speed / G
  deep_wave_length <- G * deep_wave_period^2 / (2 * pi)
  avg_depth <- mean(depths, na.rm = TRUE)

  # Determine whether we are in deep or shallow water based on depth-to-wavelength ratio
  depth_WL_ratio <- avg_depth / deep_wave_length
  if(depth_WL_ratio >= 2) {
    # Deep water condition
    gf <- (0.077 * (G * fetch / wind_speed^2)^0.25)
    wave_period <- (1.2 * tanh(gf) * wind_speed * 2 * pi / G)
  } else {
    # Shallow water condition
    gf <- (0.077 * (G * fetch / wind_speed^2)^0.25)
    gd <- (0.833 * (G * avg_depth / wind_speed^2)^0.375)
    wave_period <- (1.2 * tanh(gd) * tanh(gf / tanh(gd)) * wind_speed * 2 * pi / G)
  }

  gd <- (0.833 * (G * avg_depth / wind_speed^2)^0.375)

  # ---- Calculate wave number and wave length in average depth ----
  ko <- 2 * pi / deep_wave_length
  k1 <- newtonk(h = avg_depth, ko = ko)
  l1 <- 2 * pi / k1

  wave_number[1] <- newtonk(avg_depth, ko)
  wave_length[1] <- 2 * pi / wave_number[1]

  # ---- Calculate wave celerity ----
  if(avg_depth / l1 < 0.05){
    # Shallow water celerity approximation
    celer1 <- sqrt(G * avg_depth)
    celerity[1] <- celer1
  } else {
    # General celerity formula
    celer1 <- sqrt(G / k1 * tanh(avg_depth / k1))
    celerity[1] <- celer1
  }

  # ---- Calculate group velocity and initialize wave height ----
  n1 <- 0.5 * (1 + 2 * k1 * avg_depth / sinh(2 * k1 * avg_depth))
  group_velocity[1] <- n1

  gf <- (0.0125 * (G * distances[1] / wind_speed^2)^0.42)
  wave_height[1] <- 0.283 * tanh(gf) * wind_speed^2 / G

  k <- 0  # factor used in adjusting fetch after first segment

  # ---- Loop through segments to grow and transform wave ---
  for(j in 2:length(depths)){
    if(depths[j] <= 0){
      # Handle invalid depths by setting a small default value
      gf <- (0.0125 * (G * distances[j-1] / wind_speed^2)^0.42)
      wave_height[j] <- 0.283 * tanh(gf) * wind_speed^2 / G
      depths[j] <- 0.001
    }

    # ---- Estimate equivalent fetch from previous wave height ----
    # get the fetch equivalent to the wave ht and use that as a fetch
    # also check to make sure tanh doesn't lie in the singularity based on H/U^2 values
    crit_ratio[j] <- wave_height[j-1] / (wind_speed^2)
    if(crit_ratio[j] < 0.0288){
      new_fetch <- round(wind_speed^2 / G * (atanh(G * crit_ratio[j] / 0.283) / 0.0125)^(1 / 0.42), 2)
    } else {
      new_fetch <- round(wind_speed^2 / G * (atanh(G * 0.0288 / 0.283) / 0.0125)^(1 / 0.42), 2)
    }

    # ---- Calculate wave properties in current segment ----
    k2 <- newtonk(depths[j], ko)
    wave_number[j] <- k2
    l2 <- 2 * pi / k2
    wave_length[j] <- l2

    # ---- Calculate celerity ----
    if(depths[j] / l2 < 0.05){
      celer2 <- sqrt(G * depths[j])
    } else {
      celer2 <- sqrt(G / k2 * tanh(depths[j] / k2))
    }
    celerity[j] <- celer2

    # Calculate wave group velocity
    n2 <- 1 / 2 * (1 + 2 * k2 * depths[j] / sinh(2 * k2 * depths[j]))
    group_velocity[j] <- n2

    # ---- Shoaling and wave growth ----
    if(depths[j] / deep_wave_length < 0.5){
      # cat('j =', j, 'SHALLOW WAVE GROWING \n')
      # wave grows according to CERC curve
      gf <- (0.0125 * (G * (new_fetch + k * distances[j-1]) / wind_speed ^ 2) ^ 0.42)
      wave_height[j] <- (0.283 * tanh(gf) * wind_speed ^ 2 / G)

      # shoaling coefficient to describe the local wave height relative to the deep-water wave height
      shoaling_coef <- ((n1 * celer1) / (n2 * celer2)) ^ 0.5
      shoaling_coef <- ((group_velocity[j-1] * celerity[j-1]) / (group_velocity[j] * celerity[j])) ^ 0.5
      # shoaling_coef <- (n1 * celer1 / (n2 * celer2)) ^ 0.5
      # cat('shoaling coef = ', round(shoaling_coef, 3), ', wave height prior = ', round(wave_height[j], 3), '')

      wave_height[j] <- shoaling_coef * wave_height[j]
      # cat('wave height post = ', round(wave_height[j], 3), '\n')

      # Goda's Formula
      breaking_ht <- deep_wave_length * A_goda * (1 - exp(-1.5 * pi * depths[j] / deep_wave_length * (1 + 15 * (tan(alpha * pi / 180) ^ 1.333))))

      if(breaking_ht > 0.7 * depths[j]){
        breaking_ht <- 0.7 * depths[j]
      }
      # checking for breaking
      if(wave_height[j] > breaking_ht ){
        # breaking occurred
        wave_height[j] <- breaking_ht
        breaking_flag[j] <- 1
        # cat('breaking occured, breaking ht:', round(breaking_ht, 3), '')
      }
      # cat('depth =', round(depths[j], 3), 'wave height:', round(wave_height[j], 3), '\n')
    } else{
      # cat('j =', j, 'DEEP WAVE GROWING; ')
      # wave grows according to CERC curve
      gf <- (0.0125 * (G * (new_fetch + distances[j-1]) / wind_speed ^ 2) ^ 0.42)
      wave_height[j] <- (0.283 * tanh(gf) * wind_speed ^ 2 / G)

      # # shoaling coefficient to describe the local wave height relative to the deep-water wave height
      # shoaling_coef <- ((n1 * celer1)/(n2 * celer2))^0.5
      # # shoaling_coef <- (n1 * celer1 / (n2 * celer2)) ^ 0.5
      # cat('\nshoaling coef = ', round(shoaling_coef, 3), ', wave height prior = ', round(wave_height[j], 3), '')
      #
      # wave_height[j] <- shoaling_coef * wave_height[j]
      # cat('wave height post = ', round(wave_height[j], 3), '\n')

      # deep wave breaking criteria
      if((2 * wave_height[j] / deep_wave_length) > 0.3){
        # cat('Wave breaks ')
        wave_height[j] <- deep_wave_length / 2 * 0.3
      }
      # cat('depth =', round(depths[j], 3), 'wave height:', round(wave_height[j], 3), '\n')
    }

    k <- 1
    avg_k <- avg_k + k1
    n1 <- n2
    celer1 <- celer2
    k1 <- k2
    l1 <- l2

  }
  wave_height_final <- round(wave_height[j], 5)
  wave_number <- avg_k / j
  WEI <- Con * G * wave_height_final ^ 2 * wave_period^2 * tanh(wave_number * depths[j])

  df <- data.frame(
    # fetch,
    # wind_speed,
    # wind_proportion,
    wave_height_final = wave_height_final,
    WEI = WEI,
    wave_period = wave_period,
    wave_number = wave_number,
    celerity_final = celer2,
    nnumber_final = n2
  )

  # df$heights = list(round(wave_height, 5))
  # df$depths = list(depths)
  # df$distances = list(distances)


  return(df)
}
