#' Calculate Speed Likelihood for Particle Filter (Simplified)
#'
#' Calculates a weight for each proposed particle based on the speed required
#' to reach it from the mean location of the previous particle cloud.
#' Uses a Gaussian likelihood around optimal speeds.
#'
#' @param proposed_particles_dt data.table of potential particles for the current step.
#'                              Must include 'lon', 'lat', 'dtime', 'prob.dry'.
#' @param prev_cloud_summary list containing summary of the previous cloud:
#'                           'mean_lon', 'mean_lat', 'mean_dtime'.
#' @param speed.dry Speed parameters (optimal, sd, max) for dry state.
#' @param speed.wet Speed parameters (optimal, sd, max) for wet state.
#' @param distance.method 'ellipsoid' or 'spherical'.
#' @return A numeric vector of speed weights/likelihoods for proposed_particles_dt.
#' @import data.table
#' @import sp
#' @import sf # if distance.method == 'spherical'
#' @importFrom stats dnorm

calculate_speed_likelihood_PF <- function(proposed_particles_dt, prev_cloud_summary,
                                          speed.dry, speed.wet, distance.method = "ellipsoid") {

    if (nrow(proposed_particles_dt) == 0) return(numeric(0))
        if (is.null(prev_cloud_summary) || is.na(prev_cloud_summary$mean_lon)) {
            # Handle the very first step where there's no previous cloud
            # Assign equal weight (or weight based only on environmental data later)
            return(rep(1, nrow(proposed_particles_dt)))
        }

    # Calculate distance from each proposed particle to the previous mean location
    proposed_coords <- as.matrix(proposed_particles_dt[, .(lon, lat)])
    prev_mean_coords <- matrix(c(prev_cloud_summary$mean_lon, prev_cloud_summary$mean_lat), ncol = 2)

    if (distance.method == "ellipsoid") {
        dists_m <- sp::spDistsN1(pts = proposed_coords, pt = prev_mean_coords, longlat = TRUE) * 1000
    } else { # spherical
        suppressMessages(sf_use_s2(TRUE))
        proposed_sf <- st_as_sf(proposed_particles_dt, coords = c("lon", "lat"), crs = 4326)
        prev_mean_sf <- st_as_sf(data.frame(lon=prev_mean_coords[1], lat=prev_mean_coords[2]), coords=c("lon","lat"), crs=4326)
        dists_m <- as.numeric(st_distance(proposed_sf, prev_mean_sf))
        # suppressMessages(sf_use_s2(FALSE))
    }

    # Calculate time difference
    time_diff_sec <- abs(as.numeric(difftime(proposed_particles_dt$dtime,
                                              prev_cloud_summary$mean_dtime, units = "secs")))

    # Calculate speed
    speeds_ms <- ifelse(time_diff_sec > 0, dists_m / time_diff_sec, 0)

    # --- Calculate weight based on speed ---
    # Using a Gaussian likelihood around the optimal speed, adjusted by prob.dry
    # Determine optimal speed and SD based on prob.dry (assuming prob.dry exists)
    if (!"prob.dry" %in% names(proposed_particles_dt)) {
        warning("Column 'prob.dry' not found in proposed particles. Assuming always wet for speed calculation.")
        proposed_particles_dt[, prob.dry := 0] # Assume wet if no data
    }


    optimal_speed <- proposed_particles_dt$prob.dry * speed.dry[1] + (1 - proposed_particles_dt$prob.dry) * speed.wet[1]
    sd_speed <- proposed_particles_dt$prob.dry * speed.dry[2] + (1 - proposed_particles_dt$prob.dry) * speed.wet[2]
    max_speed <- proposed_particles_dt$prob.dry * speed.dry[3] + (1 - proposed_particles_dt$prob.dry) * speed.wet[3]

    # Calculate Gaussian likelihood (weight)
    # Set weight to near zero if speed exceeds max speed
    speed_weights <- dnorm(speeds_ms, mean = optimal_speed, sd = sd_speed)
    speed_weights[speeds_ms > max_speed] <- 1e-9 # Penalize exceeding max speed heavily
    speed_weights[is.na(speed_weights)] <- 1e-9 # Handle potential NAs

    return(speed_weights)
}


#' Systematic Resampling
#'
#' Performs systematic resampling for particle filters.
#'
#' @param weights A numeric vector of particle weights (do not need to sum to 1).
#' @param N Number of particles to resample.
#' @return A numeric vector of indices of the resampled particles.

systematic_resample <- function(weights, N) {
    if (sum(weights, na.rm = TRUE) == 0) {
        # Handle case with all zero weights - resample uniformly
        warning("All particle weights are zero or NA. Resampling uniformly.")
        return(sample.int(length(weights), size = N, replace = TRUE))
    }
    # Normalize weights
    norm_weights <- weights / sum(weights, na.rm = TRUE)
    norm_weights[is.na(norm_weights)] <- 0 # Ensure NAs don't break cumsum

    # Cumulative sum
    cum_weights <- cumsum(norm_weights)

    # Generate starting point and subsequent points
    u1 <- runif(1) / N
    u <- u1 + (0:(N - 1)) / N

    # Find indices
    indices <- integer(N)
    j <- 1
    for (i in 1:N) {
        while (u[i] > cum_weights[j]) {
            j <- j + 1
        }
        indices[i] <- j
    }
    return(indices)
}


#' Calculate Particle Weights for Particle Filter
#'
#' Computes weights based on environmental factors and pre-computed speed likelihood.
#'
#' @param step_data_dt data.table for the current step particles. Must include
#'   environmental data (sat.sst, tag.sst etc.) and 'weight_land' if applicable.
#' @param speed_likelihood Pre-calculated speed likelihood/weight vector for these particles.
#' @param sst.sd Standard deviation for SST likelihood.
#' @param max.sst.diff Max allowable SST difference.
#' @param ice.conc.cutoff Max allowable ice concentration.
#' @param land.mask Logical, is land mask applied?
#' @return Original data.table with added weight columns (weight_sst, weight_ice, rel_weight).
#' @import data.table

calculate_particle_weights_PF <- function(step_data_dt, speed_likelihood, sst.sd,
                                          max.sst.diff, ice.conc.cutoff, land.mask) {

    # Add speed likelihood directly
    step_data_dt[, weight_speed := speed_likelihood]

    # --- Calculate other weights ---

    # SST Weight (using external function fun_weight_sst - assumed available)
    if ("tag.sst" %in% names(step_data_dt) && !all(is.na(step_data_dt$tag.sst))) {
        step_data_dt[, weight_sst := fun_weight_sst(.SD, y = sst.sd, z = max.sst.diff),
                    .SDcols = c("sst.diff")] # Pass necessary columns
        step_data_dt[is.na(weight_sst), weight_sst := 1] # Neutral weight if SST calc fails/NA
    } else {
        step_data_dt[, weight_sst := 1] # Neutral weight if no SST data
    }

    # Ice Weight (using external function fun_weight_ice - assumed available)
    step_data_dt[, weight_ice := fun_weight_ice(.SD, y = ice.conc.cutoff),
                .SDcols = c("sat.ice")] # Pass necessary columns
    step_data_dt[is.na(weight_ice), weight_ice := 1] # Neutral weight if Ice is NA

    # --- Combine weights ---
    step_data_dt[, rel_weight := weight_speed * weight_sst * weight_ice]

    # Apply land mask weight if it was calculated
    if (!is.null(land.mask) && "weight_land" %in% names(step_data_dt)) {
        step_data_dt[, rel_weight := rel_weight * weight_land]
    }

    # Weights do NOT need to be normalized here; normalization happens before resampling.
    # Ensure no NA weights before resampling step
    step_data_dt[is.na(rel_weight), rel_weight := 0]

    return(step_data_dt)
}


#' Run Particle Filter Iterations
#'
#' Implements a particle filter using probGLS proposals.
#'
#' @param filtered_particles_dt Data.table from `filter_particle_clouds` (proposal distribution).
#' @param particle.number Number of particles to maintain in the filter.
#' @param tagging.location sf object for the start location.
#' @param act Wet/dry activity data.frame.
#' @param sensor Tag temperature sensor data.frame.
#' @param speed.dry Speed parameters (dry).
#' @param speed.wet Speed parameters (wet).
#' @param sst.sd SST likelihood SD.
#' @param max.sst.diff Max SST difference.
#' @param ice.conc.cutoff Max ice concentration.
#' @param land.mask Is land mask used?
#' @param wetdry.resolution Resolution of wet/dry sensor.
#' @param NOAA.OI.location Path to environmental data.
#' @param distance.method 'ellipsoid' or 'spherical'.
#' @param backward Run backwards?
#' @return A data.table containing the summary (mean/median) of the particle cloud at each step.
#' @import data.table
#' @import sf
#' @import sp
#' @importFrom stats median sd
#' @importFrom utils txtProgressBar setTxtProgressBar

run_particle_filter <- function(filtered_particles_dt, particle.number, tagging.location,
                                act, sensor, speed.dry, speed.wet, sst.sd,
                                max.sst.diff, ice.conc.cutoff, land.mask,
                                wetdry.resolution, NOAA.OI.location,
                                distance.method, backward) {

    # --- Initialization ---
    # Start with particles clustered at the tagging location
    # Convert tagging location sf to data.table
    tag_loc_coords <- st_coordinates(tagging.location)
    initial_particles <- data.table(
        lon = rnorm(particle.number, mean = tag_loc_coords[1], sd = 0.01), # Add small initial spread
        lat = rnorm(particle.number, mean = tag_loc_coords[2], sd = 0.01),
        dtime = tagging.location$dtime[1] # Use the start time
        # Add other necessary columns with default/NA values if needed downstream
    )
    # Ensure initial particles are within bounds/land mask if necessary (optional refinement)

    current_particles_dt <- initial_particles
    prev_cloud_summary <- list(mean_lon = mean(current_particles_dt$lon),
                               mean_lat = mean(current_particles_dt$lat),
                               mean_dtime = current_particles_dt$dtime[1])

    # Convert sensor data to data.table for faster lookups
    if (!is.null(sensor)) setDT(sensor)
    if (!is.null(act)) setDT(act) # Also convert activity data if used


    # Determine step order
    if (backward) {
        # Backward filtering is more complex, usually requiring modifications (e.g., particle smoothing)
        # Sticking to forward filtering for this implementation.
        if (backward) warning("Backward particle filtering not fully implemented, proceeding forward.")
        # steps <- sort(unique(filtered_particles_dt$step), decreasing = TRUE)
        steps <- sort(unique(filtered_particles_dt$step), decreasing = FALSE)
    } else {
        steps <- sort(unique(filtered_particles_dt$step), decreasing = FALSE)
    }

    # List to store summary results from each step
    results_summary_list <- vector("list", length = length(steps))
    names(results_summary_list) <- steps
    # Optional: Store full clouds if memory allows
    # full_clouds_list <- vector("list", length = length(steps))
    # names(full_clouds_list) <- steps

    progress_bar <- txtProgressBar(min = 0, max = length(steps), style = 3, char = '=')
    step_counter <- 0

    # --- Main Loop through Time Steps ---
    for (ts in steps) {
        step_counter <- step_counter + 1
        # Get the proposal distribution for this step
        potential_particles_t <- filtered_particles_dt[step == ts]

        if (nrow(potential_particles_t) == 0) {
            warning(paste("No proposal particles available for step", ts, "- skipping."))
            # Handle skipped step
            results_summary_list[[as.character(ts)]] <- data.table(step = ts, lon_mean=NA, lat_mean=NA, lon_median=NA, lat_median=NA, time_median=NA) # Store NA summary
            setTxtProgressBar(progress_bar, value = step_counter)
            next
        }

        # --- Prepare Environmental Data for Proposals ---
        # (This block is similar to the previous version)
        current_date <- potential_particles_t$date[1]
        current_year <- potential_particles_t$year[1]
        ls_files <- list.files(NOAA.OI.location)

        fname.sst <- paste0(NOAA.OI.location, '/', ls_files[grep(paste0('sst.day.mean.', current_year), ls_files)])[1]
        fname.err <- paste0(NOAA.OI.location, '/', ls_files[grep(paste0('sst.day.err.', current_year), ls_files)])[1]
        fname.ice <- paste0(NOAA.OI.location, '/', ls_files[grep(paste0('icec.day.mean.', current_year), ls_files)])[1]

        tryCatch({
            potential_particles_t[, sat.sst := load_NOAA_OISST_V2(FILE_NAME = fname.sst, LONS = lon, LATS = lat, DATE = current_date, extract.value = 'sst')]
            potential_particles_t[, sat.ice := load_NOAA_OISST_V2(FILE_NAME = fname.ice, LONS = lon, LATS = lat, DATE = current_date, extract.value = 'icec')]
            potential_particles_t[, sat.sst.err := load_NOAA_OISST_V2(FILE_NAME = fname.err, LONS = lon, LATS = lat, DATE = current_date, extract.value = 'err')]
        }, error = function(e) {
            warning(paste("PF: Could not load env data for date:", current_date, "Error:", e$message))
             potential_particles_t[, `:=`(sat.sst = NA_real_, sat.ice = NA_real_, sat.sst.err = NA_real_)]
        })
        potential_particles_t[is.na(sat.ice), sat.ice := 0]
        if (!(as.Date(current_date) >= as.Date("2012-08-11") & as.Date(current_date) <= as.Date("2012-08-16"))) {
           potential_particles_t[sat.ice > ice.conc.cutoff, sat.sst := NA]
        }
        
        # browser()
        
        current_jday <- potential_particles_t$jday_int[1]
        tag_temp_today <- NA_real_
        
        if(!is.null(sensor)) {
          sensor$jday <- as.numeric(julian(sensor$date))
        }

        if (nrow(sensor[jday == current_jday]) > 0) {
           tag_temp_today <- sensor[jday == current_jday, SST][1]
        }
        potential_particles_t[, tag.sst := tag_temp_today]
        potential_particles_t[, sst.diff := sat.sst - tag.sst]

        # --- Calculate Prob Dry ---
        if (!is.null(act) && !is.null(prev_cloud_summary)) {
            potential_particles_t[, prob.dry := 0.5]
        } else {
            potential_particles_t[, prob.dry := 1.0]
        }


        # --- Weighting Step ---
        # 1. Calculate Speed Likelihood
        speed_weights <- calculate_speed_likelihood_PF(
            proposed_particles_dt = potential_particles_t,
            prev_cloud_summary = prev_cloud_summary,
            speed.dry = speed.dry, speed.wet = speed.wet,
            distance.method = distance.method
        )

        # 2. Calculate Environmental Weights and Combine
        weighted_particles_t <- calculate_particle_weights_PF(
            step_data_dt = potential_particles_t,
            speed_likelihood = speed_weights,
            sst.sd = sst.sd, max.sst.diff = max.sst.diff,
            ice.conc.cutoff = ice.conc.cutoff, land.mask = land.mask
        )

        # --- Resampling Step ---
        resampling_indices <- systematic_resample(weights = weighted_particles_t$rel_weight,
                                                  N = particle.number)

        # Create the new particle cloud for this step
        current_particles_dt <- weighted_particles_t[resampling_indices, ]

        # --- Store Results & Update Summary ---
        if (nrow(current_particles_dt) > 0) {
            # Calculate summary statistics for the current cloud
            step_summary <- current_particles_dt[, .(
                step = ts,
                lon_mean = mean(lon), lat_mean = mean(lat),
                lon_median = median(lon), lat_median = median(lat),
                lon_sd = sd(lon), lat_sd = sd(lat), # Add uncertainty estimate
                time_median = median(dtime) # Store representative time
                # Add other relevant median/mean solar angle
                # median.solar.angle = median(solar.angle, na.rm=TRUE)
            )]
            results_summary_list[[as.character(ts)]] <- step_summary

            # Update summary for the *next* step's speed calculation
            prev_cloud_summary <- list(mean_lon = step_summary$lon_mean,
                                    mean_lat = step_summary$lat_mean,
                                    mean_dtime = step_summary$time_median)

            # Optional: Store full cloud
            # current_particles_dt[, step := ts] # Add step identifier
            # full_clouds_list[[as.character(ts)]] <- copy(current_particles_dt)

        } else {
            # Handle case where resampling failed (e.g., all weights zero)
            warning(paste("Resampling failed for step", ts, "- no particles generated."))
            results_summary_list[[as.character(ts)]] <- data.table(step = ts, lon_mean=NA, lat_mean=NA, lon_median=NA, lat_median=NA, time_median=NA) # Store NA summary
            prev_cloud_summary <- NULL # Reset summary so next step starts fresh if possible
        }

        setTxtProgressBar(progress_bar, value = step_counter)
    } # End loop through steps (ts)

    close(progress_bar)

    # Combine summary results from all steps
    final_summary_track <- rbindlist(results_summary_list, use.names = TRUE, fill = TRUE)
    final_summary_track <- final_summary_track[!is.na(lon_mean)] # Remove steps that failed

    # Convert times back to POSIXct if they became numeric
    if (is.numeric(final_summary_track$time_median)) {
      final_summary_track[, time_median := as.POSIXct(time_median, origin = "1970-01-01", tz = "UTC")]
    }

    # Optional: Combine full clouds if stored
    # all_particle_clouds <- rbindlist(full_clouds_list, use.names = TRUE, fill = TRUE)


    # Return the summary track (mean/median locations per step)
    # Or return the list containing both summary and full clouds if desired
    return(final_summary_track)
}


#' Particle Filter algorithm for geolocation data (probGLS based proposals)
#'
#' Applies a particle filter using probGLS generated proposals, incorporating
#' resampling and environmental weighting.
#'
#' @param particle.number Number of particles for the filter.
#' @param ... (Other parameters identical to prob_algorithm_refactored, EXCEPT iteration.number which is ignored)
#' @return A list containing: [1] Estimated track (summary statistics), [2] All proposal particles generated, [3] Input parameters, [4] Run time.
#' @import data.table
#' @import sf
#' @import sp
#' @export

prob_algorithm_PF <- function(
    particle.number             = 500, # PF usually needs more particles
    # iteration.number parameter is ignored here
    loess.quartile              = NULL,
    trn                         = trn,
    sensor                      = NULL,
    act                         = NULL,
    tagging.date                = NULL,
    retrieval.date              = NULL,
    tol                         = 0.08,
    tagging.location            = c(-36.816,-54.316),
    boundary.box                = c(-180,180,-90,90),
    sunrise.sd                  = c(1,1,0),
    sunset.sd                   = c(1,1,0),
    range.solar                 = c(-7,-1),
    speed.wet                   = c(1,1.3,5),
    speed.dry                   = c(12,6,45),
    sst.sd                      = 0.5,
    max.sst.diff                = 3,
    ice.conc.cutoff             = 1,
    wetdry.resolution           = 1,
    east.west.comp              = TRUE,
    land.mask                   = TRUE,
    med.sea                     = TRUE,
    black.sea                   = TRUE,
    baltic.sea                  = TRUE,
    caspian.sea                 = TRUE,
    backward                    = FALSE, # Note: Backward PF not implemented
    distance.method             = "ellipsoid",
    NOAA.OI.location            = getwd()
){

    start.time <- Sys.time()
    oldw <- getOption("warn")
    options(warn = -1)

    # --- Input Validation & Setup ---
    # (Similar to prob_algorithm_refactored)
    if (is.null(trn)) stop("Twilight data 'trn' must be provided.")
    if (is.null(tagging.date)) stop("Tagging date must be provided.")
    if (is.null(retrieval.date)) stop("Retrieval date must be provided.")
    # Check packages...

    # Store input parameters (iteration.number removed/ignored)
    model.input <- data.frame(parameter=c('particle.number','loess.quartile','tagging.location',
                                          'tagging.date','retrieval.date','sunrise.sd','sunset.sd','range.solar','speed.wet',
                                          'speed.dry','sst.sd','max.sst.diff','ice.conc.cutoff','boundary.box','med.sea','black.sea',
                                          'baltic.sea','caspian.sea','east.west.comp','wetdry.resolution','NOAA.OI.location','backward','sensor.data','act.data','land.mask','distance.method', 'method'),
                              chosen=c(paste(particle.number,collapse=" "),paste(loess.quartile,collapse=" "),paste(tagging.location,collapse=" "),
                                       paste(tagging.date,collapse=" "),paste(retrieval.date,collapse=" "),paste(sunrise.sd,collapse=" "),paste(sunset.sd,collapse=" "),paste(range.solar,collapse=" "),paste(speed.wet,collapse=" "),
                                       paste(speed.dry,collapse=" "),paste(sst.sd,collapse=" "),paste(max.sst.diff,collapse=" "),paste(ice.conc.cutoff,collapse=" "),paste(boundary.box,collapse=" "),paste(med.sea,collapse=" "),paste(black.sea,collapse=" "),
                                       paste(baltic.sea,collapse=" "),paste(caspian.sea,collapse=" "),paste(east.west.comp,collapse=" "),paste(wetdry.resolution,collapse=" "),paste(NOAA.OI.location,collapse=" "),
                                       paste(backward,collapse=" "),!is.null(sensor),!is.null(act),paste(land.mask,collapse=" "), distance.method, 'ParticleFilter'))


    # --- 1. Prepare Twilight Data (Proposals) ---
    message("Step 1: Preparing twilight data (for proposals)...")
    
    # browser()
    
    # Note: particle.number here defines the size of the proposal cloud *per step*
    # This might be different from the number of particles maintained in the filter.
    # Let's keep them the same for simplicity now (generate particle.number proposals per step).
    prep_twilight_dt <- prepare_twilight_data(
        trn = copy(trn),
        tagging.date = tagging.date, retrieval.date = retrieval.date,
        loess.quartile = loess.quartile, east.west.comp = east.west.comp,
        # Generate proposals equal to the number of particles in the filter
        particle.number = particle.number
    )

    # --- 2. Generate Initial Particle Cloud (Proposals) ---
    message("Step 2: Generating proposal particle cloud...")
    # This generates the base set of potential locations for *all* steps
    proposal_particles_dt <- generate_particle_cloud(
        trn_dt = prep_twilight_dt,
        sunrise.sd = sunrise.sd, sunset.sd = sunset.sd,
        range.solar = range.solar, tol = tol, boundary.box = boundary.box
    )
    all_generated_proposals <- copy(proposal_particles_dt) # Keep for output


    # --- 3. Filter Particle Clouds (Proposals) ---
    message("Step 3: Filtering proposal particle clouds...")
    filtered_proposals_dt <- filter_particle_clouds(
        trn_dt = proposal_particles_dt,
        # Filter based on the number generated per step
        particle.number = particle.number,
        boundary.box = boundary.box, land.mask = land.mask,
        med.sea = med.sea, black.sea = black.sea, baltic.sea = baltic.sea,
        caspian.sea = caspian.sea, NOAA.OI.location = NOAA.OI.location
    )

    if (nrow(filtered_proposals_dt) == 0) {
        stop("No proposal particles remaining after filtering steps.", call. = FALSE)
    }


    # --- 4. Prepare Tagging Location SF object ---
    # Minimal sf object needed for initialization
    tag.loc_df <- data.frame(lon = tagging.location[1], lat = tagging.location[2],
                            dtime = as.POSIXct(tagging.date, tz="UTC")) # Use tagging date for start
    tagging_sf <- st_as_sf(tag.loc_df, coords = c("lon", "lat"), crs = 4326, remove=FALSE)


    # --- 5. Run Particle Filter ---
    message("Step 4: Running Particle Filter...")
    # Use the filtered proposals as input
    summary_track_dt <- run_particle_filter(
        filtered_particles_dt = filtered_proposals_dt,
        particle.number = particle.number, # Number of particles in the filter
        tagging.location = tagging_sf,
        act = act, sensor = sensor,
        speed.dry = speed.dry, speed.wet = speed.wet,
        sst.sd = sst.sd, max.sst.diff = max.sst.diff,
        ice.conc.cutoff = ice.conc.cutoff, land.mask = land.mask,
        wetdry.resolution = wetdry.resolution,
        NOAA.OI.location = NOAA.OI.location,
        distance.method = distance.method,
        backward = backward
    )

    # --- 6. Final Output ---
    message("Step 5: Finalizing output...")

    # Convert summary track to sf object
    # Use median track as the primary output
    summary_track_sf <- st_as_sf(summary_track_dt, coords = c("lon_median", "lat_median"), crs = 4326, remove = FALSE)
    # Rename geometry column if desired
    # names(summary_track_sf)[names(summary_track_sf) == "geometry"] <- "median_location"

    # Convert all proposals to sf for output consistency
    all_proposals_sf <- st_as_sf(all_generated_proposals, coords=c("lon","lat"), crs=4326, remove=FALSE)

    end.time   <- Sys.time()
    time.taken <- difftime(end.time, start.time, units = "mins")

    message(paste('PF Algorithm run time:', round(as.numeric(time.taken), 1), 'min'))

    list.all <- list(
        'estimated track' = summary_track_sf, # Contains mean/median/sd
        # Optional: Could return the full final particle clouds if stored
        'all proposal particles' = all_proposals_sf, # All particles proposed at step 2
        'input parameters' = model.input,
        'model run time' = time.taken
    )

    options(warn = oldw)
    return(list.all)
}
