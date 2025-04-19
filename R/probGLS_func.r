#' Prepare Twilight Data for probGLS
#'
#' Cleans, filters, applies EWC, adds time variables, and duplicates
#' twilight data for particle generation.
#'
#' @param trn Raw twilight data.frame.
#' @param tagging.date Start date/time.
#' @param retrieval.date End date/time.
#' @param loess.quartile Loess filter parameter (or NULL).
#' @param east.west.comp Logical, apply East-West compensation?
#' @param particle.number Number of particles.
#' @return A data.table object ready for particle generation.
#' @import data.table
#' @importFrom GeoLight loessFilter coord
#' @importFrom stats median predict loess sd na.omit weighted.mean quantile
#' @importFrom lubridate as_datetime days seconds ymd_hms parse_date_time

prepare_twilight_data <- function(trn, tagging.date, retrieval.date,
                                  loess.quartile = NULL, east.west.comp = TRUE,
                                  particle.number) {
    
    # Ensure data.table is used
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' needed for this function to work.", call. = FALSE)
    }
    setDT(trn) # Convert to data.table by reference
    
    # Ensure datetime columns are POSIXct
    trn[, tFirst := as.POSIXct(tFirst, tz = "UTC")]
    trn[, tSecond := as.POSIXct(tSecond, tz = "UTC")]
    
    # Filter by date range
    trn <- trn[tFirst >= as.POSIXct(tagging.date, tz = "UTC") &
                             tSecond <= as.POSIXct(retrieval.date, tz = "UTC")]
    trn <- trn[!is.na(tFirst) & !is.na(tSecond)]
    
    if (nrow(trn) == 0) {
        stop('No twilight data points between selected tagging and retrieval dates.', call. = FALSE)
    }
    
    # Loess filter
    if (!is.null(loess.quartile)) {
        # loessFilter requires POSIXct - ensure they are
        trn_loess <- trn[, .(tFirst, tSecond, type)] # Select necessary columns
        loess_result <- GeoLight::loessFilter(trn_loess, plot = FALSE, k = loess.quartile)
        trn <- trn[loess_result == TRUE]
        if (nrow(trn) == 0) {
            stop('No twilight data points remaining after loess filter.', call. = FALSE)
        }
    }
    
    # East-West movement compensation
    if (east.west.comp) {
        # Ensure columns are numeric for calculations
        trn[, dtime_numeric := as.numeric(tFirst + as.numeric(difftime(tSecond,tFirst,units='sec'))/2)]
        
        # Calculate temporary longitude using a standard method (e.g., zenith=-6)
        # Note: coord requires POSIXct
        temp_coords <- GeoLight::coord(trn$tFirst, trn$tSecond, degElevation = -6, note = FALSE, method = 'NOAA')
        trn[, temp_lon := temp_coords[, 1]]
        
        # Shift requires numeric time; use the numeric version calculated earlier
        trn[, lon2 := -temp_lon]
        trn[, lon1 := shift(lon2, type = "lag")] # Use data.table's shift
        trn[, timedate2 := dtime_numeric]
        trn[, timedate1 := shift(timedate2, type = "lag")]
        trn[, twilight_length_sec := abs(as.numeric(difftime(tSecond, tFirst, units = "sec")))]
        
        # Calculate correction factor 'm' in degrees/hour (converted from sec)
        trn[, m := ((lon2 - lon1) / (timedate2 - timedate1)) * (twilight_length_sec / (15 * 2)) * 3600, by = seq_len(nrow(trn))] # Calculate row-wise if needed, careful with division by zero
        
        # Handle potential division by zero or NA from shift
        trn[is.na(m) | is.infinite(m), m := 0]
        
        # Apply correction (m is in degrees, convert to seconds: * 3600 / 15 ?)
        # The original formula `* (length / 15) / 2` seems complex. Let's assume 'm' directly represents the time shift in hours as intended.
        # Original code adds/subtracts m*3600 seconds.
        trn[, tFirst_corrected := tFirst + m * 3600]
        trn[, tSecond_corrected := tSecond - m * 3600]
        
        # Apply corrections, keep original if correction failed
        trn[!is.na(tFirst_corrected), tFirst := tFirst_corrected]
        trn[!is.na(tSecond_corrected), tSecond := tSecond_corrected]
        
        # Clean up temporary columns
        trn[, c("dtime_numeric", "temp_lon", "lon1", "lon2", "timedate1", "timedate2", "twilight_length_sec", "m", "tFirst_corrected", "tSecond_corrected") := NULL]
    }
    
    
    # Add essential time variables
    trn[, dtime := tFirst + as.numeric(difftime(tSecond, tFirst, units = 'sec')) / 2]
    trn[, doy := as.numeric(strftime(dtime, format = "%j", tz = "UTC"))]
    trn[, month := as.numeric(strftime(dtime, format = "%m", tz = "UTC"))]
    trn[, year := as.numeric(strftime(dtime, format = "%Y", tz = "UTC"))]
    trn[, jday := as.numeric(julian(dtime))] # Define julian for POSIXct
    
    # Keep only necessary columns before duplication
    trn <- trn[, .(tFirst, tSecond, type, dtime, doy, jday, year, month)]
    
    # Duplicate rows for particle number
    # Efficient duplication using data.table indexing
    trn <- trn[rep(seq_len(nrow(trn)), each = particle.number)]
    
    # Assign unique step identifier based on original twilight pair
    # Create grouping factor before duplication, then replicate it
    trn[, step_group := .GRP, by = .(floor(julian(tFirst)), floor(julian(tSecond)), type)] # Group before duplication
    trn[, step := rleid(step_group)] # Assign unique run-length ID to unique twilight events
    
    # Clean up group column
    trn[, step_group := NULL]
    
    
    return(trn)
}



#' Generate Initial Particle Cloud
#'
#' Adds random errors to twilights, calculates coordinates, handles equinox,
#' and filters by boundary box.
#'
#' @param trn_dt Prepared twilight data.table from `prepare_twilight_data`.
#' @param sunrise.sd Parameters for sunrise error distribution.
#' @param sunset.sd Parameters for sunset error distribution.
#' @param range.solar Min/max solar angle range.
#' @param tol Tolerance for SGAT::thresholdEstimate.
#' @param boundary.box Min/max lon/lat boundaries.
#' @return A data.table object containing the initial particle cloud.
#' @import data.table
#' @importFrom SGAT thresholdEstimate
#' @importFrom stats rlnorm runif

generate_particle_cloud <- function(trn_dt, sunrise.sd, sunset.sd, range.solar,
                                    tol, boundary.box) {
    
    # Add random solar angle
    solar.angle.steps <- seq(range.solar[1], range.solar[2], 0.01)
    trn_dt[, solar.angle := sample(solar.angle.steps, size = .N, replace = TRUE)]
    
    # Add twilight errors (vectorized using data.table's :=)
    # Calculate errors separately first for clarity if needed
    n_type1 <- nrow(trn_dt[type == 1])
    n_type2 <- nrow(trn_dt[type == 2])
    
    # Generate all errors at once
    sunrise_errors_type1 <- 60 * (rlnorm(n_type1, meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3])
    sunset_errors_type1    <- -60 * (rlnorm(n_type1, meanlog = sunset.sd[1], sdlog = sunset.sd[2]) + sunset.sd[3])
    sunset_errors_type2    <- -60 * (rlnorm(n_type2, meanlog = sunset.sd[1], sdlog = sunset.sd[2]) + sunset.sd[3])
    sunrise_errors_type2 <- 60 * (rlnorm(n_type2, meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3])
    
    # Assign errors based on type
    trn_dt[type == 1, tFirst.err := sunrise_errors_type1]
    trn_dt[type == 1, tSecond.err := sunset_errors_type1]
    trn_dt[type == 2, tFirst.err := sunset_errors_type2]
    trn_dt[type == 2, tSecond.err := sunrise_errors_type2]
    
    # Apply errors
    trn_dt[, tFirst := tFirst + tFirst.err]
    trn_dt[, tSecond := tSecond + tSecond.err]
    
    # Calculate coordinates
    # Ensure inputs to thresholdEstimate are POSIXct
    pos <- SGAT::thresholdEstimate(
        trise = ifelse(trn_dt$type == 1, trn_dt$tFirst, trn_dt$tSecond),
        tset = ifelse(trn_dt$type == 1, trn_dt$tSecond, trn_dt$tFirst),
        zenith = 90 - trn_dt$solar.angle,
        tol = tol
    )
    trn_dt[, lon := pos[, 1]]
    trn_dt[, lat := pos[, 2]]
    
    # Handle Equinox periods / NA latitudes (Original logic: random latitudes)
    # Note: This assumes NA latitudes primarily occur around equinoxes.
    # Consider if a date-based equinox definition is more robust if needed.
    na_lat_indices <- which(is.na(trn_dt$lat))
    if (length(na_lat_indices) > 0) {
        trn_dt[na_lat_indices, lat := runif(.N, min = boundary.box[3], max = boundary.box[4])]
        # Also nullify solar angle influence where lat was NA originally
        trn_dt[na_lat_indices, solar.angle := NA]
    }
    
    
    # Filter by boundary box
    min_lon <- boundary.box[1]
    max_lon <- boundary.box[2]
    min_lat <- boundary.box[3]
    max_lat <- boundary.box[4]
    
    if (min_lon > max_lon) { # Crosses dateline
        trn_dt <- trn_dt[(lon > min_lon | lon < max_lon) & lat > min_lat & lat < max_lat]
    } else {
        trn_dt <- trn_dt[lon > min_lon & lon < max_lon & lat > min_lat & lat < max_lat]
    }
    
    return(trn_dt)
}


#' Filter Particle Clouds
#'
#' Filters particle clouds based on minimum particle count per step,
#' latitude range checks, and applies land mask weighting.
#'
#' @param trn_dt Particle cloud data.table from `generate_particle_cloud`.
#' @param particle.number Original number of particles per step.
#' @param boundary.box Min/max lon/lat boundaries.
#' @param land.mask Logical or NULL, apply land mask? T=ocean, F=land.
#' @param med.sea Logical, mask Mediterranean Sea?
#' @param black.sea Logical, mask Black Sea?
#' @param baltic.sea Logical, mask Baltic Sea?
#' @param caspian.sea Logical, mask Caspian Sea?
#' @param NOAA.OI.location Path to land mask file.
#' @return Filtered data.table with land mask weights.
#' @import data.table
#' @importFrom utils download.file

filter_particle_clouds <- function(trn_dt, particle.number, boundary.box,
                                    land.mask, med.sea, black.sea,
                                    baltic.sea, caspian.sea,
                                    NOAA.OI.location) {
    
    min_particles_per_step <- particle.number / 5
    min_particles_land         <- particle.number * 0.1
    
    # Remove steps with too few particles overall
    trn_dt[, n_in_step := .N, by = step]
    trn_dt <- trn_dt[n_in_step >= min_particles_per_step]
    
    if (nrow(trn_dt) == 0) {
        stop("No particle clouds remaining after minimum count filter.", call. = FALSE)
    }
    
    # Remove clouds completely outside latitude bounds
    trn_dt[, max_lat_step := max(lat), by = step]
    trn_dt[, min_lat_step := min(lat), by = step]
    trn_dt <- trn_dt[max_lat_step >= boundary.box[3] & min_lat_step <= boundary.box[4]]
    
    if (nrow(trn_dt) == 0) {
        stop("No particle clouds remaining after latitude range filter.", call. = FALSE)
    }
    
    # Apply land mask if requested
    if (!is.null(land.mask)) {
        landmask_file <- list.files(path = NOAA.OI.location, pattern = "lsmask.oisst.nc", recursive = TRUE, full.names = TRUE)
        if (length(landmask_file) == 0) {
            cat('\nNo land mask file found - attempting download...\n')
            tryCatch({
                download.file('https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.nc',
                                            destfile = file.path(NOAA.OI.location, 'lsmask.oisst.nc'), mode = "wb")
                landmask_file <- file.path(NOAA.OI.location, 'lsmask.oisst.nc')
                if (!file.exists(landmask_file)) stop("Download failed.")
                cat('Land mask downloaded.\n')
            }, error = function(e) {
                stop(paste('Land mask file not found in', NOAA.OI.location, 'and download failed. Please download lsmask.oisst.nc manually.'), call. = FALSE)
            })
            
        } else {
            landmask_file <- landmask_file[1] # Take the first one if multiple found
        }
        
        # Assuming load_landmask returns 0 for land, 1 for water
        trn_dt[, landmask_val := load_landmask(FILE_NAME = landmask_file, LONS = lon, LATS = lat)]
        
        # Apply specific sea masks (modify landmask_val)
        if (baltic.sea)    trn_dt[lon > 14 & lon < 33.5 & lat > 51.4 & lat < 66.2, landmask_val := 0]
        if (med.sea) {
            trn_dt[(lon >= 0 & lon <= 27 & lat > 30 & lat < 48) |
                             (lon >= 27 & lon < 40 & lat > 30 & lat < 40) |
                             (lon > 355 & lon <= 360 & lat > 30 & lat < 42), landmask_val := 0]
        }
        if (black.sea)     trn_dt[lon > 27 & lon < 45 & lat > 40 & lat < 48, landmask_val := 0]
        if (caspian.sea) trn_dt[lon > 45 & lon < 62 & lat > 35 & lat < 48, landmask_val := 0]
        
        # Assign weight based on desired habitat (ocean or land)
        trn_dt[, weight_land := 0]
        if (land.mask == TRUE) { # Ocean habitat
            trn_dt[landmask_val == 1, weight_land := 1]
            # Filter steps with too few particles in water
            trn_dt[, n_water := sum(weight_land), by = step]
            trn_dt <- trn_dt[n_water >= min_particles_land]
        } else { # Land habitat
            trn_dt[landmask_val == 0, weight_land := 1]
            # Filter steps with too few particles on land
            trn_dt[, n_land := sum(weight_land), by = step]
            trn_dt <- trn_dt[n_land >= min_particles_land]
        }
        
        if (nrow(trn_dt) == 0) {
            stop("No particle clouds remaining after land mask filtering.", call. = FALSE)
        }
        trn_dt[, c("landmask_val", "n_water", "n_land") := NULL] # Clean up
        
    } else {
        # If no land mask, assign neutral weight
        trn_dt[, weight_land := 1]
    }
    
    # Clean up intermediate columns
    trn_dt[, c("n_in_step", "max_lat_step", "min_lat_step") := NULL]
    
    # Add final date/time columns needed for iteration
    trn_dt[, dtime := tFirst + as.numeric(difftime(tSecond, tFirst, units = 'sec')) / 2]
    trn_dt[, jday_int := floor(jday)] # Julian day (integer)
    trn_dt[, date := as.Date(dtime, tz = "UTC")] # Date object
    
    # Ensure order for stepping
    setorder(trn_dt, dtime)
    
    return(trn_dt)
}
