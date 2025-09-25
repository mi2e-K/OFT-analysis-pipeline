library(sp)
library(imputeTS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(cowplot)
library(stringr)
library(paletteer)
library(tools)
library(svglite)  # Add for SVG support

################################################################################
#### STEP 1: DATA MERGING FUNCTIONS ####
################################################################################

merge_tracking_and_annotation <- function(tracking_file, annotation_file, output_file = NULL, output_dir = "tracking_csv") {
  # Merge DeepLabCut tracking data with arena corner annotations
  # Creates a single file compatible with DLCAnalyzer
  
  cat("Merging tracking and annotation data...\n")
  
  # Create tracking_csv directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  Created directory:", output_dir, "\n")
  }
  
  # Read tracking data (mouse body parts) - keep as character to preserve structure
  tracking_lines <- readLines(tracking_file)
  
  # Read annotation data (arena corners)
  annotation <- read_csv(annotation_file, col_names = FALSE, show_col_types = FALSE)
  
  # Process annotation to extract corner information
  corner_names <- as.character(annotation[2, -(1:3)])
  corner_names <- corner_names[!is.na(corner_names) & corner_names != ""]
  unique_corners <- unique(corner_names)
  
  cat("  Found corners:", paste(unique_corners, collapse = ", "), "\n")
  
  # Get the coordinate values from row 4
  corner_values <- as.numeric(annotation[4, -(1:3)])
  
  # Parse the tracking file header
  header1 <- strsplit(tracking_lines[1], ",")[[1]]  # scorer
  header2 <- strsplit(tracking_lines[2], ",")[[1]]  # bodyparts
  header3 <- strsplit(tracking_lines[3], ",")[[1]]  # coords
  
  # Find where body parts end
  body_parts <- c("nose", "bodycenter", "bodycentre", "tailbase", "head", "neck", 
                  "forepaw", "hindpaw", "tail", "ear", "paw", "body")
  
  last_body_idx <- 1
  for (i in seq_along(header2)) {
    if (header2[i] != "" && !is.na(header2[i])) {
      if (any(grepl(paste(body_parts, collapse = "|"), header2[i], ignore.case = TRUE))) {
        if (i + 2 <= length(header2) && 
            header2[i] == header2[i+1] && 
            header2[i] == header2[i+2]) {
          last_body_idx <- i + 2
        }
      }
    }
  }
  
  cat("  Keeping columns 1 to", last_body_idx, "from original tracking\n")
  
  # Keep only body part columns
  header1_body <- header1[1:last_body_idx]
  header2_body <- header2[1:last_body_idx]
  header3_body <- header3[1:last_body_idx]
  
  # Build corner headers
  header1_corners <- c()
  header2_corners <- c()
  header3_corners <- c()
  
  for (corner in unique_corners) {
    header1_corners <- c(header1_corners, rep(header1[2], 3))
    header2_corners <- c(header2_corners, corner, corner, corner)
    header3_corners <- c(header3_corners, "x", "y", "likelihood")
  }
  
  # Combine headers
  new_header1 <- c(header1_body, header1_corners)
  new_header2 <- c(header2_body, header2_corners)
  new_header3 <- c(header3_body, header3_corners)
  
  # Process data rows
  n_data_rows <- length(tracking_lines) - 3
  corner_data_matrix <- matrix(nrow = n_data_rows, ncol = length(unique_corners) * 3)
  
  col_idx <- 1
  value_idx <- 1
  for (corner in unique_corners) {
    x_val <- corner_values[value_idx]
    y_val <- corner_values[value_idx + 1]
    
    corner_data_matrix[, col_idx] <- x_val
    corner_data_matrix[, col_idx + 1] <- y_val
    corner_data_matrix[, col_idx + 2] <- 1  # likelihood = 1 for corners
    
    col_idx <- col_idx + 3
    value_idx <- value_idx + 2
  }
  
  # Combine tracking data with corner data
  new_lines <- c()
  new_lines[1] <- paste(new_header1, collapse = ",")
  new_lines[2] <- paste(new_header2, collapse = ",")
  new_lines[3] <- paste(new_header3, collapse = ",")
  
  for (i in 4:length(tracking_lines)) {
    row_data <- strsplit(tracking_lines[i], ",")[[1]]
    body_data <- row_data[1:last_body_idx]
    corner_data <- corner_data_matrix[i - 3, ]
    combined_row <- c(body_data, corner_data)
    new_lines[i] <- paste(combined_row, collapse = ",")
  }
  
  # Save merged file in tracking_csv directory
  if (is.null(output_file)) {
    output_file <- file.path(output_dir, gsub(".*DLC", "", basename(tracking_file)))
    output_file <- gsub("_filtered\\.csv", "_tracking.csv", output_file)
  } else {
    output_file <- file.path(output_dir, output_file)
  }
  
  writeLines(new_lines, output_file)
  
  cat("  Merged file saved:", basename(output_file), "\n")
  
  # Verify the merge
  verify_merged_file(output_file)
  
  return(output_file)
}

fix_and_merge_files <- function(tracking_file, annotation_file, output_file = NULL, output_dir = "tracking_csv") {
  # Alternative merge function that handles empty columns properly
  
  cat("Processing files with empty column handling...\n")
  
  # Create tracking_csv directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  Created directory:", output_dir, "\n")
  }
  
  # Read the raw files
  tracking_lines <- readLines(tracking_file)
  annotation <- read_csv(annotation_file, col_names = FALSE, show_col_types = FALSE)
  
  # Parse headers
  header1_parts <- strsplit(tracking_lines[1], ",")[[1]]
  header2_parts <- strsplit(tracking_lines[2], ",")[[1]]
  header3_parts <- strsplit(tracking_lines[3], ",")[[1]]
  
  # Find the actual end of body part data
  actual_body_end <- 1
  i <- 2
  while (i <= length(header2_parts)) {
    if (!is.na(header2_parts[i]) && header2_parts[i] != "") {
      if (i + 2 <= length(header2_parts) &&
          header2_parts[i] == header2_parts[i+1] &&
          header2_parts[i] == header2_parts[i+2]) {
        actual_body_end <- i + 2
        i <- i + 3
      } else {
        break
      }
    } else {
      break
    }
  }
  
  cat("  Body part columns: 1 to", actual_body_end, "\n")
  
  # Extract only valid body part headers
  header1_clean <- header1_parts[1:actual_body_end]
  header2_clean <- header2_parts[1:actual_body_end]
  header3_clean <- header3_parts[1:actual_body_end]
  
  # Get corner data from annotation
  corner_names <- as.character(annotation[2, -(1:3)])
  corner_names <- corner_names[!is.na(corner_names) & corner_names != ""]
  unique_corners <- unique(corner_names)
  
  cat("  Corners to add:", paste(unique_corners, collapse = ", "), "\n")
  
  # Get corner coordinates
  corner_values <- as.numeric(annotation[4, -(1:3)])
  
  # Build new headers with corners
  for (corner in unique_corners) {
    header1_clean <- c(header1_clean, rep(header1_parts[2], 3))
    header2_clean <- c(header2_clean, corner, corner, corner)
    header3_clean <- c(header3_clean, "x", "y", "likelihood")
  }
  
  # Create new file content
  new_lines <- character()
  new_lines[1] <- paste(header1_clean, collapse = ",")
  new_lines[2] <- paste(header2_clean, collapse = ",")
  new_lines[3] <- paste(header3_clean, collapse = ",")
  
  n_frames <- length(tracking_lines) - 3
  
  for (i in 4:length(tracking_lines)) {
    row_parts <- strsplit(tracking_lines[i], ",")[[1]]
    body_data <- row_parts[1:actual_body_end]
    
    corner_data <- c()
    coord_idx <- 1
    for (corner in unique_corners) {
      corner_data <- c(corner_data, 
                       corner_values[coord_idx],
                       corner_values[coord_idx + 1],
                       1)  # likelihood
      coord_idx <- coord_idx + 2
    }
    
    new_lines[i] <- paste(c(body_data, corner_data), collapse = ",")
  }
  
  # Save the cleaned and merged file in tracking_csv directory
  if (is.null(output_file)) {
    output_file <- file.path(output_dir, gsub(".*DLC", "", basename(tracking_file)))
    output_file <- gsub("_filtered\\.csv", "_tracking.csv", output_file)
  } else {
    output_file <- file.path(output_dir, output_file)
  }
  
  writeLines(new_lines, output_file)
  cat("  Saved clean merged file:", basename(output_file), "\n")
  
  # Verify
  verify_merged_file(output_file)
  
  return(output_file)
}

################################################################################
#### STEP 2: ARENA TILT CORRECTION FUNCTION ####
################################################################################

correct_arena_tilt <- function(t, corner_points = NULL) {
  # Automatically detect and correct arena tilt
  
  cat("  Correcting arena tilt...\n")
  
  # If corner points not specified, try to find them
  if (is.null(corner_points)) {
    corner_patterns <- c("top_Left", "top_Right", "bottom_Right", "bottom_Left",
                         "UpperLeft", "UpperRight", "LowerRight", "LowerLeft",
                         "UL", "UR", "LR", "LL",
                         "topleft", "topright", "bottomright", "bottomleft")
    
    available_points <- names(t$data)
    
    # Try to find a complete set of corners
    for (i in seq(1, length(corner_patterns), by = 4)) {
      test_corners <- corner_patterns[i:(i+3)]
      if (all(test_corners %in% available_points)) {
        corner_points <- test_corners
        break
      }
    }
  }
  
  if (is.null(corner_points) || length(corner_points) < 4) {
    cat("    Warning: Could not find corner points for tilt correction\n")
    return(t)
  }
  
  # Get corner coordinates
  corners <- t$median.data[corner_points[1:4], c("x", "y")]
  
  # Calculate the angle of the top edge (should be horizontal)
  # Use first two corners (assumed to be top-left and top-right)
  dx <- corners$x[2] - corners$x[1]
  dy <- corners$y[2] - corners$y[1]
  
  # Calculate angle in radians
  angle_rad <- atan2(dy, dx)
  angle_deg <- angle_rad * 180 / pi
  
  cat("    Detected tilt angle:", round(angle_deg, 2), "degrees\n")
  
  # Apply rotation if tilt is significant (> 1 degree)
  if (abs(angle_deg) > 0.3) {
    # Rotate the tracking data to correct the tilt
    t <- RotateTrackingData(t, theta = -angle_deg, center.of = corner_points)
    cat("    Applied rotation correction\n")
  } else {
    cat("    Arena is already aligned (tilt < 1 degree)\n")
  }
  
  return(t)
}

################################################################################
#### STEP 3: MODIFIED OFT ANALYSIS FUNCTION ####
################################################################################

analyze_single_OFT <- function(tracking_file, config, output_dir = NULL) {
  # Complete OFT analysis for a single file
  
  # Extract file ID
  file_id <- gsub("_tracking.*", "", basename(tracking_file))
  file_id <- gsub("DLC.*", "", file_id)
  
  cat("\nAnalyzing:", file_id, "\n")
  
  # Create output directories
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(tracking_file), "..", config$output_dir)
  }
  
  dirs <- list(
    main = output_dir,
    results = file.path(output_dir, "Results"),
    plots = file.path(output_dir, "Plots"),
    temporal = file.path(output_dir, "Temporal")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  # Read tracking data
  t <- ReadDLCDataFromCSV(tracking_file, fps = 1)
  
  # Auto-calculate FPS if needed
  if (config$fps == 0) {
    total_frames <- length(t$frames)
    t$fps <- round(total_frames / (config$video_length * 60))
    cat("  Auto-calculated FPS:", t$fps, "\n")
  } else {
    t$fps <- config$fps
  }
  
  # Update seconds based on FPS
  t$seconds <- t$frames / t$fps
  
  # Clean tracking data if requested
  if (config$interpolate_jumps) {
    t_clean <- clean_tracking_jumps(tracking_file)
    t$data <- t_clean$data
  }
  
  # Find corner points for calibration
  corner_patterns <- c("top_Left", "top_Right", "bottom_Right", "bottom_Left",
                       "UpperLeft", "UpperRight", "LowerRight", "LowerLeft",
                       "UL", "UR", "LR", "LL",
                       "topleft", "topright", "bottomright", "bottomleft")
  
  available_points <- names(t$data)
  corner_points <- NULL
  
  # Try to find a complete set of corners
  for (i in seq(1, length(corner_patterns), by = 4)) {
    test_corners <- corner_patterns[i:(i+3)]
    if (all(test_corners %in% available_points)) {
      corner_points <- test_corners
      break
    }
  }
  
  if (is.null(corner_points)) {
    cat("  Warning: Standard corner names not found.\n")
    cat("  Available points:", paste(available_points, collapse = ", "), "\n")
    
    # Try to find corners by pattern matching
    corner_candidates <- available_points[grepl("Left|Right|top|bottom|Upper|Lower|corner", 
                                                available_points)]
    
    if (length(corner_candidates) >= 4) {
      corner_points <- corner_candidates[1:4]
      cat("  Using corners:", paste(corner_points, collapse = ", "), "\n")
    } else {
      stop("Cannot find arena corners for calibration!")
    }
  }
  
  # Apply arena tilt correction if enabled
  if (config$correct_arena_tilt) {
    t <- correct_arena_tilt(t, corner_points)
  }
  
  # Calibrate tracking data
  arena_area <- config$arena_size * config$arena_size
  t <- CalibrateTrackingData(t, 
                             method = "area",
                             in.metric = arena_area,
                             points = corner_points)
  
  cat("  Calibration complete (px to cm ratio:", round(t$px.to.cm, 4), ")\n")
  
  # Add OFT zones
  t <- AddOFTZones(t, 
                   points = corner_points,
                   scale_center = config$center_scale,
                   scale_periphery = config$periphery_scale,
                   scale_corners = config$corner_scale)
  
  # Clean tracking data
  t <- CleanTrackingData(t, 
                         likelihoodcutoff = config$likelihood_cutoff,
                         maxdelta = config$arena_size * 0.5)
  
  # Calculate movement
  t <- CalculateMovement(t, 
                         movement_cutoff = config$movement_cutoff,
                         integration_period = config$integration_period)
  
  # Find tracking point
  tracking_point <- "bodycenter"
  if (!tracking_point %in% names(t$data)) {
    if ("bodycentre" %in% names(t$data)) {
      tracking_point <- "bodycentre"
    } else if ("body" %in% names(t$data)) {
      tracking_point <- "body"
    } else {
      non_corner_points <- setdiff(names(t$data), corner_points)
      if (length(non_corner_points) > 0) {
        tracking_point <- non_corner_points[1]
        cat("  Using tracking point:", tracking_point, "\n")
      }
    }
  }
  
  # Run OFT analysis
  t <- OFTAnalysis(t,
                   movement_cutoff = config$movement_cutoff,
                   integration_period = config$integration_period,
                   points = tracking_point)
  
  # Extract results
  results <- data.frame(
    ID = file_id,
    Duration_min = length(t$frames) / t$fps / 60,
    
    # General locomotion
    Total_Distance_cm = t$Report[[paste0(tracking_point, ".raw.distance")]],
    Distance_Moving_cm = t$Report[[paste0(tracking_point, ".distance.moving")]],
    Average_Speed_cm_s = t$Report[[paste0(tracking_point, ".raw.speed")]],
    Speed_Moving_cm_s = t$Report[[paste0(tracking_point, ".speed.moving")]],
    
    # Time metrics
    Time_Moving_s = t$Report[[paste0(tracking_point, ".time.moving")]],
    Time_Stationary_s = t$Report[[paste0(tracking_point, ".time.stationary")]],
    Percent_Moving = t$Report[[paste0(tracking_point, ".percentage.moving")]],
    
    # Center zone
    Center_Time_s = t$Report[[paste0(tracking_point, ".center.total.time")]],
    Center_Distance_cm = t$Report[[paste0(tracking_point, ".center.raw.distance")]],
    Center_Transitions = t$Report[[paste0(tracking_point, ".center.transitions")]] / 2,
    Center_Percent = (t$Report[[paste0(tracking_point, ".center.total.time")]] / 
                        t$Report[[paste0(tracking_point, ".total.time")]]) * 100,
    
    # Periphery zone
    Periphery_Time_s = t$Report[[paste0(tracking_point, ".periphery.total.time")]],
    Periphery_Distance_cm = t$Report[[paste0(tracking_point, ".periphery.raw.distance")]],
    Periphery_Transitions = t$Report[[paste0(tracking_point, ".periphery.transitions")]] / 2,
    Periphery_Percent = (t$Report[[paste0(tracking_point, ".periphery.total.time")]] / 
                           t$Report[[paste0(tracking_point, ".total.time")]]) * 100,
    
    # Corners combined
    Corners_Time_s = t$Report[[paste0(tracking_point, ".corners.total.time")]],
    Corners_Distance_cm = t$Report[[paste0(tracking_point, ".corners.raw.distance")]],
    Corners_Transitions = t$Report[[paste0(tracking_point, ".corners.transitions")]] / 2,
    Corners_Percent = (t$Report[[paste0(tracking_point, ".corners.total.time")]] / 
                         t$Report[[paste0(tracking_point, ".total.time")]]) * 100,
    
    # Anxiety metrics
    Center_Periphery_Ratio = t$Report[[paste0(tracking_point, ".center.total.time")]] / 
      (t$Report[[paste0(tracking_point, ".periphery.total.time")]] + 0.001),
    Thigmotaxis_Index = t$Report[[paste0(tracking_point, ".periphery.total.time")]] / 
      (t$Report[[paste0(tracking_point, ".center.total.time")]] + 
         t$Report[[paste0(tracking_point, ".periphery.total.time")]])
  )
  
  # Save results
  write.csv(results, 
            file.path(dirs$results, paste0(file_id, "_results.csv")),
            row.names = FALSE)
  
  cat("  Results saved\n")
  
  # Temporal analysis
  if (config$analyze_temporal) {
    cat("  Running temporal analysis...\n")
    
    temporal_data <- data.frame()
    segment_length <- 60  # 1 minute segments
    total_seconds <- length(t$frames) / t$fps
    n_segments <- ceiling(total_seconds / segment_length)
    
    for (seg in 1:n_segments) {
      start_frame <- (seg - 1) * segment_length * t$fps + 1
      end_frame <- min(seg * segment_length * t$fps, length(t$frames))
      
      # Check if segment has enough frames for analysis (need at least 2 for movement calculation)
      if (end_frame - start_frame < 1) {
        cat("    Skipping segment", seg, "(insufficient frames)\n")
        next
      }
      
      t_seg <- CutTrackingData(t, keep.frames = start_frame:end_frame)
      
      # Only calculate movement if we have enough frames
      if (length(t_seg$frames) >= 2) {
        t_seg <- CalculateMovement(t_seg, 
                                   config$movement_cutoff, 
                                   config$integration_period)
        
        # Get zone times for segment
        seg_center <- sum(IsInZone(t_seg, tracking_point, "center")) / t$fps
        seg_periphery <- sum(IsInZone(t_seg, tracking_point, "periphery", invert = TRUE)) / t$fps
        seg_corners <- sum(IsInZone(t_seg, tracking_point, t$corner.names)) / t$fps
        
        # Calculate distance in segment
        seg_distance <- sum(t_seg$data[[tracking_point]]$speed, na.rm = TRUE)
      } else {
        # For segments with only 1 frame, set movement-based metrics to 0
        seg_center <- sum(IsInZone(t_seg, tracking_point, "center")) / t$fps
        seg_periphery <- sum(IsInZone(t_seg, tracking_point, "periphery", invert = TRUE)) / t$fps
        seg_corners <- sum(IsInZone(t_seg, tracking_point, t$corner.names)) / t$fps
        seg_distance <- 0
      }
      
      temporal_row <- data.frame(
        Minute = seg,
        Distance_cm = seg_distance,
        Center_Time_s = seg_center,
        Periphery_Time_s = seg_periphery,
        Corners_Time_s = seg_corners,
        Center_Percent = (seg_center / min(segment_length, length(t_seg$frames) / t$fps)) * 100
      )
      
      temporal_data <- rbind(temporal_data, temporal_row)
    }
    
    write.csv(temporal_data,
              file.path(dirs$temporal, paste0(file_id, "_temporal.csv")),
              row.names = FALSE)
    
    cat("  Temporal analysis saved\n")
  }
  
  # Generate plots
  if (config$save_plots) {
    cat("  Generating plots...\n")
    
    # Trajectory plot with custom color palette
    p <- PlotDensityPaths(t, points = tracking_point, Title = file_id, 
                          color_palette = config$heatmap_palette,
                          color_by_speed = ifelse(is.null(config$color_path_by_speed), FALSE, config$color_path_by_speed))
    
    # Add zone boundaries (center zone)
    if (!is.null(t$zones$center)) {
      center_zone <- t$zones$center
      center_zone_closed <- rbind(center_zone, center_zone[1,])
      
      p[[tracking_point]] <- p[[tracking_point]] + 
        geom_path(data = center_zone_closed,
                  aes(x = x, y = y),
                  size = 1,
                  color = "#ff4b00",
                  linetype = "longdash")
    }
    
    # Save plot in specified format
    if (config$plot_format == "svg") {
      ggsave(file.path(dirs$plots, paste0(file_id, "_trajectory.svg")),
             plot = p[[tracking_point]],
             width = 8, height = 8,
             device = svglite)
    } else {
      ggsave(file.path(dirs$plots, paste0(file_id, "_trajectory.", config$plot_format)),
             plot = p[[tracking_point]],
             width = 8, height = 8,
             dpi = config$plot_dpi)
    }
    
    # Zone timeline (if debug mode)
    if (config$debug_mode) {
      p_timeline <- PlotZoneVisits(t, points = tracking_point)
      
      if (config$plot_format == "svg") {
        ggsave(file.path(dirs$plots, paste0(file_id, "_timeline.svg")),
               plot = p_timeline,
               width = 12, height = 4,
               device = svglite)
      } else {
        ggsave(file.path(dirs$plots, paste0(file_id, "_timeline.", config$plot_format)),
               plot = p_timeline,
               width = 12, height = 4,
               dpi = config$plot_dpi)
      }
    }
    
    cat("  Plots saved\n")
  }
  
  cat("  Analysis complete for", file_id, "\n")
  
  return(list(
    results = results,
    temporal = if(config$analyze_temporal) temporal_data else NULL,
    tracking = t
  ))
}

################################################################################
#### STEP 4: MODIFIED BATCH PROCESSING ####
################################################################################

process_OFT_batch <- function(config) {
  # Complete pipeline: merge data, clean, and analyze all files
  
  cat("\n========================================\n")
  cat("     OFT BATCH ANALYSIS PIPELINE\n")
  cat("========================================\n\n")
  
  # Set working directory
  setwd(config$working_dir)
  
  # Create main output directory
  output_dir <- file.path(config$working_dir, config$output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create tracking_csv directory
  tracking_csv_dir <- file.path(config$working_dir, "tracking_csv")
  if (!dir.exists(tracking_csv_dir)) {
    dir.create(tracking_csv_dir, recursive = TRUE)
  }
  
  all_results <- data.frame()
  processed_files <- c()
  
  # Check if merging should be skipped
  if (config$skip_merging) {
    cat("Step 1: Skipping merging step (using existing tracking files)...\n")
    
    # Look for tracking files in tracking_csv directory
    tracking_files <- list.files(tracking_csv_dir, 
                                 pattern = "_tracking\\.csv$", 
                                 full.names = TRUE)
    
    cat("  Found", length(tracking_files), "tracking files in tracking_csv\n")
    
  } else {
    # Step 1: Find and merge files
    cat("Step 1: Finding and merging data files...\n")
    
    tracking_files <- list.files(config$raw_data_dir, 
                                 pattern = "_filtered\\.csv$", 
                                 full.names = TRUE)
    
    if (length(tracking_files) == 0) {
      tracking_files <- list.files(config$working_dir,
                                   pattern = "_filtered\\.csv$",
                                   full.names = TRUE)
    }
    
    cat("  Found", length(tracking_files), "raw tracking files\n")
    
    # Process merging
    merged_files <- c()
    
    for (tracking_file in tracking_files) {
      file_id <- basename(tracking_file)
      file_id <- gsub("DLC.*", "", file_id)
      file_id <- trimws(file_id)
      
      # Find annotation file
      anno_patterns <- c(
        file.path(config$annotation_dir, file_id, paste0(file_id, "_annotation.csv")),
        file.path(config$annotation_dir, file_id, paste0(file_id, ".csv")),
        file.path(config$annotation_dir, paste0(file_id, "_annotation.csv")),
        file.path(config$working_dir, paste0(file_id, "_annotation.csv"))
      )
      
      annotation_file <- NULL
      for (pattern in anno_patterns) {
        if (file.exists(pattern)) {
          annotation_file <- pattern
          break
        }
      }
      
      if (is.null(annotation_file)) {
        cat("  Warning: No annotation file found for", file_id, ", skipping merge\n")
        next
      }
      
      # Merge and save to tracking_csv directory
      cat("  Merging", file_id, "...\n")
      merged_file <- fix_and_merge_files(
        tracking_file, 
        annotation_file,
        output_dir = tracking_csv_dir
      )
      merged_files <- c(merged_files, merged_file)
    }
    
    tracking_files <- merged_files
  }
  
  # Step 2: Process each tracking file
  for (i in 1:length(tracking_files)) {
    
    tryCatch({
      
      tracking_file <- tracking_files[i]
      
      # Extract ID
      file_id <- basename(tracking_file)
      file_id <- gsub("_tracking.*", "", file_id)
      file_id <- trimws(file_id)
      
      cat("\n[", i, "/", length(tracking_files), "] Processing:", file_id, "\n")
      
      # Run OFT analysis
      result <- analyze_single_OFT(tracking_file, config, output_dir)
      
      # Collect results
      all_results <- rbind(all_results, result$results)
      processed_files <- c(processed_files, file_id)
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
    })
  }
  
  # Step 3: Save combined results
  if (nrow(all_results) > 0) {
    cat("\n========================================\n")
    cat("Saving combined results...\n")
    
    # Save all results
    write.csv(all_results,
              file.path(output_dir, "Combined_Results.csv"),
              row.names = FALSE)
    
    cat("  Combined results saved\n")
    cat("  Successfully processed", length(processed_files), "files\n")
  }
  
  cat("\n========================================\n")
  cat("Analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  cat("========================================\n\n")
  
  return(list(
    results = all_results,
    files = processed_files
  ))
}

################################################################################
#### MAIN EXECUTION ####
################################################################################

# Example usage:
run_analysis <- function() {
  # Option 1: Run full batch analysis with default config
  results <- process_OFT_batch(config)
  
  return(results)
}

# Quick single file analysis
analyze_single_file <- function(tracking_file, annotation_file = NULL) {
  
  # If annotation provided, merge first
  if (!is.null(annotation_file)) {
    tracking_file <- fix_and_merge_files(tracking_file, annotation_file, 
                                         output_dir = "tracking_csv")
  }
  
  # Analyze
  result <- analyze_single_OFT(tracking_file, config)
  
  return(result)
}