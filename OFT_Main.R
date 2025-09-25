source(r"(H:\OFT\DLCAnalyzer_Functions_final_v7.R)")
source(r"(H:\OFT\OFT_Analysis_Pipeline.R)")

setwd(r"(H:\OFT\honli)")

config <- list(
  # Directories
  working_dir = getwd(),
  raw_data_dir = "filtered_csv",
  annotation_dir = "Annotation",
  output_dir = "OFT_Results",
  
  # Video parameters
  fps = 0,                    # Frames per second (will be auto-calculated if set to 0)
  video_length = 5,           # Video length in minutes
  
  # Arena parameters
  arena_size = 50,            # Arena size in cm (for square arena)
  
  # OFT zone scales
  center_scale = 0.6,         # Scale factor for center zone (0.6 = 9/25 of area)
  periphery_scale = 0.8,      # Scale factor for periphery zone
  corner_scale = 0.25,        # Scale factor for corner zones
  
  # Movement parameters
  movement_cutoff = 3,        # cm/s threshold for movement detection
  integration_period = 3,     # Frames for smoothing
  likelihood_cutoff = 0.9,    # DLC confidence threshold
  
  # Analysis options
  interpolate_jumps = FALSE,  # Clean tracking jumps (WIP)
  save_plots = TRUE,          # Save trajectory plots
  plot_dpi = 300,             # Plot resolution
  analyze_temporal = TRUE,    # Perform per-minute analysis
  skip_merging = TRUE,        # Skip merging step if TRUE (use existing tracking files)
  correct_arena_tilt = TRUE,  # Automatically correct arena tilt
  color_path_by_speed = FALSE,
  debug_mode = FALSE,          # Additional debug outputs
  
  # Visualization options
  heatmap_palette = "viridis::turbo",
  plot_format = "png"         # Options: "png", "svg", "eps",
)

# Run this code
results <- run_analysis()
