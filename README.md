# OFT Analysis Pipeline for DeepLabCut Data

An R-based pipeline for analyzing Open Field Test (OFT) behavioral data from DeepLabCut tracking output.

### Required Files
- **DLCAnalyzer_Functions_final_v7.R** - Core DLCAnalyzer functions
- **OFT_Analysis_Pipeline.R** - Main analysis functions
- **OFT_Main.R** - Configuration and execution script

```

3. Place `DLCAnalyzer_Functions_final_v7.R` in your R code directory

4. Update file paths in `OFT_Main.R`:
```r
source("path/to/DLCAnalyzer_Functions_final_v7.R")
source("path/to/OFT_Analysis_Pipeline.R")
setwd("path/to/your/data")
```

## Data Structure

### Input Directory Structure
```
Working Directory/
├── filtered_csv/           # DLC output files (*_filtered.csv)
└── Annotation/             # Arena corner annotations
    ├── Subject1/
    │   └── Subject1_annotation.csv
    └── Subject2/
        └── Subject2_annotation.csv
```

### Output Directory Structure
```
Working Directory/
├── tracking_csv/           # Merged tracking files
│   └── *_tracking.csv
└── OFT_Results/
    ├── Results/           # Individual CSV results
    ├── Plots/             # Trajectory heatmaps
    ├── Temporal/          # Per-minute analysis
    └── Combined_Results.csv
```

## Usage

### Basic Analysis
```r
source("OFT_Main.R")
results <- run_analysis()
```

## Output Files

### Combined_Results.csv
Main output file containing all behavioral metrics:

- **Locomotion**: Total distance, speed, movement time
- **Zone metrics**: Time, distance, and transitions for each zone

### Trajectory Plots
Density heatmaps showing spatial occupancy with zone boundaries

### Temporal Analysis
Per-minute breakdown of distance traveled and zone occupancy

## Behavioral Metrics

The Combined_Results.csv file contains comprehensive behavioral metrics from the Open Field Test analysis. Each row represents one subject/video file.
Metric Descriptions
Basic Information

ID: Subject/file identifier (extracted from filename)
Duration_min: Total recording duration in minutes

General Locomotion Metrics

Total_Distance_cm: Total distance traveled during the entire session (cm)
Distance_Moving_cm: Distance traveled only during movement periods (cm)
Average_Speed_cm_s: Mean speed across entire session including stationary periods (cm/s)
Speed_Moving_cm_s: Mean speed calculated only during movement periods (cm/s)
Time_Moving_s: Total time spent in motion (seconds)
Time_Stationary_s: Total time spent stationary/immobile (seconds)
Percent_Moving: Percentage of session time spent moving (%)

Center Zone Metrics
The center zone typically represents 36% of the arena area (60% × 60% of dimensions)

Center_Time_s: Total time spent in the center zone (seconds)
Center_Distance_cm: Distance traveled within the center zone (cm)
Center_Transitions: Number of entries into the center zone (count)
Center_Percent: Percentage of total time spent in center zone (%)

Periphery Zone Metrics
The periphery is the outer ring of the arena, excluding corners

Periphery_Time_s: Total time spent in the periphery zone (seconds)
Periphery_Distance_cm: Distance traveled within the periphery zone (cm)
Periphery_Transitions: Number of entries into the periphery zone (count)
Periphery_Percent: Percentage of total time spent in periphery zone (%)

Corner Zone Metrics
The four corner regions of the arena

Corners_Time_s: Total time spent in all corner zones combined (seconds)
Corners_Distance_cm: Distance traveled within corner zones (cm)
Corners_Transitions: Number of entries into corner zones (count)
Corners_Percent: Percentage of total time spent in corners (%)

Anxiety-Related Indices

Center_Periphery_Ratio: Ratio of center time to periphery time (higher = less anxious behavior)
Thigmotaxis_Index: Wall-hugging behavior index; periphery time/(center + periphery time) (0-1 scale; higher = more thigmotactic/anxious behavior)

Interpretation Notes
Anxiety Assessment

Lower anxiety indicators: Higher center time/transitions, higher Center_Periphery_Ratio, lower Thigmotaxis_Index
Higher anxiety indicators: More time in periphery/corners, fewer center entries, higher Thigmotaxis_Index

Movement Analysis

Hyperactivity: High Total_Distance, high Percent_Moving
Hypoactivity: Low Total_Distance, high Time_Stationary
Movement efficiency: Compare Speed_Moving to Average_Speed

Zone Preferences

Exploration: Balanced time across zones with multiple transitions
Avoidance: Minimal center time with few transitions
Thigmotaxis: High periphery/corner time with minimal center exploration

Configuration Dependencies
These metrics are calculated based on your config settings:

movement_cutoff: Threshold for determining movement vs. stationary (default: 3 cm/s)
center_scale: Size of center zone (default: 0.6 = 60% of arena width/height)
periphery_scale: Size of periphery boundary (default: 0.8)
corner_scale: Size of corner zones (default: 0.25)
