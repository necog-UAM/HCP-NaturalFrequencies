# HCP-NaturalFrequencies
# HCP MEG Natural Frequency Analysis Pipeline

Matlab code employed to conduct the analysis and create the figures of the following paper:

Herrera-Morueco,J.J., Stern, E., Arana, L., Capilla, A. [Preprint] Globally stable, locally flexible: Dynamic reconfiguration of brain natural frequencies during cognitive processing. 

The main function is BATCH_HCP
The main folder as well as subfolders should be added to the Matlab path.

## 📖 Project Overview
This repository contains an automated MATLAB pipeline for analyzing MEG data from the Human Connectome Project (HCP). The workflow extracts subject-specific and group-level natural frequency maps, performs voxel-wise statistical comparisons across experimental conditions, applies non-parametric permutation testing for multiple comparison correction, and generates publication-ready spatial maps and spectral visualizations.

## 🛠️ Setup & Dependencies
1. **Verify Dependencies First:** Before running the main pipeline, you should either execute `get_dependencies.m` and/or genpath/addpath the 'BATCH_HCP_Script_Dependencies' folder to ensure you can run the main pipeline (BATCH_HCP.m) without problems. This script (`get_dependencies.m`) uses MATLAB's static code analyzer to scan `BATCH_HCP.m`, verify that all required toolboxes (FieldTrip, Signal Processing, Statistics, etc.) and custom functions are accessible, and report any missing files.
2. The code has been developed using Matlab 2022b and fieldtrip-20230118.

## 📂 Configuration & Folder Management
⚠️ **Important:** The `BATCH_HCP.m` script contains hardcoded base directories. **To replicate results or run the pipeline on a new system, you must modify the path variables** (`p.rawpath`, `p.datapath`, `p.kmeanspath`, `p.FNAT`, etc.) at the beginning of **Part 1** and **Part 2** to match your local directory structure. You might need to change some other root directories along the script to accurately find the data and plots generated throughout the script. If you don't mind the already existing structure, the script Will generate most of the needed folders automatically.

📁 **Automatic Folder Creation:** You do not need to manually create output directories. Throughout execution, the script automatically generates the necessary subdirectories (e.g., `Kmeans/`, `FNAT_max/`, `img/MRICroN_nii/`, `img/test_img/`, etc.) using `mkdir()`. Ensure the base paths you specify have write permissions.

## 🧠 What `BATCH_HCP.m` Does
The pipeline is split into two self-contained sections that can be run independently. If Part 1 outputs already exist, you may skip directly to Part 2.

### 🔹 Part 1: Preprocessing & Natural Frequency Estimation
- **Subject & Task Validation:** Scans raw data directories and builds a validated subject-task matrix.
- **Source Reconstruction & Frequency Analysis:** Reconstructs source activity per session and computes frequency-domain data separated by experimental condition.
- **K-Means Clustering:** Pools power spectra across subjects, applies K-means clustering (`K=25`) to identify reproducible spectral profiles, and assigns each voxel to a cluster.
- **Natural Frequency Mapping:** Computes subject-specific natural frequency maps based on cluster similarity, aggregates them into group-level averages, and exports the final dataset as a structured MATLAB object (`fnat_struct`).
- ATTENTION: No preprocessing (bad channel, bad segment and artifact removal) was needed since HCP data comes already clean and mostly free of artifacts.

### 🔹 Part 2: Statistical Analysis & Hypothesis Testing
- **Voxel-Wise Comparisons:** Computes effect sizes (Cohen’s *d*, Cliff’s *δ*), T-sum statistics, and normality tests across predefined condition contrasts (e.g., `Motort_Fixation vs Left_hand`, `Wrkmem_0_back vs 2_back`, `Story vs Math`, `Rest vs Story`).
- **Permutation Testing:** Runs 1,000 non-parametric permutations to generate max-statistic null distributions, establishing robust significance thresholds (95th percentile) for multiple comparison correction.
- **Visualization & Export:** 
  - Generates thresholded NIfTI brain maps for MRIcroN visualization.
  - Extracts local maxima and ranks voxels by maximun effect size across frequencies.
  - Produces condition-specific spectral difference plots for significant and manually curated voxels.
  - Organizes all figures into structured task/condition subdirectories for manuscript preparation.

## 🚀 Execution Instructions
1. Open `BATCH_HCP.m` in MATLAB.
2. Update all `p.*` path variables to your local working directory.
3. Run the script section-by-section using `Ctrl + Enter` (or run sequentially). Each `%%` block is logically isolated and can be paused, inspected, or resumed.
4. Outputs will be saved in the `FNAT_max/` and `img/` directories specified in the path variables.

## 📝 Notes
- **Processing Time:** Permutation testing and high-resolution figure export are computationally intensive. Consider running Part 2 on a machine with ≥16 GB RAM.
- **Dynamic Functions:** MATLAB's static dependency checker cannot resolve functions called via `eval`, `feval`, or string-based handles. Ensure all `HCP*_*.m` and `Omega_*.m` files are in the added paths.
- **License Compliance:** Toolbox files (e.g., FieldTrip, Statistics) are referenced but not redistributed. Ensure your MATLAB installation includes required toolboxes.
