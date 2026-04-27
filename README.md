# DAS Project

This project contains all research code for DAS (Distributed Acoustic Sensing) seismology conducted using fiber-optic cables as arrays of seismometers. Two cables are used: **TERRA** and **KKFLS**, both deployed in lower Cook Inlet, Alaska.

---
## Project Details:
**Task #1:** Parallelization/Concurrency  
**Task #2:** State Saving

## Installation

**Prerequisites:** Conda, Miniconda, Mamba, or equivalent.

```bash
# 1. Clone / download the project, then cd into it
cd Final_Project

# 2. Create the environment (first time only — takes a few minutes)
conda env create -f environment.yaml

# 3. Activate it
conda activate dasenv

# 4. Verify the key packages loaded correctly
python -c "import dascore, obspy, cartopy; print('All good!')"
```

To update the environment after changing `environment.yaml`:
```bash
conda env update -f environment.yaml --prune
```

To fully remove the environment:
```bash
conda deactivate
conda env remove -n dasenv
```

---

## Directory Layout

```
Final_Project/
├── README.md                         ← this file
├── environment.yaml                  ← conda environment definition
└── das/                              ← main project directory
    ├── ms-das.py                     ← primary research pipeline script
    ├── event_metadata_cache.json     ← event lat/lon/depth cache used
    │                                    by mapview_TERRA_events.py
    ├── das_coords_bathymetry/        ← cable geometry and seafloor data
    │   ├── TERRA_coords.xycz         ← TERRA cable coordinates + depth
    │   ├── KKFLS_coords.xycz         ← KKFLS cable coordinates + depth
    │   └── PW24.xyz                  ← regional seafloor bathymetry grid
    ├── das_records/
    │   └── good-events-3.2-up/
    │       └── <event_id>_TERRA.h5   ← HDF5 event files (user-downloaded due 
    │                                    to extremely large filesizes)
    ├── das_figures/                  ← auto-created on first run,
    │                                    stores figures & caches
    └── minor_das_scripts/
        ├── URL_table.py
        ├── USGS_file_renaming.py
        ├── das_bathymetry.py
        ├── raw_noise_files.py
        ├── mapview_TERRA_events.py
        ├── ETI2_slide_maker.py
        ├── noise_analysis.py
        ├── taup_spaghetti_maps.py
        ├── uberfigure_TERRA_only.py
        └── uberfigure_combined.py
```

---

## Scripts

### `ms-das.py` — Main Research Pipeline
The primary working script. Loads DAS event files, computes TauP-predicted P and S arrival windows, plots observed strain data against those windows, and exports results as a multi-page PDF and CSV. This is the main testing ground for core analysis functionality, particularly for selecting time windows used to calculate wave amplitudes. The name is a nod to MS-DOS — because the correct pronunciation of DAS rhymes with the OS.

### `minor_das_scripts/`
                
                    

| Script | What it does | How | Why |
|---|---|---|---|
| `URL_table.py` | Prints a formatted table of USGS events for a given URL query | Queries the USGS FDSN event service for M ≥ 3.2 earthquakes within 150 km of lower Cook Inlet, Jun–Dec 2023 | Useful as a manual check on query results and as a reference when renaming downloaded files |
| `USGS_file_renaming.py` | Renames downloaded DAS files from the old date-based naming scheme to USGS event ID-based names | Maps old filenames to their USGS IDs (e.g. `ak0237eejw69_TERRA.h5`) using a hardcoded lookup table; set `DRY_RUN = False` to apply | Standardizes filenames to match keys in `event_metadata_cache.json` |
| `das_bathymetry.py` | 3-D wireframe and surface bathymetry map with both cable tracks overlaid | Reads `.xycz` coordinate files and renders a 3-D matplotlib figure of the seafloor with KKFLS and TERRA routes plotted | Visualizes the physical layout of both cables relative to the seafloor of lower Cook Inlet |
| `raw_noise_files.py` | Per-event raw noise profiles for both cables | Loads the final 45 s of each event file and plots absolute median phase/strain amplitude vs. distance along the cable | Used to characterize the raw background noise floor of each cable independently |
| `mapview_TERRA_events.py` | Albers Equal Area map of the TERRA cable and earthquake epicenters | Reads `event_metadata_cache.json` and the TERRA coordinate file, then renders a cartopy map with an Alaska overview inset | Originally produced for a conference presentation |
| `ETI2_slide_maker.py` | Single-channel waveform plot for a specific event and distance | Loads one event file, isolates the channel nearest a target distance, high-pass filters above 1 Hz, and saves a PDF | A quick script for generating introductory waveform figures for presentations |
| `noise_analysis.py` | Heatmap of pre-event noise floors sorted by earthquake magnitude | Extracts median absolute strain amplitude for the first 5 s and last 20 s of each event, fetches live magnitudes from the USGS API, and renders a log-scale heatmap with bathymetry | Active working script for understanding noise patterns and deviations across different event magnitudes |
| `taup_spaghetti_maps.py` | Per-channel TauP moveout curves for each event, saved to a multi-page PDF | Computes P and S travel times from the ak135 model for every cable channel and plots predicted arrival time vs. channel number alongside a map of the event geometry | Visualizes expected moveout across the cable — the resulting plots have a tendency to look like spaghetti |
| `uberfigure_TERRA_only.py` | Full TERRA reference wiggle plot with labelled zoom insets and a red-highlighted channel | Renders a decimated reference wiggle for the entire TERRA cable, then draws dashed zoom boxes and expands each as a separate high-resolution sub-panel | Produced the primary visual outputs during early data exploration |
| `uberfigure_combined.py` | Dual-cable reference wiggle combining KKFLS (negative distances) and TERRA (positive distances) with zoom insets for both | Merges the two cables into a single sorted distance axis, renders the combined reference wiggle, and produces labeled zoom panels for user-defined regions of both cables | Used to compare waveforms across both cables in a single figure |

---

## Typical Run Order

```bash
conda activate dasenv

# 1. Browse available events and identify files to download
python minor_das_scripts/URL_table.py

# 2. Download files from Data Archives (under Event Data) and place them
# in the folder das_records/good-events-3.2-up/
# https://dasway.ess.washington.edu/gci

# 3. Rename downloaded files to USGS ID format (set DRY_RUN = False first)
python minor_das_scripts/USGS_file_renaming.py

# 4. Explore geometry and data quality
python minor_das_scripts/das_bathymetry.py
python minor_das_scripts/raw_noise_files.py
python minor_das_scripts/mapview_TERRA_events.py

# 5. Main analysis
python das/ms-das.py
python minor_das_scripts/noise_analysis.py
python minor_das_scripts/taup_spaghetti_maps.py
python minor_das_scripts/uberfigure_TERRA_only.py
python minor_das_scripts/uberfigure_combined.py
python minor_das_scripts/ETI2_slide_maker.py
```

---

## Notes

- **Caching** — `ms-das.py`, `noise_analysis.py`, and `taup_spaghetti_maps.py` all write `.npz` and `.json` cache files into `das_figures/`. Delete the relevant cache files to force a full reprocess.
- **`das_records/`** — HDF5 event files are not included in the repository and must be downloaded separately from the DAS data portal at https://dasway.ess.washington.edu/gci and placed in `das_records/good-events-3.2-up/`.