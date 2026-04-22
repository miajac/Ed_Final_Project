# research_code
My scripts, files, etc. from my research, and for my GEOS694 final project. Current subfolder is `das`.  

## What, Why, and How:

`das` contains all of my DAS (Distributed Acoustic Sensing) research. This is, in other words, seismology conducted using fiber-optic cables as arrays of seismometers. As this folder is currently in use it is a bit messy, but here are the main features:

- `das_coords_bathymetry`  
  This is where the `.xcyz` bathymetry files for the cable coordinates are stored.

- `das_freq_plots`  
  This contains early frequency-filtering trials of the DAS data. The code to generate these figures exists, in some form, in `minor_das_scripts/Uberfigure.ipynb`.

- `minor_das_scripts`  
  This folder is where various scripts not currently in use are stored. There is often a figure or two in this folder as well.
    - Within `minor_das_scripts` are the majority of scripts under this project.  

    -  `DAS_Bathymetry.ipynb` generates local bathymetry to `das_seafloor_bathymetry.png`, and reads the files in `das_coords_bathymetry/` to render a 3D plot of the seafloor with both the KKFL-S and TERRA cable routes mapped onto it.

    - `ETI2.0 Slide Maker.ipynb` is a quickly-made script to create introductory figures for a conference presentation.  It loads a single event file to create figures of a few waveforms over a small distance.

    - `DAS_Raw_Noise_Files.ipynb` plots absolute median strain across channels in both cables for a set of low-magnitude events, loading the final 45 seconds of each 2-minute long event and plotting their noise profiles.

    - `event_metadata_cache.json` is a JSON cache of seismic events spanning June–December 2023, with the USGS event ID serving as the key, and storing latitude, longitude, depth, and origin time.  This JSON is used in various scripts throughout the project to avoid repeated API calls.

    - `das_seafloor_bathymetry.png` is the output of `DAS_Bathymetry.ipynb`, and shows the two cables (KKFL-S in Red, TERRA in Blue) overlaid on a 3D seafloor bathymetry map of the lower Cook Inlet region of Alaska.  Properly, this file would be in `das_figures/` but I haven't rerouted that script's output.

    - `Mapview_TERRA_events.ipynb` plots the TERRA cable and all cached seismic events onto a regional Albers Equal Area map using Cartopy, made for a conference presentation.

    - `DAS_TauP_Spaghetti_Maps.ipynb` contains two related scripts, the first uses TauP to compute P and S-wave travel times between the first and last cable channels for each event with the ak135 velocity model and outputs a summary table.  The second generates plots of predicted arrival times across all cable channels for every event simultaneously, and these outputs look a bit like spaghetti thrown at a wall, hence the filename.

    - `Noise_Analysis.ipynb` is the primary noise characterization script, and has become a 'working file' similar to `MS-DAS.ipynb`.  As it is in active development, this description may not be fully up to date.  For each event, this script extracts median absolute strain amplitude over the first 5 seconds and last 20 seconds for comparison purposes, then displays noise amplitude readings for each event, sorted by event magnitude, which are fetched live from the USGS API.  Contains multiple plotting variants.

    - `DAS_URL_table.ipynb` Queries the USGS FDSN event service for M3.2+ earthquakes within 150 km of the cable center during 2023 and prints a formatted table of origin time, magnitude, coordinates, and depth. Used for manually identifying events to download from the DAS data portal, which does not support automated crawling.
    
    - `DAS_USGS_File_Renaming.ipynb` renames raw DAS files from the website to USGS event ID format (e.g. `ak0237eejw69_TERRA.h5`) to match `event_metadata_cache.json`.  This one might not be terribly useful as it was only for renaming from my initial scheme.

    - `Uberfigure.ipynb` produces a series of wiggle plots of the complete DAS cable waveforms, with zoomed-in subplots.  This gave the primary outputs during early data exploration days.  Resulting figures are saved to `das_freq_plots/`.

#
- `MS-DAS.ipynb`  
  Named to convey the proper pronunciation of DAS by likening it to the operating system MS-DOS, this is the main active research environment, containing multiple scripts within its Jupyter Notebook code blocks. Over time, inactive or mostly finalized scripts are broken out of `MS-DAS.ipynb` and into separate notebook files in `das/minor_das_scripts/`, then reworked into functioning in their new location. The current `MS-DAS.ipynb` file contains a few minor scripts along with the main functionality script focusing on plotting calculated arrival-time windows for multiple events over observed arrival-time data.
