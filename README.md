# research_code
My scripts, files, etc. from my research, and for my GEOS694 final project. Current subfolder is `das`.  

## What

`das` contains all of my DAS (Distributed Acoustic Sensing) research. This is, in other words, seismology conducted using fiber-optic cables as arrays of seismometers. As this folder is currently in use it is a bit messy, but here are the main features:

- `das_coords_bathymetry`  
  This is where the `.xcyz` bathymetry files for the cable coordinates are stored.

- `das_freq_plots`  
  This contains early frequency-filtering trials of the DAS data. The code to generate these figures exists, in some form, in `minor_das_scripts/Uberfigure.ipynb`.

- `minor_das_scripts`  
  This folder is where various scripts not currently in use are stored. There is often a figure or two in this folder as well.
  -  Within `minor_das_scripts` are the majority of scripts under this project.  
  -  `DAS_Bathymetry.ipynb`
        
  -  `ETI2.0 Slide Maker.ipynb`
        
  -  `DAS_Raw_Noise_Files.ipynb`
        
  -  `event_metadata_cache.json`
        
  -  `das_seafloor_bathymetry.png`
        
  -  `Mapview_TERRA_events.ipynb`
        
  -  `DAS_TauP_Spaghetti_Maps.ipynb`
        
  -  `Noise_Analysis.ipynb`
        
  -  `DAS_URL_table.ipynb`
        
  -  `Uberfigure.ipynb`
        
  -  `DAS_USGS_File_Renaming.ipynb`

#
- `MS-DAS.ipynb`  
  Named to convey the proper pronunciation of DAS by likening it to the operating system MS-DOS, this is the main active research environment, containing multiple scripts within its Jupyter Notebook code blocks. Over time, inactive or mostly finalized scripts are broken out of `MS-DAS.ipynb` and into separate notebook files in `das/minor_das_scripts/`, then reworked into functioning in their new location. The current `MS-DAS.ipynb` file contains a few minor scripts along with the main functionality script focusing on plotting calculated arrival-time windows for multiple events over observed arrival-time data.

## How

## Why
