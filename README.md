# research_code
My scripts, files, etc. from my research.  Current subfolder is das.  
Revised from local.  

What:
das contains all of my DAS (Distributed Acoustic Sensing) research.  This is, 
in other words, seismology conducted using fiber-optic cables as arrays of 
seismometers.  As this folder is currently in use (and more alarmingly, in use 
by the author) it is a bit messy, but here are the main features: 
    das_coords_bathymetry
        This is where the .xcyz bathymetry files for the cable coordinates are 
        stored.  

    das_freq_plots
        This contains a bunch of early frequency filtering trials of the DAS 
        data.  The code to generate these figures exists (in some form) in 
        minor_das_scripts/Uberfigure.ipynb.  

    minor_das_scripts
        This folder is where various scripts not currently in use are stored.  
        There is often a figure or two hanging out in this folder due to my 
        organizational shortcomings.
    
    MS-DAS.ipynb
        Named to convey the proper pronunciation of DAS (with the A in 
        Acoustic) by likening it to the operating system MS-DOS, this is the 
        main active research environment, containing multiple scripts within 
        its Jupyter Notebook code blocks.  Over time I break the inactive or 
        mostly-finalized scripts out of MS-DAS.ipynb and into separate notebook 
        files in /das/minor-das-scripts/ and rework them into functioning in 
        their new location.  The current MS-DAS file has a few minor scripts in 
        it along with the main functionality script focusing on plotting 
        calculated arrival time windows for multiple events over observed 
        arrival time data.

How:


Why: