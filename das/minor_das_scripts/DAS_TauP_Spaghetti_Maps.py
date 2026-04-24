import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

warnings.filterwarnings("ignore", category=UserWarning)

model = TauPyModel(model="ak135")


def load_coords(filepath):
    """Load and format cable coordinates, calculating cumulative distance."""
    df = pd.read_csv(filepath, sep=r'\s+', header=None,
        names=['lon', 'lat', 'cha', 'dep']).dropna()
    
    df['lon'] = df['lon'].apply(lambda x: x - 360 if x > 180 else x)
    
    R = 6371.0
    lat = np.radians(df['lat'].values)
    lon = np.radians(df['lon'].values)
    
    dlat = np.diff(lat)
    dlon = np.diff(lon)
    
    a = (np.sin(dlat / 2)**2 + 
         np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(dlon / 2)**2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    df['dist_km'] = np.insert(R * c, 0, 0).cumsum()
    return df


def analyze_cable_moveout(coords_df, event_dict):
    """Analyze travel time differences between start and end of the cable."""
    p1 = coords_df.iloc[0]
    p2 = coords_df.iloc[-1]
    results = []
    
    print(f"{'Event ID':<32} | {'Phase':<5} | "
          f"{'T_Start':<8} | {'T_End':<8} | {'Diff (s)':<8}")
    print("-" * 75)

    for eid, info in event_dict.items():
        z_km = info['dep'] / 1000.0
        dist1 = locations2degrees(info['lat'], info['lon'], p1['lat'], p1['lon'])
        dist2 = locations2degrees(info['lat'], info['lon'], p2['lat'], p2['lon'])
        
        p_phases = ["p", "P", "Pn", "Pg"]
        arr1_p = model.get_travel_times(z_km, dist1, phase_list=p_phases)
        arr2_p = model.get_travel_times(z_km, dist2, phase_list=p_phases)
        
        t1, t2, phase_used = None, None, None

        if len(arr1_p) > 0 and len(arr2_p) > 0:
            t1, t2 = arr1_p[0].time, arr2_p[0].time
            phase_used = "P"
        else:
            s_phases = ["s", "S", "Sn", "Sg"]
            arr1_s = model.get_travel_times(z_km, dist1, phase_list=s_phases)
            arr2_s = model.get_travel_times(z_km, dist2, phase_list=s_phases)
            if len(arr1_s) > 0 and len(arr2_s) > 0:
                t1, t2 = arr1_s[0].time, arr2_s[0].time
                phase_used = "S"

        if t1 is not None and t2 is not None:
            diff = abs(t1 - t2)
            results.append({'id': eid, 'phase': phase_used, 't1': t1, 't2': t2, 
                'diff': diff, 'dist1': dist1, 'z': z_km})
            print(f"{eid:<32} | {phase_used:<5} | "
                  f"{t1:>8.2f} | {t2:>8.2f} | {diff:>8.2f}")
        else:
            print(f"{eid:<32} | {'None':<5} | "
                  f"{'N/A':>8} | {'N/A':>8} | {'N/A':>8}")
    
    return results


def compute_event_arrivals(args):
    """Worker function to compute travel times for all channels."""
    ev, kkfls_coords, ev_data = args
    
    # This imports inside function to avoid multiprocessing issues
    from obspy.taup import TauPyModel  
    from obspy.geodetics import locations2degrees
    
    worker_model = TauPyModel(model="ak135")
    
    eid = ev['id']
    eq_lat = ev_data[eid]['lat']
    eq_lon = ev_data[eid]['lon']
    eq_depth = ev['z']
    
    if ev['phase'] == "P":
        target_phases = ["p", "P", "Pn", "Pg"]
    else:
        target_phases = ["s", "S", "Sn", "Sg"]
    
    arrival_times = []
    for _, row in kkfls_coords.iterrows():
        dist = locations2degrees(eq_lat, eq_lon, row['lat'], row['lon'])
        arrivals = worker_model.get_travel_times(eq_depth, dist,
            phase_list=target_phases)
        arrival_times.append(arrivals[0].time if arrivals else None)
        
    return eid, arrival_times


if __name__ == '__main__':
    event_data = {
        "2023-12-31-02-11-56.537000ML1.5": {
            "lat": 59.832, 
            "lon": -151.569, 
            "dep": 61400.0
        },
        "2023-12-10-06-49-27.068000ML1.4": {
            "lat": 59.929, 
            "lon": -152.260, 
            "dep": 73300.0
        },
        "2023-11-18-06-42-18.743000ML1.8": {
            "lat": 60.450, 
            "lon": -152.067, 
            "dep": 81100.0
        },
        "2023-10-26-02-36-59.823000ML1.3": {
            "lat": 59.976, 
            "lon": -152.272,
            "dep": 75100.0
        },
        "2023-10-03-19-50-59.443000ML1.4": {
            "lat": 60.090, 
            "lon": -152.540, 
            "dep": 92900.0
        },
        "2023-09-10-03-19-28.336000ML1.4": {
            "lat": 59.297, 
            "lon": -152.916, 
            "dep": 81200.0
        },
        "2023-08-18-01-19-46.681000ML1.3": {
            "lat": 59.845, 
            "lon": -150.637, 
            "dep": 31600.0
        },
        "2023-07-26-06-07-07.039000ML1.2": {
            "lat": 59.103, 
            "lon": -152.437, 
            "dep": 64100.0
        },
        "2023-07-03-09-42-42.845000ML1.5": {
            "lat": 60.133, 
            "lon": -151.232, 
            "dep": 46300.0
        },
        "2023-06-10-07-33-00.219000ML1.9": {
            "lat": 59.262, 
            "lon": -152.349, 
            "dep": 71500.0
        }
    }

    base_dir = "/Users/ed/research_code/das"
    coords_path = os.path.join(
        base_dir, 'das_coords_bathymetry', 'KKFLS_coords.xycz'
    )
    
    figures_dir = os.path.join(base_dir, 'das_figures')
    os.makedirs(figures_dir, exist_ok=True)
    pdf_out_path = os.path.join(figures_dir, 'KKFLS_predicted_moveouts.pdf')
    
    kkfls_coords = load_coords(coords_path)
    moveout_data = analyze_cable_moveout(kkfls_coords, event_data)
    events_to_process = moveout_data[:10]
    args_list = [(ev, kkfls_coords, event_data) for ev in events_to_process]

    print("\nCalculating arrival times (parallelized)...")
    with ProcessPoolExecutor() as ex:
        all_arrivals = list(ex.map(compute_event_arrivals, args_list))
    
    # Map event IDs back to their computed arrival arrays
    arrivals_dict = {eid: times for eid, times in all_arrivals}

    # Plotting and saving to PDF    
    with PdfPages(pdf_out_path) as pdf:
        for ev in events_to_process:
            eid = ev['id']
            eq_lat = event_data[eid]['lat']
            eq_lon = event_data[eid]['lon']
            eq_depth = ev['z']
            
            arrival_times = arrivals_dict[eid]
            channels = kkfls_coords['cha'].values
            
            fig, (ax_map, ax_time) = plt.subplots(1, 2, figsize=(16, 6))
            
            ax_map.plot(kkfls_coords['lon'], kkfls_coords['lat'], 'k-',
                linewidth=3, label='KKFLS Cable', alpha=0.5)
            ax_map.scatter(kkfls_coords['lon'].iloc[0],
                kkfls_coords['lat'].iloc[0], c='green', marker='s', s=100,
                label='Start (Ch 0)', zorder=5)
            ax_map.scatter(kkfls_coords['lon'].iloc[-1],
                kkfls_coords['lat'].iloc[-1], c='red', marker='s', s=100,
                label='End', zorder=5)
            ax_map.scatter(eq_lon, eq_lat, marker='*', s=300, c='gold',
                edgecolors='k', label=f'Hypocenter ({eq_depth:.1f}km depth)',
                zorder=6)
            
            ax_map.set_title(f"Event: {eid}")
            ax_map.set_xlabel("Longitude")
            ax_map.set_ylabel("Latitude")
            ax_map.legend(loc='best')
            ax_map.grid(True, alpha=0.3)
            ax_map.set_aspect('equal', adjustable='datalim')

            ax_time.plot(channels, arrival_times, color='blue', linewidth=2,
                alpha=0.7, label='Predicted Arrival')
            
            if arrival_times[0] is not None:
                ax_time.scatter(channels[0], arrival_times[0], c='green',
                    marker='s', s=100, label='Ch 0 Arrival', zorder=5)
            if arrival_times[-1] is not None:
                ax_time.scatter(channels[-1], arrival_times[-1], c='red',
                    marker='s', s=100, label='End Arrival', zorder=5)
                
            ax_time.set_title(f"Predicted {ev['phase']}-wave Moveout")
            ax_time.set_xlabel("Channel Number")
            ax_time.set_ylabel("Travel Time (s)")
            ax_time.legend(loc='best')
            ax_time.grid(True, alpha=0.3)

            plt.tight_layout()
            
            pdf.savefig(fig)
            plt.close(fig)
            
    print(f"\nSuccessfully saved all moveout figures to: {pdf_out_path}")