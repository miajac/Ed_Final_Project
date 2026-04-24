import os
import json
import warnings
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use('Agg') 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import hilbert

from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.clients.fdsn import Client

import dascore as dc
import dascore.proc as dcp

warnings.filterwarnings("ignore")

base_project_dir = os.getcwd()
records_path = os.path.join(base_project_dir,
                            "das_records",
                            "good-events-3.2-up")
coords_file = os.path.join(base_project_dir,
                           "das_coords_bathymetry",
                           "TERRA_coords.xycz")
figures_dir = os.path.join(base_project_dir, "das_figures")
cache_file = os.path.join(base_project_dir, "event_metadata_cache.json")
os.makedirs(figures_dir, exist_ok=True)
pdf_path = os.path.join(figures_dir, "TERRA_good_events_arrivals.pdf")
csv_path = os.path.join(figures_dir, "TERRA_event_amplitudes.csv")

good_terra_events = [
    "ak0237eejw69_TERRA.h5", "ak023aw5mbdk_TERRA.h5", "ak023eccchzy_TERRA.h5",
    "ak0237q2shdo_TERRA.h5", "ak023bebgmhd_TERRA.h5", "ak023em715sv_TERRA.h5",
    "ak02381ibekf_TERRA.h5", "ak023bhlw02w_TERRA.h5", "ak023f7eyaqg_TERRA.h5",
    "ak0238ghnzxp_TERRA.h5", "ak023bkx215x_TERRA.h5", "ak023fhgggc6_TERRA.h5",
    "ak0238qkcxek_TERRA.h5", "ak023btef8mo_TERRA.h5", "ak023fnzshe1_TERRA.h5",
    "ak0239af45c3_TERRA.h5", "ak023bzqw7a7_TERRA.h5", "ak023frilkvn_TERRA.h5",
    "ak0239lyp68s_TERRA.h5", "ak023c3206y0_TERRA.h5", "ak023g35jxin_TERRA.h5",
    "ak0239qzbmym_TERRA.h5", "ak023cgc5fmi_TERRA.h5", "ak023gbgys2j_TERRA.h5",
    "ak0239saxy95_TERRA.h5", "ak023ctiuyia_TERRA.h5", "ak023gjh7z4b_TERRA.h5",
    "ak0239vxdtm6_TERRA.h5", "ak023d3dyqv0_TERRA.h5", "ak023godcr3i_TERRA.h5",
    "ak023a0j9eo0_TERRA.h5", "ak023gqcxl3z_TERRA.h5", "ak023a7ds9th_TERRA.h5",
    "ak023dif8i7c_TERRA.h5", "ak023djxyhod_TERRA.h5",
    "us6000lggw_TERRA.h5", "ak023ah9zcn5_TERRA.h5", "ak023dk7iyjo_TERRA.h5",
    "ak023asybefb_TERRA.h5", "ak023ds94esw_TERRA.h5"
]


def get_and_cache_metadata(file_list):
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as j:
            return json.load(j)

    usgs_client = Client("USGS")
    data = {}
    for f in file_list:
        eid = f.split('_')[0]
        try:
            cat = usgs_client.get_events(eventid=eid)
            origin = cat[0].origins[0]
            data[f] = {"lat": origin.latitude,
                       "lon": origin.longitude,
                       "dep": origin.depth,
                       "time": str(origin.time)}
        except Exception:
            if eid == "ak023gjh7z4b":
                data[f] = {"lat": 60.335,
                           "lon": -151.754,
                           "dep": 71300.0,
                           "time": "2023-12-26T06:30:02.910"}
            elif eid == "ak023d3dyqv0":
                data[f] = {"lat": 60.320,
                           "lon": -150.699,
                           "dep": 79300.0,
                           "time": "2023-10-12T03:02:49.457"}
    with open(cache_file, 'w') as j:
        json.dump(data, j, indent=4)
    return data


def load_coords(filepath):
    df = pd.read_csv(filepath,
                     sep=r'\s+',
                     header=None,
                     names=['lon', 'lat', 'cha', 'dep']).dropna()
    df['lon'] = np.where(df['lon'] > 180, df['lon'] - 360, df['lon'])
    return df


def process_single_event(args):
    """Worker function to process a single event independently"""
    fname, info, records_path, terra_coords, figures_dir = args

    try:
        model = TauPyModel(model="ak135")

        patch = dc.read(os.path.join(records_path, fname))[0]
        patch = dcp.pass_filter(patch, time=(0.3, 10.0))
        patch = dcp.decimate(patch, distance=20)
        data_array = patch.get_array()
        envelope = np.abs(hilbert(data_array, axis=0))

        origin_time = np.datetime64(info['time'].replace('Z', ''))
        times_sec = (
            patch.coords.get_array('time') - origin_time
        ) / np.timedelta64(1, 's')
        decimated_coords = terra_coords.iloc[::20]
        limit = min(len(decimated_coords), data_array.shape[1])
        z_km = info['dep'] / 1000.0 

        p_win_low, p_win_high = [], []
        s_win_low, s_win_high = [], []
        max_amps = []
        p_max_amps = []
        pick_times, pick_chans = [], []
        slant_km_list = []
        inside_s_count = 0
        valid_s_channels = 0

        for i in range(limit):
            row = decimated_coords.iloc[i]
            dist_deg = locations2degrees(
                info['lat'], info['lon'], row['lat'], row['lon']
            )
            slant_km = np.hypot(dist_deg * 111.19, z_km)
            slant_km_list.append(slant_km)

            arrivals = model.get_travel_times(
                z_km, dist_deg, phase_list=["p", "P", "s", "S"]
            )
            p_arr = [a.time for a in arrivals if a.name.lower() in ['p', 'pg']]
            s_arr = [a.time for a in arrivals if a.name.lower() in ['s', 'sg']]
            p_t = p_arr[0] if p_arr else np.nan
            s_t = s_arr[0] if s_arr else np.nan

            pw_l, pw_h = p_t - 2, p_t + 4
            s_low, s_high = s_t - 6, s_t + 15

            p_win_low.append(pw_l)
            p_win_high.append(pw_h)
            s_win_low.append(s_low)
            s_win_high.append(s_high)

            if not np.isnan(p_t):
                p_mask = (times_sec >= pw_l) & (times_sec <= pw_h)
                p_max_amps.append(
                    np.max(envelope[p_mask, i]) if np.any(p_mask) else 0)
            else:
                p_max_amps.append(0)

            max_val = np.max(envelope[:, i])
            max_idx_time = times_sec[np.argmax(envelope[:, i])]
            max_amps.append(max_val)
            pick_times.append(max_idx_time)
            pick_chans.append(i)

            if not np.isnan(s_t):
                valid_s_channels += 1
                if s_low <= max_idx_time <= s_high:
                    inside_s_count += 1

        if valid_s_channels > 0:
            perc_inside = (inside_s_count / valid_s_channels) * 100
        else:
            perc_inside = 0

        min_slant_idx = np.argmin(slant_km_list)
        min_slant_val = slant_km_list[min_slant_idx]

        stat_dict = {
            "event_id": fname.split('_')[0],
            "perc_inside_s": perc_inside,
            "mean_p_max": np.mean(p_max_amps),
            "mean_total_max": np.mean(max_amps),
            "min_slant_km": min_slant_val
        }

        # Plotting starts here
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10),
                                       gridspec_kw={'width_ratios': [3, 1]})
        vmax = np.nanpercentile(np.abs(data_array), 99.5) or 1.0
        chans = np.arange(limit)

        for i in range(limit):
            ax1.plot(times_sec,
                     i + (data_array[:, i] / vmax) * 2.5,
                     color='black',
                     lw=0.4,
                     alpha=0.3)

        ax1.fill_betweenx(chans, p_win_low, p_win_high, color='blue',
                          alpha=0.1, label='P-Window (-2s/+4s)')
        ax1.fill_betweenx(chans, s_win_low, s_win_high, color='red',
                          alpha=0.1, label='S-Window (-6s/+15s)')
        ax1.scatter(pick_times, pick_chans, color='magenta', s=15,
                    alpha=0.2, label='Max Energy', zorder=5)

        title_string = (f"ID: {fname.split('_')[0]} | Depth: {z_km:.1f} km | "
                     f"Straightline Distance: {min_slant_val:.2f} km\n"
                     f"{perc_inside:.1f}% of Max Amplitudes within S-window")
        ax1.set_title(title_string)
        ax1.set_xlabel("Seconds from Origin")
        ax1.set_ylabel("Channel Index (Decimated x20)")
        ax1.invert_yaxis()
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.1)

        ax2.plot(max_amps, chans, color='magenta', lw=1.5,
                 label='Global Max (Abs)')
        ax2.set_title("Max Absolute Amplitude Profile")
        ax2.set_xlabel("Amplitude")
        ax2.set_xscale('log')
        ax2.set_xlim(-1, 1000)
        ax2.invert_yaxis()
        ax2.legend(fontsize='small')
        ax2.grid(True, which='both', alpha=0.2)

        plt.tight_layout()
        temp_fig_path = os.path.join(figures_dir,
                                     f"temp_{fname.split('_')[0]}.png")
        fig.savefig(temp_fig_path, dpi=200)
        plt.close(fig)

        return temp_fig_path, stat_dict, perc_inside

    except Exception as e:
        print(f"[ERROR] {fname}: {e}")
        return None, None, None


if __name__ == '__main__':
    event_metadata = get_and_cache_metadata(good_terra_events)
    terra_coords = load_coords(coords_file)

    event_stats = []
    s_window_percentages = []
    args_list = [(fname, info, records_path, terra_coords, figures_dir)
                 for fname, info in event_metadata.items()]

    print(f"Parallel processing {len(args_list)} events...")

    # Execute in parallel
    with ProcessPoolExecutor() as ex:
        results = list(ex.map(process_single_event, args_list))

    # Assemble PDF from saved figures afterward
    with PdfPages(pdf_path) as pdf:
        for res in results:
            temp_path, stat, perc = res

            if temp_path is not None:
                fig = plt.figure(figsize=(16, 10))
                img = plt.imread(temp_path)
                plt.imshow(img)
                plt.axis('off')
                pdf.savefig(fig, bbox_inches='tight', pad_inches=0)
                plt.close(fig)

                if os.path.exists(temp_path):
                    os.remove(temp_path)

                event_stats.append(stat)
                s_window_percentages.append(perc)

    if s_window_percentages:
        mean_val = np.mean(s_window_percentages)
        median_val = np.median(s_window_percentages)
        min_val = np.min(s_window_percentages)
        max_val = np.max(s_window_percentages)

        print("\n% Max S-Amplitudes Summary")
        print(f"Mean:   {mean_val:.2f}%")
        print(f"Median: {median_val:.2f}%")
        print(f"Range:  {min_val:.2f}% to {max_val:.2f}%")

    df_results = pd.DataFrame(event_stats)
    df_results.to_csv(csv_path, index=False)

    print(f"\nPDF saved to {pdf_path}")
    print(f"CSV saved to {csv_path}")