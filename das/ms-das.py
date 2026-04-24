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
metadata_cache = os.path.join(base_project_dir, "event_metadata_cache.json")

# Save-state cache: two files, .npz and .json
# Stores all derived arrays to skip processing on repeats.
master_cache_npz = os.path.join(figures_dir, "TERRA_master_arrays.npz")
master_cache_json = os.path.join(figures_dir, "TERRA_master_stats.json")

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
    if os.path.exists(metadata_cache):
        with open(metadata_cache, 'r') as j:
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
                data[f] = {"lat": 60.335, "lon": -151.754,
                           "dep": 71300.0, "time": "2023-12-26T06:30:02.910"}
            elif eid == "ak023d3dyqv0":
                data[f] = {"lat": 60.320, "lon": -150.699,
                           "dep": 79300.0, "time": "2023-10-12T03:02:49.457"}
    with open(metadata_cache, 'w') as j:
        json.dump(data, j, indent=4)
    return data


def load_coords(filepath):
    df = pd.read_csv(filepath, sep=r'\s+', header=None,
                     names=['lon', 'lat', 'cha', 'dep']).dropna()
    df['lon'] = np.where(df['lon'] > 180, df['lon'] - 360, df['lon'])
    return df


def _plot_event_figure(fname, data_array, times_sec, p_win_low, p_win_high,
                       s_win_low, s_win_high, max_amps, p_max_amps,
                       pick_times, pick_chans, slant_km_list,
                       vmax, z_km, perc_inside, limit, figures_dir):
    """Shared plotting helper that saves a temporary PNG."""
    event_id = fname.split('_')[0]
    chans = np.arange(limit)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10),
                                   gridspec_kw={'width_ratios': [3, 1]})

    for i in range(limit):
        ax1.plot(times_sec,
                 i + (data_array[:, i] / vmax) * 2.5,
                 color='black', lw=0.4, alpha=0.3)

    ax1.fill_betweenx(chans, p_win_low, p_win_high, color='blue',
                      alpha=0.1, label='P-Window (-2s/+4s)')
    ax1.fill_betweenx(chans, s_win_low, s_win_high, color='red',
                      alpha=0.1, label='S-Window (-6s/+15s)')
    ax1.scatter(pick_times, pick_chans, color='magenta', s=15,
                alpha=0.2, label='Max Energy', zorder=5)

    min_slant_val = slant_km_list[np.argmin(slant_km_list)]
    title_string = (f"ID: {event_id} | Depth: {z_km:.1f} km | "
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
    temp_fig_path = os.path.join(figures_dir, f"temp_{event_id}.png")
    fig.savefig(temp_fig_path, dpi=200)
    plt.close(fig)
    return temp_fig_path


def process_single_event(args):
    """
    Worker function to process a single event, utilizing passed-in caches.
    If 'cached_arrays' is populated, it skips directly to plotting.
    Otherwise, it processes the data and returns the arrays for saving.
    """
    (fname, info, records_path, terra_coords, figures_dir, 
     cached_stat, cached_arrays) = args
    
    event_id = fname.split('_')[0]

    if cached_stat is not None and cached_arrays is not None:
        try:
            temp_path = _plot_event_figure(
                fname,
                data_array    = cached_arrays['data_array'],
                times_sec     = cached_arrays['times_sec'],
                p_win_low     = cached_arrays['p_win_low'],
                p_win_high    = cached_arrays['p_win_high'],
                s_win_low     = cached_arrays['s_win_low'],
                s_win_high    = cached_arrays['s_win_high'],
                max_amps      = cached_arrays['max_amps'],
                p_max_amps    = cached_arrays['p_max_amps'],
                pick_times    = cached_arrays['pick_times'],
                pick_chans    = cached_arrays['pick_chans'],
                slant_km_list = cached_arrays['slant_km_list'],
                vmax          = float(cached_arrays['vmax']),
                z_km          = float(cached_arrays['z_km']),
                perc_inside   = float(cached_arrays['perc_inside']),
                limit         = int(cached_arrays['limit']),
                figures_dir   = figures_dir,
            )
            return temp_path, cached_stat, cached_stat['perc_inside_s'], None
            
        except Exception as e:
            print(f"Cache Error {event_id}: {e} — reprocessing...")


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
        max_amps      = []
        p_max_amps    = []
        pick_times    = []
        pick_chans    = []
        slant_km_list = []
        inside_s_count  = 0
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

        perc_inside = (inside_s_count / valid_s_channels * 100
                       if valid_s_channels > 0 else 0)
        min_slant_val = slant_km_list[np.argmin(slant_km_list)]
        vmax = np.nanpercentile(np.abs(data_array), 99.5) or 1.0

        stat_dict = {
            "event_id":       event_id,
            "perc_inside_s":  perc_inside,
            "mean_p_max":     np.mean(p_max_amps),
            "mean_total_max": np.mean(max_amps),
            "min_slant_km":   min_slant_val,
        }

        arrays_dict = {
            'data_array':    data_array,
            'times_sec':     times_sec,
            'p_win_low':     np.array(p_win_low),
            'p_win_high':    np.array(p_win_high),
            's_win_low':     np.array(s_win_low),
            's_win_high':    np.array(s_win_high),
            'max_amps':      np.array(max_amps),
            'p_max_amps':    np.array(p_max_amps),
            'pick_times':    np.array(pick_times),
            'pick_chans':    np.array(pick_chans),
            'slant_km_list': np.array(slant_km_list),
            'vmax':          np.array([vmax]),
            'z_km':          np.array([z_km]),
            'perc_inside':   np.array([perc_inside]),
            'limit':         np.array([limit]),
        }

        temp_path = _plot_event_figure(
            fname,
            data_array    = data_array,
            times_sec     = times_sec,
            p_win_low     = arrays_dict['p_win_low'],
            p_win_high    = arrays_dict['p_win_high'],
            s_win_low     = arrays_dict['s_win_low'],
            s_win_high    = arrays_dict['s_win_high'],
            max_amps      = arrays_dict['max_amps'],
            p_max_amps    = arrays_dict['p_max_amps'],
            pick_times    = arrays_dict['pick_times'],
            pick_chans    = arrays_dict['pick_chans'],
            slant_km_list = arrays_dict['slant_km_list'],
            vmax          = vmax,
            z_km          = z_km,
            perc_inside   = perc_inside,
            limit         = limit,
            figures_dir   = figures_dir,
        )
        
        return temp_path, stat_dict, perc_inside, arrays_dict

    except Exception as e:
        print(f"[ERROR] {fname}: {e}")
        return None, None, None, None


if __name__ == '__main__':
    event_metadata = get_and_cache_metadata(good_terra_events)
    terra_coords = load_coords(coords_file)

    global_stats = {}
    if os.path.exists(master_cache_json):
        with open(master_cache_json, 'r') as f:
            global_stats = json.load(f)

    global_arrays = {}
    if os.path.exists(master_cache_npz):
        with np.load(master_cache_npz, allow_pickle=False) as c:
            for k in c.files:
                global_arrays[k] = c[k]

    args_list = []
    keys_needed = [
        'data_array', 'times_sec', 'p_win_low', 'p_win_high',
        's_win_low', 's_win_high', 'max_amps', 'p_max_amps',
        'pick_times', 'pick_chans', 'slant_km_list',
        'vmax', 'z_km', 'perc_inside', 'limit'
    ]

    for fname, info in event_metadata.items():
        event_id = fname.split('_')[0]
        cached_stat = global_stats.get(event_id)
        cached_arrays = None

        if cached_stat is not None:
            try:
                cached_arrays = {k: global_arrays[f"{event_id}_{k}"] 
                                 for k in keys_needed}
            except KeyError:
                cached_arrays = None
                cached_stat = None  # Force reprocess if NPZ is missing arrays

        args_list.append((fname, info, records_path, terra_coords, 
                          figures_dir, cached_stat, cached_arrays))

    print(f"Processing {len(args_list)} events..")

    with ProcessPoolExecutor() as ex:
        results = list(ex.map(process_single_event, args_list))

    event_stats = []
    s_window_percentages = []
    needs_save = False

    with PdfPages(pdf_path) as pdf:
        for res in results:
            temp_path, stat, perc, new_arrays = res

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

            if new_arrays is not None:
                event_id = stat['event_id']
                global_stats[event_id] = stat
                for k, v in new_arrays.items():
                    global_arrays[f"{event_id}_{k}"] = v
                needs_save = True

    # Caching
    if needs_save:
        print("\nUpdating master cache files...")
        np.savez_compressed(master_cache_npz, **global_arrays)
        with open(master_cache_json, 'w') as f:
            json.dump(global_stats, f, indent=4)

    if s_window_percentages:
        mean_val   = np.mean(s_window_percentages)
        median_val = np.median(s_window_percentages)
        min_val    = np.min(s_window_percentages)
        max_val    = np.max(s_window_percentages)

    df_results = pd.DataFrame(event_stats)
    df_results.to_csv(csv_path, index=False)

    print(f"\nPDF saved to {pdf_path}")
    print(f"CSV saved to {csv_path}")
    print(f"Master cache maintained at {master_cache_npz}")

    print("\n% Max S-Amplitudes Summary")
    print(f"Mean:   {mean_val:.2f}%")
    print(f"Median: {median_val:.2f}%")
    print(f"Range:  {min_val:.2f}% to {max_val:.2f}%")