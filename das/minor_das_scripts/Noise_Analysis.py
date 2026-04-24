# Noise Analysis

import os
import json
import requests
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dascore as dc
from matplotlib.colors import LogNorm


# Save-state caches: an .npz holds the amplitude matrices while a
# .json holds the magnitude-sorted labels. If both exist
# and contain the expected events, ProcessPoolExecutor block
# is skipped.
CACHE_NPZ  = "noise_analysis_cache.npz"
CACHE_JSON = "noise_analysis_cache.json"


def load_coords(filepath):
    """Loads .xycz file, and calculates cumulative distance."""
    df = pd.read_csv(filepath, sep=r'\s+', header=None,
        names=['lon', 'lat', 'cha', 'dep'],
        na_values=['', ' ', 'NaN', 'nan']).dropna(how='any')

    df['lon'] = df['lon'].apply(lambda x: x - 360 if x > 180 else x)

    R = 6371.0
    lat, lon = np.radians(df['lat'].values), np.radians(df['lon'].values)
    dlat, dlon = np.diff(lat), np.diff(lon)
    a = (np.sin(dlat / 2)**2 + np.cos(lat[:-1]) * np.cos(lat[1:]) *
         np.sin(dlon / 2)**2)

    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    df['dist_km'] = np.insert(R * c, 0, 0).cumsum()
    return df


def get_noise_rows(f_path, target_channels):
    """Extracts Median Absolute Amplitude of first 5s and last 20s."""
    if not os.path.exists(f_path):
        return None, None
    try:
        patch = dc.spool(f_path)[0].detrend("time")

        t_min = patch.coords.min('time')
        t_max = patch.coords.max('time')

        t_noise_end_5s = t_min + np.timedelta64(5, 's')
        patch_5s = patch.select(time=(t_min, t_noise_end_5s))

        t_noise_start_20s = t_max - np.timedelta64(20, 's')
        patch_20s = patch.select(time=(t_noise_start_20s, t_max))

        raw_noise_5s = np.squeeze(patch_5s.abs().median(dim="time").data)
        raw_noise_20s = np.squeeze(patch_20s.abs().median(dim="time").data)

        row_5s = np.full(len(target_channels), np.nan)
        row_20s = np.full(len(target_channels), np.nan)

        for i, cha_idx in enumerate(target_channels):
            cha_idx = int(cha_idx)
            if cha_idx >= 74 and 0 <= cha_idx < len(raw_noise_5s):
                row_5s[i] = raw_noise_5s[cha_idx]
                row_20s[i] = raw_noise_20s[cha_idx]

        return row_5s, row_20s
    except Exception as e:
        print(f"Error processing {os.path.basename(f_path)}: {e}")
        return None, None


def get_magnitude(event_id):
    """Fetches earthquake magnitude from USGS API using the event ID."""
    url = (f"https://earthquake.usgs.gov/fdsnws/event/1/query?"
           f"eventid={event_id}&format=geojson")
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            return data['properties']['mag']
    except Exception as e:
        print(f"Error fetching magnitude for {event_id}: {e}")
    return -1.0


def process_event(args):
    """Worker function to process a single event."""
    fname, noise_dir, target_channels = args
    event_id = fname.split('_')[0]
    f_path = os.path.join(noise_dir, fname)

    mag = get_magnitude(event_id)
    row_5s, row_20s = get_noise_rows(f_path, target_channels)

    return mag, event_id, row_5s, row_20s


def plot_and_save(terra_coords, data_matrix, valid_labels, title,
    output_path, vmin, vmax, avg_func=np.nanmean):
    """Generates the dual-axes plot (Bathymetry + Heatmap)."""
    avg_row = avg_func(data_matrix, axis=0)
    data_matrix_with_avg = [avg_row] + data_matrix
    labels_with_avg = ["Event Average"] + valid_labels

    Z = np.array(data_matrix_with_avg)
    Z = np.abs(Z)
    Z[Z <= 0] = 1e-10

    X_orig = terra_coords['dist_km'].values
    X_edges = np.zeros(len(X_orig) + 1)
    X_edges[1:-1] = (X_orig[:-1] + X_orig[1:]) / 2.0
    X_edges[0] = X_orig[0] - (X_orig[1] - X_orig[0]) / 2.0
    X_edges[-1] = X_orig[-1] + (X_orig[-1] - X_orig[-2]) / 2.0

    Y_edges = np.arange(len(labels_with_avg) + 1)

    fig, (ax_bath, ax_data) = plt.subplots(2, 1, figsize=(14, 10),
        gridspec_kw={'height_ratios': [1, 6]}, sharex=True)
    plt.subplots_adjust(hspace=0.02, bottom=0.1)

    ax_bath.plot(X_orig, -terra_coords['dep'], color='black', lw=1.5)
    ax_bath.axhline(0, color='blue', linestyle='--', lw=1)
    ax_bath.text(X_orig[0] + 5, -1, " Sea Level", color='blue',
        va='bottom', fontsize=7.5)
    ax_bath.invert_yaxis()
    ax_bath.set_ylabel("Depth (m)", fontsize=10, fontweight='bold')
    ax_bath.set_title(title, fontsize=14, pad=10)
    ax_bath.grid(True, alpha=0.15, linestyle=':')

    ax_data.pcolormesh(X_edges, Y_edges, Z, norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap='turbo', shading='flat')

    ax_data.set_yticks(np.arange(len(labels_with_avg)) + 0.5)
    ax_data.set_yticklabels(labels_with_avg, fontsize=7.5)

    yticklabels = ax_data.get_yticklabels()
    yticklabels[0].set_fontweight('bold')
    yticklabels[0].set_color('red')

    ax_data.invert_yaxis()
    ax_data.set_ylabel("Catalog Event ID", fontsize=11, fontweight='bold')
    ax_data.set_xlabel("Distance along cable (km)",
        fontsize=11, fontweight='bold', labelpad=5)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    print(f"Saving figure to: {output_path}")
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)


def _try_load_cache(cache_npz, cache_json, expected_events):
    """
    Attempts to load cached matrices and labels.

    Returns (matrix_5s, matrix_20s, valid_labels) if the cache exists,
    is readable, and contains all expected events. Returns None otherwise.
    """
    if not (os.path.exists(cache_npz) and os.path.exists(cache_json)):
        return None

    try:
        with open(cache_json) as f:
            meta = json.load(f)

        cached_ids = {entry['event_id'] for entry in meta}
        expected_ids = {fname.split('_')[0] for fname in expected_events}

        if not expected_ids.issubset(cached_ids):
            missing = expected_ids - cached_ids
            print(f"Stale Cache: Missing {len(missing)} event(s), "
                  f"reprocessing all...")
            return None

        c = np.load(cache_npz, allow_pickle=False)
        matrix_5s  = c['matrix_5s']
        matrix_20s = c['matrix_20s']
        valid_labels = [entry['label'] for entry in meta]

        print(f"[CACHE HIT] Loaded {len(valid_labels)} events from cache.")
        return matrix_5s, matrix_20s, valid_labels

    except Exception as e:
        print(f"Cache corrupted {e} reprocessing...")
        return None


if __name__ == "__main__":
    base_dir = "/Users/ed/research_code/das"
    noise_game_dir = os.path.join(base_dir, "das_records/good-events-3.2-up")
    coords_path = os.path.join(
        base_dir, 'das_coords_bathymetry/TERRA_coords.xycz'
    )
    figures_dir = os.path.join(base_dir, "das_figures")
    os.makedirs(figures_dir, exist_ok=True)

    cache_npz  = os.path.join(figures_dir, CACHE_NPZ)
    cache_json = os.path.join(figures_dir, CACHE_JSON)

    good_terra_events = [
        "ak0237eejw69_TERRA.h5", "ak023aw5mbdk_TERRA.h5",
        "ak023eccchzy_TERRA.h5", "ak0237q2shdo_TERRA.h5",
        "ak023bebgmhd_TERRA.h5", "ak023em715sv_TERRA.h5",
        "ak02381ibekf_TERRA.h5", "ak023bhlw02w_TERRA.h5",
        "ak023f7eyaqg_TERRA.h5", "ak0238ghnzxp_TERRA.h5",
        "ak023bkx215x_TERRA.h5", "ak023fhgggc6_TERRA.h5",
        "ak0238qkcxek_TERRA.h5", "ak023btef8mo_TERRA.h5",
        "ak023fnzshe1_TERRA.h5", "ak0239af45c3_TERRA.h5",
        "ak023bzqw7a7_TERRA.h5", "ak023frilkvn_TERRA.h5",
        "ak0239lyp68s_TERRA.h5", "ak023c3206y0_TERRA.h5",
        "ak023g35jxin_TERRA.h5", "ak0239qzbmym_TERRA.h5",
        "ak023cgc5fmi_TERRA.h5", "ak023gbgys2j_TERRA.h5",
        "ak0239saxy95_TERRA.h5", "ak023ctiuyia_TERRA.h5",
        "ak023gjh7z4b_TERRA.h5", "ak0239vxdtm6_TERRA.h5",
        "ak023d3dyqv0_TERRA.h5", "ak023godcr3i_TERRA.h5",
        "ak023a0j9eo0_TERRA.h5", "ak023gqcxl3z_TERRA.h5",
        "ak023a7ds9th_TERRA.h5", "ak023dif8i7c_TERRA.h5",
        "ak023a91bmgs_TERRA.h5", "ak023djxyhod_TERRA.h5",
        "us6000lggw_TERRA.h5", "ak023ah9zcn5_TERRA.h5",
        "ak023dk7iyjo_TERRA.h5", "ak023asybefb_TERRA.h5",
        "ak023ds94esw_TERRA.h5"
    ]

    terra_coords = load_coords(coords_path)
    target_channels = terra_coords['cha'].values
    
    # Caching to .json
    cached = _try_load_cache(cache_npz, cache_json, good_terra_events)

    if cached is not None:
        matrix_5s, matrix_20s, valid_labels = cached
    else:

        args_list = [(fname, noise_game_dir, target_channels)
                     for fname in good_terra_events]
        print(f"Parallel processing {len(args_list)} events...")

        with ProcessPoolExecutor() as ex:
            results = list(ex.map(process_event, args_list))

        extracted_data = [res for res in results
                          if res[2] is not None and res[3] is not None]

        if not extracted_data:
            print("Error: No data successfully extracted.")
            exit(1)

        extracted_data.sort(key=lambda x: x[0], reverse=True)

        matrix_5s    = []
        matrix_20s   = []
        valid_labels = []
        cache_meta   = []

        for mag, event_id, row_5s, row_20s in extracted_data:
            matrix_5s.append(row_5s)
            matrix_20s.append(row_20s)
            label = (f"M {mag:.1f} - {event_id}" if mag != -1.0
                     else f"M ? - {event_id}")
            valid_labels.append(label)
            cache_meta.append({'event_id': event_id, 'label': label})

        matrix_5s  = np.array(matrix_5s)
        matrix_20s = np.array(matrix_20s)

        # Caching to .npz
        np.savez_compressed(cache_npz,
                            matrix_5s=matrix_5s,
                            matrix_20s=matrix_20s)
        with open(cache_json, 'w') as f:
            json.dump(cache_meta, f, indent=4)
        print(f"Matrices saved to: {cache_npz}")
        print(f"Labels saved to: {cache_json}")


    median_5s       = np.nanmedian(matrix_5s, axis=0)
    norm_matrix_5s  = matrix_5s - median_5s
    median_20s      = np.nanmedian(matrix_20s, axis=0)
    norm_matrix_20s = matrix_20s - median_20s

    plot_and_save(terra_coords, matrix_5s.tolist(), valid_labels,
        "TERRA: Median Amplitude (First 5s Noise)",
        os.path.join(figures_dir, "TERRA_noise_floor_first_5.png"),
        vmin=1, vmax=1000, avg_func=np.nanmedian)

    plot_and_save(terra_coords, matrix_20s.tolist(), valid_labels,
        "TERRA: Median Amplitude (Last 20s Noise)",
        os.path.join(figures_dir, "TERRA_noise_floor_last_20.png"),
        vmin=1, vmax=1000, avg_func=np.nanmedian)

    plot_and_save(terra_coords, norm_matrix_5s.tolist(), valid_labels,
        "TERRA: Median Deviation (First 5s)",
        os.path.join(figures_dir, "TERRA_noise_deviation_first_5.png"),
        vmin=1, vmax=1000, avg_func=np.nanmean)

    plot_and_save(terra_coords, norm_matrix_20s.tolist(), valid_labels,
        "TERRA: Median Deviation (Last 20s)",
        os.path.join(figures_dir, "TERRA_noise_deviation_last_20.png"),
        vmin=1, vmax=1000, avg_func=np.nanmean)