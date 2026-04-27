# uberfigure_combined.py
# Plots both KKFLS (negative distances) and TERRA (positive distances)
# on a shared reference wiggle plot, with labelled zoom insets for each
# cable pulled from the subset dictionaries below.

import os
from pathlib import Path

import dascore as dc
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from dascore.units import Hz

HERE = Path(__file__).parent

# ── Configurable paths ────────────────────────────────────────────────
RECORDS_DIR = HERE / "das_records" / "good-events-3.2-up"
FIGURES_DIR = HERE / "das_figures"
# Update these to point at the two HDF5 files for the event you want.
KKFLS_FILE  = "2023-08-03KKFLS.h5"
TERRA_FILE  = "2023-08-03TERRA.h5"
OUTPUT_FILE = "uberfigure_combined.pdf"
# ─────────────────────────────────────────────────────────────────────


def matplotlib_wiggle(ax, data, times, scale, linewidth, alpha):
    """
    Wiggle plot for a DAS patch.  Each channel is drawn as a waveform
    centred on its channel index, scaled by the 99th-percentile amplitude.
    """
    if data.size == 0:
        print("matplotlib_wiggle: empty data, skipping.")
        return

    n_samples_from_coords = len(times)
    if data.shape[0] != n_samples_from_coords:
        print(
            f"Shape mismatch: data {data.shape} vs "
            f"{n_samples_from_coords} time samples."
        )
        return

    _n_samples, n_channels = data.shape

    global_max = np.nanpercentile(np.abs(data), 99)
    if global_max == 0 or not np.isfinite(global_max):
        global_max = 1.0

    for i in range(n_channels):
        trace = data[:, i]
        if np.all(np.isnan(trace)):
            continue
        scaled = (np.nan_to_num(trace) / global_max) * scale
        ax.plot(times, i + scaled, color="black",
                linewidth=linewidth, alpha=alpha)

    ax.set_ylim(-1, n_channels)


def plot_zoomed_wiggle(patch, ax, title):
    """High-resolution wiggle plot for a zoomed DAS sub-patch."""
    try:
        processed = (
            patch.dropna(dim="time", how="all")
                 .pass_filter(time=(1 * Hz, None))
        )
        matplotlib_wiggle(
            ax, processed.data,
            processed.coords.get_array("time"),
            scale=2, linewidth=1, alpha=0.5,
        )
        ax.set_title(title, fontsize=10)
        ax.invert_yaxis()

    except Exception as exc:
        print(f"Failed to plot wiggle for {title}: {exc}")
        ax.set_title(f"Failed to plot: {title}")


def draw_zoom_box(ax_ref, ref_dist_coords, time_range,
                  distance_range, label_text):
    """
    Draw a labelled dashed rectangle on *ax_ref* that marks the region
    expanded in the corresponding zoom panel.
    """
    i0 = int(np.argmin(np.abs(ref_dist_coords - distance_range[0])))
    i1 = int(np.argmin(np.abs(ref_dist_coords - distance_range[1])))
    y1, y2 = min(i0, i1), max(i0, i1)

    t0, t1 = time_range
    rect = patches.Rectangle(
        (t0, y1), t1 - t0, y2 - y1,
        linewidth=1.5, edgecolor="red",
        facecolor="none", linestyle="--",
    )
    ax_ref.add_patch(rect)

    dt = t1 - t0
    ax_ref.text(
        t0 - dt / 2, y1 + (y2 - y1) / 2, label_text,
        color="red", ha="center", va="center",
        fontsize=12, fontweight="bold",
        bbox=dict(
            facecolor="white", alpha=0.9,
            pad=0.1, boxstyle="square",
        ),
    )


def set_meter_ticks(ax, dist_coords, tick_interval_m=1_000):
    """Label the y-axis in metres using the supplied distance array."""
    try:
        min_d = np.min(dist_coords)
        max_d = np.max(dist_coords)
        start = np.ceil(min_d / tick_interval_m) * tick_interval_m
        end   = np.floor(max_d / tick_interval_m) * tick_interval_m

        if end < start:
            ax.set_ylabel("Distance (m)")
            return

        labels  = np.arange(start, end + tick_interval_m, tick_interval_m)
        indices = []
        for m in labels:
            indices.append(int(np.argmin(np.abs(dist_coords - m))))

        ax.set_yticks(indices)
        ax.set_yticklabels([f"{int(m)}" for m in labels])
        ax.set_ylabel("Distance (m)", fontsize=12)

    except Exception as exc:
        print(f"Failed to set metre ticks: {exc}")
        ax.set_ylabel("Channel Index")


def combine_patches(kkfls_patch, terra_patch):
    """
    Merge KKFLS (distances negated) and TERRA patches into a single
    array sorted by distance.

    Returns
    -------
    sorted_data : ndarray, shape (n_times, n_channels)
    times       : ndarray
    sorted_dist : ndarray
    """
    data_k = kkfls_patch.data
    dist_k = kkfls_patch.coords.get_array("distance") * -1
    data_t = terra_patch.data
    dist_t = terra_patch.coords.get_array("distance")
    times  = kkfls_patch.coords.get_array("time")

    # Drop duplicate zero channel from KKFLS
    zero_k = np.where(dist_k == 0)[0]
    zero_t = np.where(dist_t == 0)[0]
    if zero_k.size > 0 and zero_t.size > 0:
        print("Removing duplicate zero channel from KKFLS.")
        data_k = np.delete(data_k, zero_k, axis=1)
        dist_k = np.delete(dist_k, zero_k)

    combined_data = np.concatenate((data_k, data_t), axis=1)
    combined_dist = np.concatenate((dist_k, dist_t))

    order       = np.argsort(combined_dist)
    sorted_data = combined_data[:, order]
    sorted_dist = combined_dist[order]

    print("Combined patch created.")
    return sorted_data, times, sorted_dist


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    kkfls_path  = RECORDS_DIR / KKFLS_FILE
    terra_path  = RECORDS_DIR / TERRA_FILE
    output_path = FIGURES_DIR / OUTPUT_FILE

    try:
        kkfls_patch = dc.spool(kkfls_path)[0]
        terra_patch = dc.spool(terra_path)[0]
        print("Files loaded successfully.")
    except Exception as exc:
        print(f"Could not load files: {exc}")
        return

    # ── Zoom region definitions ───────────────────────────────────────
    # 'dist' = (min_m, max_m); other keys = (start_s, end_s).
    # KKFLS distances are positive here; they are negated internally.
    kkfls_subsets = {
        "Dist_2k_3k": {
            "dist":      (2_000, 3_000),
            "P-Arrival": (10, 14),
            "S-Arrival": (24, 28),
        },
        "Dist_30k_31k": {
            "dist":      (30_000, 31_000),
            "P-Arrival": (8, 12),
            "S-Arrival": (22, 26),
        },
        "Dist_34k_35k": {
            "dist":      (34_000, 35_000),
            "P-Arrival": (8, 12),
            "S-Arrival": (22, 26),
        },
        "Dist_45k_46k": {
            "dist":      (45_000, 46_000),
            "P-Arrival": (8, 12),
            "S-Arrival": (21, 25),
        },
    }

    terra_subsets = {
        "Dist_8k_9k": {
            "dist":      (8_000, 9_000),
            "P-Arrival": (10, 14),
            "S-Arrival": (23, 27),
        },
        "Dist_22k_23k": {
            "dist":      (22_000, 23_000),
            "P-Arrival": (9, 13),
            "S-Arrival": (22, 26),
        },
        "Dist_30k_31k": {
            "dist":      (30_000, 31_000),
            "P-Arrival": (9, 13),
            "S-Arrival": (21, 25),
        },
        "Dist_42k_43k": {
            "dist":      (42_000, 43_000),
            "P-Arrival": (7, 11),
            "S-Arrival": (20, 24),
        },
        "Dist_53k_54k": {
            "dist":      (53_000, 54_000),
            "P-Arrival": (6, 10),
            "S-Arrival": (19, 23),
        },
    }
    # ─────────────────────────────────────────────────────────────────

    ref_h  = 40
    zoom_h = 30

    n_kkfls = sum(len(v) - 1 for v in kkfls_subsets.values())
    n_terra = sum(len(v) - 1 for v in terra_subsets.values())
    n_zooms = n_terra + n_kkfls
    n_plots = 1 + n_zooms

    total_h  = ref_h + n_zooms * zoom_h
    h_ratios = [ref_h] + [zoom_h] * n_zooms

    fig, axes = plt.subplots(
        n_plots, 1,
        figsize=(20, total_h),
        gridspec_kw={"height_ratios": h_ratios},
    )
    if n_plots == 1:
        axes = [axes]

    try:
        ax_ref = axes[0]

        t_start_k = kkfls_patch.coords.min("time")
        t_start_t = terra_patch.coords.min("time")

        dec_k = (
            kkfls_patch
            .decimate(distance=4, filter_type="iir")
            .pass_filter(time=(1 * Hz, None))
        )
        dec_t = (
            terra_patch
            .decimate(distance=4, filter_type="iir")
            .pass_filter(time=(1 * Hz, None))
        )

        ref_data, ref_times, ref_dists = combine_patches(dec_k, dec_t)

        matplotlib_wiggle(
            ax_ref, ref_data, ref_times,
            scale=5, linewidth=0.5, alpha=0.3,
        )

        ax_ref.set_title(
            "Combined Reference: KKFLS (Negative) & TERRA (Positive)",
            fontsize=26,
        )
        set_meter_ticks(ax_ref, ref_dists, 2_000)
        ax_ref.invert_yaxis()

        plot_idx    = 1
        box_counter = 1

        # ── TERRA zoom panels ─────────────────────────────────────────
        for subset_name, subset_info in terra_subsets.items():
            dist_range = subset_info["dist"]
            for key, (start_s, end_s) in subset_info.items():
                if key == "dist":
                    continue

                t0 = t_start_t + np.timedelta64(start_s, "s")
                t1 = t_start_t + np.timedelta64(end_s, "s")

                ax_zoom = axes[plot_idx]
                zoomed  = terra_patch.select(
                    time=(t0, t1), distance=dist_range,
                )
                title = (
                    f"Zoom {box_counter} (TERRA): {key} — "
                    f"{subset_name}\n"
                    f"Distance: {dist_range[0]}m–{dist_range[1]}m | "
                    f"Time: {start_s}s–{end_s}s"
                )
                plot_zoomed_wiggle(zoomed, ax_zoom, title)
                draw_zoom_box(
                    ax_ref, ref_dists, (t0, t1),
                    dist_range, str(box_counter),
                )

                plot_idx    += 1
                box_counter += 1

        # ── KKFLS zoom panels ─────────────────────────────────────────
        for subset_name, subset_info in kkfls_subsets.items():
            dist_range = subset_info["dist"]
            for key, (start_s, end_s) in subset_info.items():
                if key == "dist":
                    continue

                t0 = t_start_k + np.timedelta64(start_s, "s")
                t1 = t_start_k + np.timedelta64(end_s, "s")

                ax_zoom = axes[plot_idx]
                zoomed  = kkfls_patch.select(
                    time=(t0, t1), distance=dist_range,
                )
                title = (
                    f"Zoom {box_counter} (KKFLS): {key} — "
                    f"{subset_name}\n"
                    f"Distance: {dist_range[0]}m–{dist_range[1]}m | "
                    f"Time: {start_s}s–{end_s}s"
                )
                plot_zoomed_wiggle(zoomed, ax_zoom, title)

                # KKFLS distances are negative in the reference frame
                neg_range = (-dist_range[0], -dist_range[1])
                draw_zoom_box(
                    ax_ref, ref_dists, (t0, t1),
                    neg_range, str(box_counter),
                )

                plot_idx    += 1
                box_counter += 1

        print("\nAll panels plotted.")
        plt.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved to {output_path}")

    except Exception as exc:
        print(f"Failed to build figure: {exc}")
        if "fig" in dir():
            plt.close(fig)


if __name__ == "__main__":
    main()
