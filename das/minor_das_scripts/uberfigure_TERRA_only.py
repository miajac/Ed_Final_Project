# uberfigure_TERRA_only.py
# Plots a single TERRA DAS event as a reference wiggle plot with zoomed
# insets, highlighting a target channel (default: 17 500 m) in red.

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
EVENT_FILE  = "ak0239vxdtm6_TERRA.h5"
OUTPUT_FILE = "uberfigure_terra_highlighted.pdf"
# ─────────────────────────────────────────────────────────────────────


def matplotlib_wiggle(
    ax, data, times, dist_coords, scale, linewidth, alpha,
    highlight_dist=None,
):
    """
    Wiggle plot for a DAS patch.  Optionally highlights one trace in red.
    """
    if data.size == 0:
        return

    n_samples, n_channels = data.shape

    highlight_idx = None
    if highlight_dist is not None:
        highlight_idx = int(np.argmin(np.abs(dist_coords - highlight_dist)))

    global_max = np.nanpercentile(np.abs(data), 99)
    if global_max <= 0 or not np.isfinite(global_max):
        global_max = 1.0

    for i in range(n_channels):
        trace = data[:, i]
        if np.all(np.isnan(trace)):
            continue

        scaled = (np.nan_to_num(trace) / global_max) * scale

        if i == highlight_idx:
            color, lw, alpha_ = "red", linewidth * 1.5, 1.0
        else:
            color, lw, alpha_ = "black", linewidth, alpha

        ax.plot(times, i + scaled, color=color, linewidth=lw, alpha=alpha_)

    ax.set_ylim(-1, n_channels)


def plot_zoomed_wiggle(patch, ax, title, highlight_dist=17_500):
    """
    High-resolution wiggle plot for a zoomed DAS sub-patch with a red
    highlight at *highlight_dist* metres.
    """
    try:
        processed = (
            patch.dropna(dim="time", how="all")
                 .pass_filter(time=(1 * Hz, None))
        )
        data  = processed.data
        times = processed.coords.get_array("time")
        dists = processed.coords.get_array("distance")

        matplotlib_wiggle(
            ax, data, times, dists,
            scale=2, linewidth=3, alpha=0.5,
            highlight_dist=highlight_dist,
        )
        ax.set_title(title, fontsize=10)
        ax.invert_yaxis()

    except Exception as exc:
        print(f"Failed to plot wiggle for {title}: {exc}")
        ax.set_title(f"Failed to plot: {title}")


def set_meter_ticks(ax, dist_coords, tick_interval_m=1_000):
    """Label the y-axis in metres using the supplied distance array."""
    try:
        min_d = np.min(dist_coords)
        max_d = np.max(dist_coords)
        start = np.ceil(min_d / tick_interval_m) * tick_interval_m
        end   = np.floor(max_d / tick_interval_m) * tick_interval_m
        labels  = np.arange(start, end + tick_interval_m, tick_interval_m)
        indices = [np.argmin(np.abs(dist_coords - m)) for m in labels]
        ax.set_yticks(indices)
        ax.set_yticklabels([f"{int(m)}" for m in labels])
        ax.set_ylabel("Distance (m)", fontsize=12)
    except Exception:
        ax.set_ylabel("Channel Index")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    terra_path  = RECORDS_DIR / EVENT_FILE
    output_path = FIGURES_DIR / OUTPUT_FILE

    if not terra_path.exists():
        print(f"ERROR: File not found at {terra_path}")
        return

    try:
        terra_patch = dc.spool(terra_path)[0]
    except Exception as exc:
        print(f"Error loading file: {exc}")
        return

    # ── Zoom region definitions ───────────────────────────────────────
    # Each entry: 'dist' = (min_m, max_m); other keys = (start_s, end_s)
    terra_subsets = {
        "Dist_17.5k_Center": {
            "dist":      (17_400, 17_600),
            "P-Arrival": (9, 14),
        }
    }
    # ─────────────────────────────────────────────────────────────────

    ref_h  = 40
    zoom_h = 40
    n_zooms = sum(len(v) - 1 for v in terra_subsets.values())
    n_plots = 1 + n_zooms
    h_ratios = [ref_h] + [zoom_h] * n_zooms

    fig, axes = plt.subplots(
        n_plots, 1,
        figsize=(20, sum(h_ratios)),
        gridspec_kw={"height_ratios": h_ratios},
    )
    if n_plots == 1:
        axes = [axes]

    ax_ref  = axes[0]
    t_start = terra_patch.coords.min("time")

    # Reference plot (decimated)
    ref = (
        terra_patch
        .decimate(distance=4)
        .pass_filter(time=(1 * Hz, None))
    )
    ref_dists = ref.coords.get_array("distance")

    matplotlib_wiggle(
        ax_ref, ref.data, ref.coords.get_array("time"),
        ref_dists, scale=5, linewidth=0.5, alpha=0.3,
    )
    ax_ref.set_title(
        f"TERRA Reference Plot: {terra_path.name}", fontsize=26,
    )
    set_meter_ticks(ax_ref, ref_dists, 2_000)
    ax_ref.invert_yaxis()

    # Zoom plots + boxes
    plot_idx = 1
    for subset_info in terra_subsets.values():
        dist_range = subset_info["dist"]
        for key, (start_s, end_s) in subset_info.items():
            if key == "dist":
                continue

            t0 = t_start + np.timedelta64(start_s, "s")
            t1 = t_start + np.timedelta64(end_s, "s")

            y0 = int(np.argmin(np.abs(ref_dists - dist_range[0])))
            y1 = int(np.argmin(np.abs(ref_dists - dist_range[1])))
            rect = patches.Rectangle(
                (t0, y0), t1 - t0, y1 - y0,
                linewidth=4, edgecolor="red",
                facecolor="none", linestyle="--", zorder=10,
            )
            ax_ref.add_patch(rect)

            ax_zoom = axes[plot_idx]
            zoomed  = terra_patch.select(
                time=(t0, t1), distance=dist_range,
            )
            title = (
                f"Zoom (TERRA): {key}\n"
                f"Distance: {dist_range[0]}m–{dist_range[1]}m | "
                f"Time: {start_s}s–{end_s}s"
            )
            plot_zoomed_wiggle(zoomed, ax_zoom, title)
            plot_idx += 1

    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()
