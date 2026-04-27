# Mapview TERRA Events
# Plots TERRA cable and earthquake epicenters on an Albers Equal Area
# Projection centered on lower Cook Inlet, with an inset Alaska map.  

import json
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def add_map_features(target_ax, scale='10m'):
    target_ax.add_feature(
        cfeature.OCEAN.with_scale(scale), facecolor='#ABC9DF', zorder=0)
    target_ax.add_feature(
        cfeature.LAND.with_scale(scale), facecolor='#F5F5DC',
        edgecolor='#666666', zorder=1)
    target_ax.add_feature(
        cfeature.LAKES.with_scale(scale), facecolor='#ABC9DF',
        edgecolor='#5A9CBF', alpha=0.9, zorder=2)
    target_ax.add_feature(
        cfeature.RIVERS.with_scale(scale), edgecolor='#5A9CBF',
        linewidth=0.7, zorder=2)
    target_ax.add_feature(
        cfeature.COASTLINE.with_scale(scale), linewidth=0.8,
        color='#2F4F4F', zorder=4)


def add_cities(target_ax):
    cities = {
        'Anchorage': (-149.9003, 61.2181),
        'Homer':     (-151.5483, 59.6425),
        'Seward':    (-149.4422, 60.1042),
        'Kenai':     (-151.2583, 60.5544)
    }
    for city, (lon, lat) in cities.items():
        extent = target_ax.get_extent(crs=ccrs.PlateCarree())
        if (extent[0] <= lon <= extent[1] and extent[2] <= lat <= extent[3]):
            target_ax.plot(lon, lat, 'ko', markersize=5,
                           transform=ccrs.PlateCarree(), zorder=30)
            target_ax.text(
                lon + 0.05, lat + 0.05, city,
                transform=ccrs.PlateCarree(), fontsize=10,
                fontweight='bold', zorder=31,
                bbox=dict(facecolor='white', alpha=0.7,
                          edgecolor='none', pad=1))


def albers_terra_plot():
    center_lat, center_lon = 60.0, -152.0
    lon_buffer, lat_buffer = 3.0, 1.5

    min_lat = center_lat - lat_buffer
    max_lat = center_lat + lat_buffer
    min_lon = center_lon - lon_buffer
    max_lon = center_lon + lon_buffer

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()

    json_path = os.path.join(script_dir, 'event_metadata_cache.json')
    data_dir = os.path.join(script_dir, '..', 'das_coords_bathymetry')

    try:
        terra_path = os.path.join(data_dir, 'TERRA_coords.xycz')
        terra = pd.read_csv(
            terra_path, sep=r'\s+', header=None,
            names=['lon', 'lat', 'cha', 'dep'])

        with open(json_path, 'r') as f:
            events_data = json.load(f)
            if isinstance(events_data, dict):
                events_df = pd.DataFrame(list(events_data.values()))
            else:
                events_df = pd.DataFrame(events_data)

    except Exception as e:
        print(f"Error loading files: {e}")
        return

    albers_proj = ccrs.AlbersEqualArea(
        central_longitude=center_lon, central_latitude=center_lat,
        standard_parallels=(55, 65))

    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(1, 1, 1, projection=albers_proj)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat],
                  crs=ccrs.PlateCarree())

    add_map_features(ax, scale='10m')
    add_cities(ax)

    ax.plot(terra['lon'], terra['lat'], color='#0000FF', linewidth=3.0,
            label='TERRA Cable', transform=ccrs.PlateCarree(), zorder=20)

    if not events_df.empty:
        lat_key = 'lat' if 'lat' in events_df.columns else 'latitude'
        lon_key = 'lon' if 'lon' in events_df.columns else 'longitude'

        events_df[lat_key] = pd.to_numeric(events_df[lat_key], errors='coerce')
        events_df[lon_key] = pd.to_numeric(events_df[lon_key], errors='coerce')
        events_df = events_df.dropna(subset=[lat_key, lon_key])

        events_subset = events_df[
            (events_df[lon_key] >= min_lon) & (events_df[lon_key] <= max_lon) &
            (events_df[lat_key] >= min_lat) & (events_df[lat_key] <= max_lat)]

        ax.scatter(
            events_subset[lon_key], events_subset[lat_key],
            s=180, color='#FF0000', edgecolor='black', linewidths=1.0,
            label='Earthquake Epicenters', transform=ccrs.PlateCarree(),
            zorder=25)

    ax_overview = inset_axes(
        ax, width="30%", height="30%", loc='upper left',
        axes_class=GeoAxes,
        axes_kwargs={'map_projection': ccrs.LambertConformal(
            central_longitude=-154, central_latitude=63)})

    ax_overview.set_extent([-170, -130, 54, 72], crs=ccrs.PlateCarree())
    add_map_features(ax_overview, scale='50m')

    n_pts = 50
    lons = np.r_[
        np.linspace(min_lon, max_lon, n_pts),
        [max_lon] * n_pts,
        np.linspace(max_lon, min_lon, n_pts),
        [min_lon] * n_pts]
    lats = np.r_[
        [min_lat] * n_pts,
        np.linspace(min_lat, max_lat, n_pts),
        [max_lat] * n_pts,
        np.linspace(max_lat, min_lat, n_pts)]

    ax_overview.plot(lons, lats, color='red', linewidth=2,
                     transform=ccrs.PlateCarree(), zorder=10)

    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.4,
                      zorder=5, x_inline=False, y_inline=False, color='gray')
    gl.top_labels = False
    gl.right_labels = False

    ax.legend(loc='lower right', frameon=True, facecolor='white',
              framealpha=0.9, fontsize=12)

    plt.tight_layout()
    output_file = os.path.join(script_dir, 'terra_albers_overview_map.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Map saved to {output_file}")
    plt.show()


if __name__ == "__main__":
    albers_terra_plot()
