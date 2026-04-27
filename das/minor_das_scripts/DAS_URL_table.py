# DAS URL Table
# Queries the USGS FDSN event service and prints a formatted table of
# earthquakes in the Cook Inlet / Kenai Peninsula study area with
# Magnitude ≥ 3.2, from Jun–Dec 2023). Just a 'checking' script.

from obspy import read_events


def main():
    url = (
        "https://earthquake.usgs.gov/fdsnws/event/1/query?"
        "format=quakeml&starttime=2023-06-10&endtime=2023-12-31&"
        "latitude=59.6&longitude=-152.2&maxradiuskm=150&minmagnitude=3.2")

    catalog = read_events(url)
    catalog.events.sort(key=lambda x: x.origins[0].time)

    print(f"{'Origin Time':<25} | {'Mag':<5} | {'Lat':<8} | "
          f"{'Lon':<8} | {'Depth (km)'}")
    print("-" * 70)

    for event in catalog:
        pref_origin = event.preferred_origin() or event.origins[0]
        pref_mag = event.preferred_magnitude() or event.magnitudes[0]

        otime = pref_origin.time
        lat = pref_origin.latitude
        lon = pref_origin.longitude
        depth = pref_origin.depth / 1000.0
        mag = pref_mag.mag

        print(f"{str(otime):<25} | {mag:<5.1f} | {lat:<8.3f} | "
              f"{lon:<8.3f} | {depth:.1f}")


if __name__ == '__main__':
    main()
