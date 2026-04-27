# ETI2.0 Slide Maker
# Plots a single DAS channel trace (nearest to target_dist) from a TERRA
# event file, filtered above 1 Hz, and saves it as a PDF.  Used for a 
# presentation slide or two.

import os

import dascore as dc
import matplotlib.pyplot as plt
import numpy as np
from dascore.units import Hz


def main():
    target_dist = 17500
    file_path = "das_records/good-events-3.2-up/ak0239vxdtm6_TERRA.h5"
    output_file = "terra_single_channel_black.pdf"

    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Could not find {file_path}")

        patch = dc.spool(file_path)[0]
        processed_patch = patch.pass_filter(time=(1 * Hz, None))

        distances = processed_patch.coords.get_array('distance')
        channel_index = np.argmin(np.abs(distances - target_dist))
        actual_dist = distances[channel_index]
        print(f"Using channel {channel_index} at {actual_dist:.1f} m "
              f"(target: {target_dist} m)")

        trace_data = processed_patch.data[:, channel_index]
        times = processed_patch.coords.get_array('time')

        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(times, trace_data, color='black', linewidth=0.8)

        timebig = 19572.91452
        ax.set_xlabel("Time (hh-mm-ss)", fontsize=12)
        ax.set_ylabel("Strain (Radians)", fontsize=12)
        plt.xlim(timebig, timebig + 0.00145)
        plt.tight_layout()
        plt.savefig(output_file)
        plt.show()

        print(f"Successfully plotted channel {channel_index} to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == '__main__':
    main()
