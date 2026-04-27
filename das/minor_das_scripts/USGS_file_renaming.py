# DAS USGS File Renaming
# Renames locally downloaded DAS event files from the old date-based 
# naming convention to USGS event ID-based names.

import os

# True = preview only; False = apply changes.
DRY_RUN = True

RECORDS_PATH = os.path.join(os.getcwd(), "das_records", "good-events-3.2-up")

NAME_MAP = {
    "2023-06-10-08-44-56.996000ML3.3TERRA": "ak0237eejw69",
    "2023-06-17-14-28-52.384000ML3.4TERRA": "ak0237q2shdo",
    "2023-06-24-05-49-01.125000ML3.4TERRA": "ak02381ibekf",
    "2023-07-03-11-39-07.216000ML3.6TERRA": "ak0238ghnzxp",
    "2023-07-09-20-57-48.377000ML3.9TERRA": "ak0238qkcxek",
    "2023-07-21-21-09-05.673000ML3.5TERRA": "ak0239af45c3",
    "2023-07-28-19-11-22.434000ML5.0TERRA": "ak0239lyp68s",
    "2023-07-31-22-49-11.705000ML3.5TERRA": "ak0239qzbmym",
    "2023-08-01-02-47-42.411000ML3.5TERRA": "ak0239saxy95",
    "2023-08-03-21-56-48.577000ML3.9TERRA": "ak0239vxdtm6",
    "2023-08-06-00-17-51.206000ML4.7TERRA": "ak023a0j9eo0",
    "2023-08-10-14-24-07.870000ML3.5TERRA": "ak023a7ds9th",
    "2023-08-11-14-24-09.686000ML3.2TERRA": "ak023a91bmgs",
    "2023-08-16-12-51-42.839000ML3.7TERRA": "ak023ah9zcn5",
    "2023-08-23-18-51-52.010000ML3.7TERRA": "ak023asybefb",
    "2023-08-25-12-18-03.189000ML3.5TERRA": "ak023aw5mbdk",
    "2023-09-05-10-38-35.907000ML3.7TERRA": "ak023bebgmhd",
    "2023-09-07-09-29-58.254000ML3.4TERRA": "ak023bhlw02w",
    "2023-09-09-09-45-14.634000ML3.2TERRA": "ak023bkx215x",
    "2023-09-14-22-35-05.805000ML4.4TERRA": "ak023btef8mo",
    "2023-09-18-06-09-43.358000ML3.4TERRA": "ak023bzqw7a7",
    "2023-09-20-06-15-49.537000ML3.2TERRA": "ak023c3206y0",
    "2023-09-28-09-25-59.004000ML3.4TERRA": "ak023cgc5fmi",
    "2023-10-06-06-58-47.455000ML3.8TERRA": "ak023ctiuyia",
    "2023-10-12-03-02-49.457000ML3.4TERRA": "ak023d3dyqv0",
    "2023-10-12-03-02-50.850000ML3.3TERRA": "ak023d3dyut1",
    "2023-10-18-22-02-33.464000ML3.5TERRA": "us6000lggw",
    "2023-10-21-12-15-27.172000ML3.5TERRA": "ak023dif8i7c",
    "2023-10-22-04-06-43.215000ML3.7TERRA": "ak023djxyhod",
    "2023-10-22-20-14-00.011000ML3.8TERRA": "ak023dk7iyjo",
    "2023-10-27-06-54-57.303000ML3.3TERRA": "ak023ds94esw",
    "2023-11-08-21-28-55.242000ML3.5TERRA": "ak023eccchzy",
    "2023-11-14-17-02-01.362000ML3.4TERRA": "ak023em715sv",
    "2023-11-27-00-11-32.040000ML3.8TERRA": "ak023f7eyaqg",
    "2023-12-03-07-30-16.367000ML3.3TERRA": "ak023fhgggc6",
    "2023-12-07-02-57-27.077000ML4.3TERRA": "ak023fnzshe1",
    "2023-12-09-15-56-18.535000ML3.3TERRA": "ak023frilkvn",
    "2023-12-16-19-24-40.053000ML4.7TERRA": "ak023g35jxin",
    "2023-12-21-22-14-19.497000ML3.4TERRA": "ak023gbgys2j",
    "2023-12-26-06-30-02.910000ML3.2TERRA": "ak023gjh7z4b",
    "2023-12-29-02-13-56.903000ML3.5TERRA": "ak023godcr3i",
    "2023-12-30-22-36-05.065000ML3.4TERRA": "ak023gqcxl3z"
}


def rename_das_files(directory, mapping):
    print(f"Checking directory: {directory}\n")

    count = 0
    for old_name, usgs_id in mapping.items():
        old_path = os.path.join(directory, old_name)
        new_name = f"{usgs_id}_TERRA.h5"
        new_path = os.path.join(directory, new_name)

        if os.path.exists(old_path):
            if DRY_RUN:
                print(f"[DRY RUN] Would rename: {old_name} -> {new_name}")
            else:
                os.rename(old_path, new_path)
                print(f"Renamed: {old_name} -> {new_name}")
            count += 1
        elif os.path.exists(new_path):
            print(f"[INFO] Already renamed: {new_name}")
        else:
            print(f"[SKIP] File not found: {old_name}")

    if DRY_RUN:
        print(f"\nDry run complete. Found {count} files to rename.")
        print("Set DRY_RUN = False to apply changes.")
    else:
        print(f"\nSuccess! Renamed {count} files.")


if __name__ == "__main__":
    rename_das_files(RECORDS_PATH, NAME_MAP)
