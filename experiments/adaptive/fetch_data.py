"""Download NeuroKit2 sample recordings into data_external/ (skips existing files)."""

import os
import urllib.request

BASE = "https://raw.githubusercontent.com/neuropsychology/NeuroKit/master/data/"
FILES = ["bio_resting_5min_100hz.csv", "bio_resting_8min_100hz.csv", "bio_eventrelated_100hz.csv"]
DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_external")


def main():
    os.makedirs(DEST, exist_ok=True)
    for f in FILES:
        path = os.path.join(DEST, f)
        if os.path.exists(path):
            print(f"exists: {f}")
            continue
        print(f"downloading {f} ...")
        urllib.request.urlretrieve(BASE + f, path)


if __name__ == "__main__":
    main()
