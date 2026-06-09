#!/bin/env python3

import gzip
import argparse
import os

def get_arguments() -> str:
    """function to return the path to the sync file"""

    parser = argparse.ArgumentParser(description="Get Fst stats from sync file")
    parser.add_argument("-s", "--sync", required=True, type=str)

    args = parser.parse_args()
    sync_file = args.sync
    assert os.path.isfile(sync_file), f"Could not location {sync_file}"

    return sync_file

def parse_sync_file(sync_file: str) -> int:
    """get the stats & print to the console"""

    fh    = gzip.open(sync_file, 'rt') if sync_file.endswith(".gz") else open(sync_file, 'r')
    cnts  = [0]   * 4
    sums  = [0.0] * 4
    total = 0

    # columns
    r1_vs_y1 = 6
    r1_vs_y2 = 7
    r2_vs_y1 = 8
    r2_vs_y2 = 9

    for line in fh:
        fields = line.strip('\n').split('\t')
        fst1   = float(fields[r1_vs_y1].split('=')[1])
        fst2   = float(fields[r1_vs_y2].split('=')[1])
        fst3   = float(fields[r2_vs_y1].split('=')[1])
        fst4   = float(fields[r2_vs_y2].split('=')[1])
        if (fst1 > 0.0):
            cnts[0] += 1
            sums[0] += fst1
        if (fst2 > 0.0):
            cnts[1] += 1
            sums[1] += fst2
        if (fst3 > 0.0):
            cnts[2] += 1
            sums[2] += fst3
        if (fst4 > 0.0):
            cnts[3] += 1
            sums[3] += fst4
        if (fst1 > 0.0 or fst2 > 0.0 or fst3 > 0.0 or fst4 > 0.0):
            total += 1        

    fh.close()

    messages = ["Red 1 vs Yellow 1", "Red 1 vs Yellow 2", "Red 2 vs Yellow 1", "Red 2 vs Yellow 2"]

    print(f"Found a total of {total} variant positions")

    for i in range(4):
        msg = messages[i]
        cnt = cnts[i]
        tot = sums[i]
        avg = tot / cnt
        msg = f"{msg}: Total sites {cnt}, avg fst {avg}"
        print(msg)

    return 0

def main():
    """entry point to the pipeline"""

    # get the file
    sync_file = get_arguments()

    # parse the file
    parse_sync_file(sync_file)


if __name__ == "__main__":
    main()