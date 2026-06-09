#!/bin/env python3

import gzip
import argparse
import os

def get_arguments():
    """function to return the path to the sync file"""

    parser = argparse.ArgumentParser(description="Find homozygous vs heterozygous sites in a sync file")
    parser.add_argument("-s", "--sync", required=True, type=str)

    args = parser.parse_args()
    sync_file = args.sync
    assert os.path.isfile(sync_file), f"Could not location {sync_file}"

    return sync_file

def parse_sync_file(sync_file: str):
    """will create a results file to store sites of interest"""

    fh = gzip.open(sync_file, 'rt') if sync_file.endswith(".gz") else open(sync_file, 'r')
    outfh = open("Variants_to_consider.sync", 'w')

    # create key for easy look up of meta-data
    nucs_key = {i: nuc for i, nuc in enumerate("ATCGND")}
    samples_color = {}
    for i in range(4):
        if i < 2:
            samples_color[i] = "red"
        else:
            samples_color[i] = "yellow"

    for line in fh:
        fields = line.strip('\n').split('\t')
        allele_counts = [0, 0, 0 , 0]
        possible_error = False
        # start from first sample column
        for i, col in enumerate(fields[3:]):
            counts = col.split(":")
            for count in counts:
                count = int(count)
                if count > 0 :
                    allele_counts[i] += 1
                    if count < 5:
                        possible_error = True
        skip = False
        for i in range(4):
            if skip:
                break
            if allele_counts[i] >= 3:
                skip = True
                continue
        if skip:
            continue    
        if ((allele_counts[0] == 1 and allele_counts[1] == 1 and allele_counts[2] == 2 and allele_counts[3] == 2) 
            or (allele_counts[0] == 2 and allele_counts[1] == 2 and allele_counts[2] == 1 and allele_counts[3] == 1)):
            if possible_error:
                continue
            outfh.write(line)
            outfh.flush()

    fh.close()
    outfh.close()

def main():
    """entry point to the pipeline"""

    # get the file
    sync_file = get_arguments()

    # parse the file
    parse_sync_file(sync_file)


if __name__ == "__main__":
    main()