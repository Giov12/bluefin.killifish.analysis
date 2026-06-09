#!/bin/env python3

import gzip
import argparse
import os
import sys

def get_arguments():
    """function to return the path to the sync file"""

    parser = argparse.ArgumentParser(description="Find homozygous sites flanking a target region in a .sync file")
    parser.add_argument("-s", "--sync", required=True, type=str)
    parser.add_argument("-r", "--rightbound", type = int, default=10000, help = "right most position of target region")
    parser.add_argument("-l", "--leftbound", type = int, default=0, help = "left most position of target region")
    parser.add_argument("-f", "--flanking_distance", type = int, default=5000, help = "distance in base pairs flanking the target bounds")
    parser.add_argument("-m", "--min_count", type = int, default=5, help = "minimum number of read counts to consider a position, else 0")
    parser.add_argument("-c", "--chrom", type = str, required=True, help = "name of chromosome that the target is located at")
    parser.add_argument("-e", "--error_aware", action="store_true", help = "be more strict and don't trust sites with possible errors")

    args = parser.parse_args()
    sync_file = args.sync
    assert os.path.isfile(sync_file), f"Could not location {sync_file}"
    rightbound = args.rightbound
    leftbound = args.leftbound
    flanking_dist = args.flanking_distance
    min_count = args.min_count
    for i in [rightbound, leftbound, flanking_dist, min_count]:
        assert i >= 0, f"{i} is below 0 base pairs"
    assert leftbound < rightbound, "leftbound cannot be lower than rightbound"
    chrom = args.chrom
    error_aware = args.error_aware

    return sync_file, leftbound, rightbound, flanking_dist, min_count, chrom, error_aware

def make_output_name(sync_file: str):
    """creates the output name"""

    bname = os.path.basename(sync_file)
    if sync_file.endswith(".gz"):
        extension = ".sync.gz"
    else:
        extension = ".sync"

    outname = bname.replace(extension, "_homozygous_sites.tsv")

    return outname

def parse_alleles(allele_column: str, min_count: int, error_aware: bool):
    """parse the allele count column in a sync file"""

    alleles = []
    allele_key = ['A', 'T', 'C', 'G', 'N', 'D']

    allele_counts = allele_column.split(":")

    possible_error = False

    for i, allele_cnt in enumerate(allele_counts):
        allele_cnt = int(allele_cnt)
        if error_aware:
            if allele_cnt > 0:
                alleles.append(allele_key[i])
                if allele_cnt < min_count:
                    possible_error = True
        elif allele_cnt >= min_count:
            alleles.append(allele_key[i])

    return alleles, possible_error


def parse_sync_file(sync_file: str, leftbound: int, rightbound: int, flanking_dist: int, min_count: int, chrom: str, error_aware: bool):
    """will create a results file to store sites of interest"""

    fh = gzip.open(sync_file, 'rt') if sync_file.endswith(".gz") else open(sync_file, 'r')
    
    # do the calculation just once
    left_limit = leftbound - flanking_dist
    right_limit = rightbound + flanking_dist
    positions = [] # will store the positions to create the window

    for line in fh:
        fields = line.strip('\n').split('\t')
        if fields[0] == chrom:
            pos = int(fields[1])
            if (left_limit <= pos < leftbound) or (rightbound < pos <= right_limit):
                alleles = []
                if error_aware:
                    errors = []
                for allele_column in fields[3:]:
                    sample_alleles, possible_error = parse_alleles(allele_column, min_count, error_aware)
                    alleles.append(sample_alleles)
                    if error_aware:
                        errors.append(possible_error)
                homozygous = (alleles[0] == alleles[1] == alleles[2] == alleles[3])
                print(alleles)
                if (error_aware == False) and (homozygous == True):
                    positions.append(pos)
                elif (error_aware == True) and (homozygous == True):
                    trust = True
                    for e in errors:
                        if e == False:
                            trust = False
                            break
                    if trust:
                        positions.append(pos)
    fh.close()

    if len(positions) == 0:
        msg = "Did not find any sites to determine homozygous windows"
        sys.exit(msg)

    # now to make the windows
    outfh = open(make_output_name(sync_file), 'w')

    prev_pos = positions[0]
    start_pos = prev_pos # this will be the start

    # if a possible window can be form
    if len(positions) > 1:
        for p in positions[1:]:
            if p == prev_pos + 1:
                prev_pos = p # extend
            else:
                wind_size = 1 if start_pos == prev_pos else prev_pos - start_pos + 1
                outline = f"{chrom}\t{start_pos}\t{prev_pos}\t{wind_size}\n"
                outfh.write(outline)
                start_pos = p
                prev_pos = p
        # get last bit
        if start_pos == prev_pos:
            wind_size = 1
        else:
            wind_size = prev_pos - start_pos
        outline = f"{chrom}\t{start_pos}\t{prev_pos}\t{wind_size}\n"
        outfh.write(outline)
    # only 1 site found
    else:
        outline = f"{chrom}\t{start_pos}\t{prev_pos}\t{1}\n"
        outfh.write(outline)

    
    outfh.close()

def main():
    """entry point to the pipeline"""

    # get the inputs
    sync_file, leftbound, rightbound, \
        min_count, flanking_dist, chrom, error_aware = get_arguments()

    # parse the file and create the windows
    parse_sync_file(sync_file, leftbound, rightbound, flanking_dist,
                    min_count, chrom, error_aware)


if __name__ == "__main__":
    main()