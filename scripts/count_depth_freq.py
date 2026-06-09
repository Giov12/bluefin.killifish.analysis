#!/bin/env python3

import argparse
import gzip
import os
from   collections import defaultdict

def get_arguments():
    """get the arguments"""

    parser = argparse.ArgumentParser(
        description="Count the frequency of each depth from a .depth file from Samtools")
    parser.add_argument(
        "-d", "--depth", help="Depth file", required=True)
    parser.add_argument(
        "-o", "--outname", help="Name of output file [optional]", required=False, default='')

    args       = parser.parse_args()
    depth_file = args.depth
    outname    = args.outname

    assert os.path.isfile(depth_file), f"Could not locate {depth_file}"

    outname = outname if (outname != '') else make_outname(depth_file)

    return depth_file, outname


def make_outname(depth_file: str):
    """create an output name"""
    
    bname = os.path.basename(depth_file)
    ext   = ".depth.gz" if bname.endswith(".gz") else ".depth"
    oname = bname.replace(ext, "_depth_counts.tsv")
    return oname 

def get_freq(depth_file: str, outname: str):
    """start going through the positions"""

    fh  = gzip.open(depth_file, "rt") if (depth_file.endswith(".gz")) else open(depth_file, 'r')
    freq = defaultdict(int)

    for line in fh:
        fields = line.strip('\n').split('\t')
        freq[fields[2]] += 1

    fh.close()

    ofh = open(outname, 'w')

    for cnt, frq in freq.items():
        ofh.write(f"{cnt}\t{frq}\n")

    ofh.close()

def main():
    """Begin the pipeline here"""

    # chrom pos depth
    depth_file, outname = get_arguments()

    get_freq(depth_file, outname)

if __name__ == "__main__":
    main()
