#!/bin/python

import argparse
import os
import gzip
from   collections import defaultdict

def get_arguments():
    """get the arguments"""

    parser = argparse.ArgumentParser(
        description="Extract sites where allele distribution varies between the two phenotypes")
    parser.add_argument("-1", "--red1", help="Input depth file for red 1 pool", required=True)
    parser.add_argument("-2", "--red2", help="Input depth file for red 2 pool", required=True)
    parser.add_argument("-3", "--yellow1", help="Input depth file for yellow 1 pool", required=True)
    parser.add_argument("-4", "--yellow2", help="Input depth file for yellow 2 pool", required=True)

    args    = parser.parse_args()
    red1    = args.red1
    red2    = args.red2
    yellow1 = args.yellow1
    yellow2 = args.yellow2

    for p in [red1, red2, yellow1, yellow2]:
        assert os.path.isfile(p), f"Could not locate {p}"

    return red1, red2, yellow1, yellow2

def del_pipeline(red1: str, red2: str, yellow1: str, yellow2: str):
    """start the pipeline"""

    redmap    = make_color_map(red1, red2)
    yellowmap = make_color_map(yellow1, yellow2)

    if (len(redmap) >= len(yellow2)):
        chroms = list(yellowmap.keys())
    else:
        chroms = list(redmap.keys())

    red_fh    = open("Red_Sites_Uniq.tsv", 'w')
    yellow_fh = open("Yellow_Sites_Uniq.tsv", 'w')

    for chr in chroms:
        print("Comparing sites for", chr)
        red_sites    = redmap[chr]
        yellow_sites = yellowmap[chr]
        red_unq    = red_sites.difference(yellow_sites)
        yellow_unq = yellow_sites.difference(red_sites)
        for r in red_unq:
            red_fh.write(f"{chr}\t{r}\n")
        for y in yellow_unq:
            yellow_fh.write(f"{chr}\t{y}\n")
        del redmap[chr]
        del yellowmap[chr]
        print("Finished with", chr)
    
    red_fh.close()
    yellow_fh.close()

def make_color_map(file1: str, file2: str):
    """create a single map for all positions for all chroms"""
    
    cmap = defaultdict(set)

    for f in [file1, file2]:
        print("Starting to parse", f)
        fh = gzip.open(f, "rt") if (f.endswith(".gz")) else open(f, 'r')
        for line in fh:
            fields = line.split('\t')
            cmap[fields[0]].add(int(fields[1]))
        fh.close()
        print("Finished reading in", f)

    return cmap

def main():
    """Parse, plot, and pray for good results"""

    # get arguments
    red1, red2, yellow1, yellow2 = get_arguments()
    # start pipeline
    del_pipeline(red1, red2, yellow1, yellow2)


if __name__ == "__main__":
    main()
