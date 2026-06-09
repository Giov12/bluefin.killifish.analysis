#!/bin/env python3

import argparse
import os
import sys
import gzip
import matplotlib.pyplot as plt
from collections import defaultdict

cmh_file = ''
fai_file = ''
chroms   = dict() # str -> int

def get_arguments() -> int:
    """function to get the users input"""

    global cmh_file, fai_file

    parser = argparse.ArgumentParser(description="Plot the Fst scores from sync file")
    parser.add_argument("-c", "--cmh", required=True, type=str, help="cmh file with the 10th column containing -log(q) scores")
    parser.add_argument("-f", "--fai", required=True, type=str, help="fasta index file with the lengths of the chromosomes")

    args = parser.parse_args()

    assert os.path.isfile(args.cmh), f"Could not locate {args.cmh}"
    assert os.path.isfile(args.fai), f"Could not locate {args.fai}"

    cmh_file = args.cmh
    fai_file = args.fai

    return 0

def load_chrom_lengths() -> int:
    """read in the chromosome lengths from the index file"""

    global fai_file, chroms

    fh = open(fai_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        chrom  = fields[0]
        length = int(fields[1])
        chroms[chrom] = length
    fh.close()

    return 0

def plot_cmh(chrom: str, positions: list[int], scores: list[float]) -> int:
    """plot a single chromosome"""

    global chroms

    if (chrom not in chroms or len(positions) == 0):
        positions.clear()
        scores.clear()
        return 1
    
    chrom_len = chroms[chrom]
    plt.figure(figsize=(12, 6))
    plt.ylim(0, 32)
    plt.xlim(0, chrom_len)
    plt.scatter(positions, scores, color = "#2a52be", s=5, alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{chrom}.cmh.png", dpi = 300, bbox_inches = "tight")
    plt.savefig(f"{chrom}.cmh.svg", dpi = 300, bbox_inches = "tight")
    plt.close("all")
    plt.clf()
    # no longer needed in memory
    positions.clear()
    scores.clear()

    return 0

def process_cmh() -> int:
    """plot one chromosome at a time"""

    global cmh_file

    fh        = gzip.open(cmh_file, "rt") if cmh_file.endswith(".gz") else open(cmh_file, 'r')
    curChrom  = ''
    scores    = list()
    positions = list()

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        
        # just using red vs yellow 1
        fields = line.split('\t')
        chrom  = fields[0]
        pos    = int(fields[1])
        score  = float(fields[9]) if fields[9][0] != 'N' else 0.0

        if (curChrom == ''):
            curChrom = chrom
        if (curChrom == chrom):
            positions.append(pos)
            scores.append(score)
        else:
            plot_cmh(curChrom, positions, scores)
            curChrom = chrom
            positions.append(pos)
            scores.append(score)

    fh.close()

    # plot last chrom
    plot_cmh(curChrom, positions, scores)

    return 0

def main():
    """entry point to this application"""

    # get the two input files
    get_arguments()

    # load the chromosome lengths
    load_chrom_lengths()

    # now iteratively plot the figs
    process_cmh()


if __name__ == "__main__":
    main()
