#!/bin/env python3

import argparse
import os
import sys
import gzip
import matplotlib.pyplot as plt
from collections import defaultdict

fst_file = ''
fai_file = ''
chroms   = dict() # str -> int

def get_arguments() -> int:
    """function to get the users input"""

    global fst_file, fai_file

    parser = argparse.ArgumentParser(description="Plot the Fst scores from sync file")
    parser.add_argument("-f", "--fst", required=True, type=str, help="fst file generated from popoolation2")
    parser.add_argument("-F", "--fai", required=True, type=str, help="fasta index file with the lengths of the chromosomes")

    args = parser.parse_args()

    assert os.path.isfile(args.fst), f"Could not locate {args.fst}"
    assert os.path.isfile(args.fai), f"Could not locate {args.fai}"

    fst_file = args.fst
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

def plot_fst(chrom: str, positions: list[int], scores: list[float]) -> int:
    """plot a single chromosome"""

    global chroms

    if (chrom not in chroms or len(positions) == 0):
        positions.clear()
        scores.clear()
        return 1
    
    chrom_len = chroms[chrom]
    plt.figure(figsize=(12, 6))
    plt.ylim(0, 1.0)
    plt.xlim(0, chrom_len)
    plt.scatter(positions, scores, color = "#2a52be", s=1, alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{chrom}.fst.pdf")
    plt.close("all")
    plt.clf()
    # no longer needed in memory
    positions.clear()
    scores.clear()

    return 0

def process_fst() -> int:
    """plot one chromosome at a time"""

    #
    #             column
    #  r1_vs_y1 = 6
    #  r1_vs_y2 = 7
    #  r2_vs_y1 = 8
    #  r2_vs_y2 = 9
    #

    global fst_file

    fh        = gzip.open(fst_file, "rt") if fst_file.endswith(".gz") else open(fst_file, 'r')
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
        record = fields[6]
        fst    = float(record.split('=')[1])

        if (curChrom == ''):
            curChrom = chrom
        if (curChrom == chrom):
            positions.append(pos)
            scores.append(fst)
        else:
            plot_fst(chrom, positions, scores)
            curChrom = chrom
            positions.append(pos)
            scores.append(fst)

    fh.close()

    # plot last chrom
    plot_fst(chrom, positions, scores)

    return 0

def main():
    """entry point to this application"""

    # get the two input files
    get_arguments()

    # load the chromosome lengths
    load_chrom_lengths()

    # now iteratively plot the figs
    process_fst()


if __name__ == "__main__":
    main()
