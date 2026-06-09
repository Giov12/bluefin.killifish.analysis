#!/bin/env python3

import argparse
import os
import gzip
import matplotlib.pyplot as plt
from glob import glob

cov_path  = ''
cov_files = list()
fai_file  = ''
chroms    = dict() # str -> int
figures   = dict() # str -> plot tuple

def get_arguments() -> int:
    """function to get the users input"""

    global cov_path, fai_file

    parser = argparse.ArgumentParser(description="Plot the coverages from a set of depth files")
    parser.add_argument("-d", "--dir", required=True, type=str, help="directory containing the depth files")
    parser.add_argument("-F", "--fai", required=True, type=str, help="fasta index file with the lengths of the chromosomes")

    args = parser.parse_args()

    assert os.path.isdir(args.dir),  f"Could not locate {args.dir}"
    assert os.path.isfile(args.fai), f"Could not locate {args.fai}"

    cov_path = args.dir
    fai_file = args.fai

    return 0

def get_coverages() -> int:
    """get the depth files from the provided path"""

    global cov_path, cov_files

    if (cov_path[-1] == '/'):
        cov_files = glob(f"{cov_path}*.depth.gz")
    else:
        cov_files = glob(f"{cov_path}/*.depth.gz")

    assert len(cov_files) > 0, f"Could not locate any .depth files in {cov_path}"
    cov_files.sort()

    return 0

def load_chrom_lengths() -> int:
    """read in the chromosome lengths from the index file"""

    global fai_file, chroms, figures

    fh        = open(fai_file, 'r')
    interval = 5_000_000

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        chrom  = fields[0]
        length = int(fields[1])
        chroms[chrom] = length
        fig, ax = plt.subplots(figsize = (12, 6))
        ax.set_ylabel("Coverage (x)")
        ax.set_xlim(0, length)

        if (length >= interval):
            ticks  = list()
            labels = list()
            for tick in range(interval, length + 1, interval):
                ticks.append(tick)
                label = tick // 1_000_000
                label = f"{label} Mb"
                labels.append(label)
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)

        figures[chrom] = (fig, ax)

    fh.close()

    return 0

def plot_coverage(chrom: str, positions: list[int], scores: list[float], color: str) -> int:
    """plot a single chromosome"""

    global chroms, figures

    if (chrom not in chroms or len(positions) == 0):
        positions.clear()
        scores.clear()
        return 1
    
    fig, ax = figures[chrom]
    ax.plot(positions, scores, color = color)
    # no longer needed in memory
    positions.clear()
    scores.clear()

    return 0

def process_coverages(cov_file: str, color: str) -> int:
    """process a depth file one at a time"""

    global chroms

    fh = gzip.open(cov_file, "rt") if cov_file.endswith(".gz") else open(cov_file, 'r')

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
        depth  = float(fields[2])

        if (curChrom == ''):
            curChrom = chrom
        if (curChrom == chrom):
            positions.append(pos)
            scores.append(depth)
        else:
            plot_coverage(curChrom, positions, scores, color)
            curChrom = chrom
            positions.append(pos)
            scores.append(depth)

    fh.close()

    # plot last chrom
    plot_coverage(curChrom, positions, scores, color)

    return 0

def close_figures() -> int:
    """save and close all the figures"""

    global figures

    for chrom, figs in figures.items():
        fig = figs[0]
        fig.savefig(f"{chrom}.depth.svg", dpi = 300, bbox_inches = "tight")
        fig.savefig(f"{chrom}.depth.png", dpi = 300, bbox_inches = "tight")
        plt.close(fig)
    return 0

def main() -> int:
    """entry point to this application"""

    # get the two inputs
    get_arguments()

    # load the chromosome lengths & get the depth files
    load_chrom_lengths()
    get_coverages()

    # now iteratively plot the figs
    global cov_files
    for i, cov_file in enumerate(cov_files):
        color = "#B22222" if i < 2 else "#4169E1"
        process_coverages(cov_file, color)

    close_figures()

    return 0

if __name__ == "__main__":
    main()
