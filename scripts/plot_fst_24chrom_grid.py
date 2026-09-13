#!/bin/env python3

import argparse
import os
import gzip
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

fst_file = ''
fontfam  = "DejaVu Sans"

# str -> (chrom length, Row, Column)
chroms_map = {"chr1":  (44090661, 0, 0),
              "chr10": (45801424,2,1),
              "chr11": (39918781, 2, 2),
              "chr12": (35890316, 2, 3),
              "chr13": (28028775, 3, 0),
              "chr14": (38801762, 3, 1),
              "chr15": (42821987, 3, 2),
              "chr16": (39540539, 3, 3),
              "chr17": (39649672, 4, 0),
              "chr18": (21809455, 4, 1),
              "chr19": (48135270, 4, 2),
              "chr2":  (47733668, 0, 1),
              "chr20": (45197809, 4, 3),
              "chr21": (48642572, 5, 0),
              "chr22": (37475815, 5, 1),
              "chr23": (38664211, 5, 2),
              "chr24": (48724087, 5, 3),
              "chr3":  (44242798, 0, 2),
              "chr4":  (47097095, 0, 3),
              "chr5":  (49292181, 1, 0),
              "chr6":  (39054102, 1, 1),
              "chr7":  (45720964, 1, 2),
              "chr8":  (46627064, 1, 3),
              "chr9":  (38401359, 2, 0)
}


def get_arguments() -> int:
    """function to get the users input"""

    global fst_file

    parser = argparse.ArgumentParser(description="Plot the Fst scores from sync file")
    parser.add_argument("-f", "--fst", required=True, type=str, help="fst file generated from popoolation2")

    args = parser.parse_args()

    assert os.path.isfile(args.fst), f"Could not locate {args.fst}"

    fst_file = args.fst

    return 0


def plot_fst(chrom: str, positions: list[int], scores: list[float], axes) -> int:
    """plot a single chromosome to its corresponding axes"""

    global chroms_map, fontfam
    
    if (len(positions) == 0):
        positions.clear()
        scores.clear()
        return 1
        
    # pull out the positional information
    chrom_tuple = chroms_map[chrom]
    chrom_len   = chrom_tuple[0]
    chrom_row   = chrom_tuple[1]
    chrom_col   = chrom_tuple[2]
    color       = "#f83817" if chrom_col % 2 else "#268ac5"
    axe         = axes[chrom_row][chrom_col]
    
    # start with the actual data
    # axe.scatter(positions, scores, color = color, s = 3, alpha = 0.5, linewidths = 0)
    axe.plot(positions, scores, color=color, linewidth=0.7, alpha=0.7)
    
    # only first column gets a y-axis
    axe.set_ylim(0, 0.06)
    axe.yaxis.set_major_locator(MultipleLocator(0.01)) # every 0.01
    if (chrom_col != 0):
        axe.set_yticklabels([])
        axe.tick_params(axis = 'y', length = 0)
    else:
        axe.tick_params(axis = 'y', labelsize = 12)
        axe.set_ylabel(r"$F_{ST}$", fontsize=13, fontfamily=fontfam)
    
    # now to add chrom length info
    axe.set_xlim(0, chrom_len)
    axe.xaxis.set_major_locator(MultipleLocator(10_000_000)) # every 10 Mb
    axe.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x / 1_000_000)} Mb"))
    axe.tick_params(axis="x", labelsize=11)
     
    # chromosome name on top of the subplot
    axe.set_title(chrom, fontsize=15, fontweight="bold", pad=4, fontfamily = fontfam)
    
    # change font for x & y axis
    for label in axe.get_xticklabels() + axe.get_yticklabels():
        label.set_fontfamily(fontfam)
    
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)
    
    # no longer needed in memory
    scores.clear()

    return 0

def process_fst() -> int:
    """plot one chromosome at a time"""

    #
    #   columns
    #  r1_vs_y1 = 2
    #  r1_vs_y2 = 3
    #  r2_vs_y1 = 4
    #  r2_vs_y2 = 95
    #

    global fst_file

    fh          = gzip.open(fst_file, "rt") if fst_file.endswith(".gz") else open(fst_file, 'r')
    fst_columns = {"r1_vs_y1": 2, "r1_vs_y2": 3, "r2_vs_y1": 4, "r2_vs_y2": 5}
    curChrom    = ''
    scores      = {p_comp : list() for p_comp in fst_columns.keys()} # pairwise comparisons
    positions   = list()

    # each pairwise comparison will be plotted onto its
    # own figure
    fig_axes_map = {}
    for p_comp in fst_columns.keys():
        fig_axes_map[p_comp] = plt.subplots(nrows=6, ncols=4, figsize=(16, 18), constrained_layout = True)

    for line in fh:
        if (len(line) == 0 or line[0] != 'c'):
            continue
        
        # just using red vs yellow 1
        fields = line.split('\t')
        chrom  = fields[0]
        pos    = int(fields[1])

        if (curChrom == ''):
            curChrom = chrom
        if (curChrom == chrom):
            positions.append(pos)
            # get each comparison's Fst value
            for name, col in fst_columns.items():
                # fst = float(fields[col].split('=')[1])
                fst = float(fields[col])
                scores[name].append(fst)
        else:
            # pass in each axes one-by-one
            for p_comp, (fig, axes) in fig_axes_map.items():
                plot_fst(curChrom, positions, scores[p_comp], axes)
            positions.clear() # clear out here
            curChrom = chrom
            positions.append(pos)

            # start collecting scores for the new current chrom
            for name, col in fst_columns.items():
                # fst = float(fields[col].split('=')[1])
                fst = float(fields[col])
                scores[name].append(fst)

    fh.close()

    # plot last chrom
    for p_comp, (fig, axes) in fig_axes_map.items():
        plot_fst(curChrom, positions, scores[p_comp], axes)
    positions.clear()

    # close out all the figures
    for p_comp, (fig, axes) in fig_axes_map.items():
        fig.savefig(f"all_chroms.{p_comp}.fst.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    return 0

def main():
    """entry point to this application"""

    # get the two input files
    get_arguments()

    # now iteratively plot each of the figures
    process_fst()

if __name__ == "__main__":
    main()
