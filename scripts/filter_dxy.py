#!/bin/python3

from glob import glob
import gzip
import os

sdir = "/projects/catchenlab/gm33/poolseq/Popoolation2/20230208_RedvsYellow"

dxy_files = glob(f"{sdir}/Dxy*sync.gz")

def inspect_col(col1: str, col2: str, min_val = 4):
    
    cnt1 = col1.split(':')
    cnt2 = col2.split(':')
    vals1 = [int(c) for c in cnt1 if c != '0']
    vals2 = [int(c) for c in cnt2 if c != '0']

    if len(vals1) == 1:
        if len(vals2) == 1:
            keep1 = False
            keep2 = False
        elif len(vals2) == 2:
            for i in vals2:
                if i < min_val:
                    keep2 = False
                    break
                else:
                    keep2 = True
            if vals1[0] < min_val:
                keep1 = False
            else:
                keep1 = True
    elif len(vals2) == 1:
        for i in vals1:
            if i < min_val:
                keep1 = False
                break
            else:
                keep1 = True
        if vals2[0] < min_val:
            keep2 = False
        else:
            keep2 = True

    return keep1, keep2


for d_file in dxy_files:
    dfname = os.path.basename(d_file)
    outfile = "Filtered_" + dfname
    outfh = open(outfile, 'w')
    with gzip.open(d_file, 'rt') as fh:
        for line in fh:
            fields = line.strip('\n').split('\t')
            col1 = fields[2]
            col2 = fields[3]
            keep1, keep2 = inspect_col(col1, col2, 4)
            if keep1 and keep2:
                outfh.write(line)
    outfh.close()