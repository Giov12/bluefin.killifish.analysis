#!/bin/env python3

import gzip
import argparse
import os
import glob

path    = ''
depths  = list()
samples = dict()

EXONS = [
        (25607727,25607913), (25608083,25608212), (25610157,25610203), (25610290,25610366), (25612291,25612383), 
        (25612467,25612606), (25614250,25614338), (25614422,25614565), (25614823,25614930), (25616043,25616226), 
        (25616355,25616478), (25616584,25616679)
        ]

class Exon:
    def __init__(self) -> None:
        self.depth = 0
        self.count = 0

class Sample:
    def __init__(self, filename: str) -> None:
        self.id     = os.path.basename(filename).split('.')[0]
        self.tissue = "fin" if self.id[-1] == 'F' else "eye"
        self.exons  = {f'Exon-{i}': Exon() for i in range(1, 13)}

def get_arguments() -> int:
    """get the path to the depth files"""

    global path

    parser = argparse.ArgumentParser(description="Get the Average Coverages for each depth file for the exonic regions")
    parser.add_argument("-p", "--path", required=True, type=str)

    args = parser.parse_args()
    path = args.path
    assert os.path.isdir(path), f"Could not locate {path}"

    return 0

def get_depths() -> int:
    """get the *.depth.gz files in the directory"""

    global path, depths

    if (path[-1] == '/'):
        depths = glob.glob(f"{path}*depth.gz")
    else:
        depths = glob.glob(f"{path}/*depth.gz")

    assert len(depths) > 0, f"Could not locate depth files in {path}"

    return 0

def process_depth(depth_file: str) -> int:
    """process each sample"""

    global samples, EXONS

    sample = Sample(depth_file)

    fh = gzip.open(depth_file, "rt")

    for line in fh:
        fields   = line.strip().split('\t')
        depth    = int(fields[2])
        position = int(fields[1])

        if (depth == 0):
            continue
        for i, exon in enumerate(EXONS, start=1):
            if (exon[0] <= position <= exon[1]):
                sample.exons[f"Exon-{i}"].count += 1
                sample.exons[f"Exon-{i}"].depth += depth
                break

    fh.close()

    samples[sample.id] = sample

    return 0

def write_output() -> int:
    """write a tsv file"""

    global samples

    ids = list(samples.keys())
    ids.sort()
    keys = [f"Exon-{i}" for i in range(1, 13)]

    fh     = open("Coverage.stats.tsv", 'w')
    header = f"Sample.ID\tTissue\t" + '\t'.join(keys) + "\tOutside.Deletion\tInside.Deletion\n"
    fh.write(header)

    for id_ in ids:
        sample   = samples[id_]
        outline  = [sample.id, sample.tissue]
        sections = [[0, 0], [0,0]]
        for i, key in enumerate(keys, start=1):
            exon = sample.exons[key]
            if (exon.count == 0):
                outline.append('0')
                continue
            avg = round(exon.depth / exon.count, 2)
            outline.append(str(avg))
            if (i <= 4):
                sections[0][0] += 1
                sections[0][1] += avg
            else:
                sections[1][0] += 1
                sections[1][1] += avg
        outside = 0 # outside vs inside the deletion
        inside  = 0 
        if (sections[0][0] > 1):
            outside = round(sections[0][1] / sections[0][0], 2)
        if (sections[1][0] > 1):
            inside = round(sections[1][1] / sections[1][0], 2)
        outline.append(str(outside))
        outline.append(str(inside))
        outline = '\t'.join(outline) + '\n'
        fh.write(outline)

    fh.close()

    return 0

def main():
    """entry point to the pipeline"""

    # get the path
    get_arguments()

    # retrieve the input files
    get_depths()

    # loop through and process each sample
    global depths
    for depth in depths:
        process_depth(depth)

    # write the stats
    write_output()

if __name__ == "__main__":
    main()