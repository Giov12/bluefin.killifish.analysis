#!/bin/env python3

import gzip
import argparse
import os
import sys

vcf          =  ''
chrom        = ''
target_samps = list()
left         = -1
right        = float("inf")
target_alles = list()
   
def get_arguments() -> int:
    """function to set the global variables"""

    global vcf, chrom, left, right, target_samps

    parser = argparse.ArgumentParser(description="Find fixed or private alleles for killifish samples of interest")
    parser.add_argument("-v", "--vcf",   required=True, type=str, help="vcf file with common and white stickleback")
    parser.add_argument("-s", "--samples", required=True, type=str, help="sample IDs for the samples of interest", nargs='*')
    parser.add_argument("-c", "--chrom", required=True, type=str, help="locus ID of region of interest")
    parser.add_argument("-l", "--left", type=float, help="left base pair position of region of interest", default=left)
    parser.add_argument("-r", "--right", type=float, help="right base pair position of region of interest", default=right)

    args    = parser.parse_args()
    vcf     = args.vcf
    samples = args.samples
    chrom   = args.chrom
    left    = args.left
    right   = args.right

    # make sure input files exist
    assert os.path.isfile(vcf), f"Could not locate {vcf}"
    assert left > 0, "--left must be greater than 0"
    assert right > left, "--right must be greater than --left"

    for sample in samples:
        sample = sample.strip().strip("\",")
        target_samps.append(sample)

    assert len(target_samps) > 0, "at least one sample must be provided"

    return 0

def get_sites() -> int:
    """Find the load the fixed differences between commons vs whites"""

    global vcf, left, right, chrom, target_samps

    fh      = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')
    targets = list()

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields = line.strip().split('\t')
                cnt    = 0
                for i in range(9, len(fields)):
                    targets.append(fields[i] in target_samps)
                    cnt += 1 if targets[-1] else 0
                if (len(target_samps) != cnt):
                    msg = "Did not find all samples in vcf file. Terminating.."
                    sys.exit(msg)
            continue
        fields   = line.strip().split('\t')
        if (chrom != fields[0]):
            continue
        pos = int(fields[1])
        if ((left <= pos <= right) == False):
            continue

        tAlleles = set() # target alleles
        oAlleles = set() # other alleles
        ref_AL   = fields[3]
        alt_AL   = fields[4]
        has_miss = False

        for i in range(9, len(fields)):
            is_target = targets[i - 9]
            info      = fields[i]
            idx       = info.find(':')
            geno      = info[:idx]
            alleles   = geno.split('/')
            if (len(alleles) != 2 and geno.count('|') > 0):
                alleles = geno.split('|')
            if (len(alleles) != 2):
                msg = f"Encountered an ilformated entry in the vcf\n" + \
                      f"Offending line: {line}"
                sys.exit(msg)
            if ('.' in alleles):
                has_miss = True
                break
            allele1 = ref_AL if (alleles[0] == '0') else alt_AL
            allele2 = ref_AL if (alleles[1] == '0') else alt_AL
            if (is_target):
                tAlleles.add(allele1 + '_' + allele2)
            else:
                oAlleles.add(allele1 + '_' + allele2)

        if (has_miss):
            continue

        valid = False

        # now check that they don't match
        if (len(tAlleles) == 1):
            # if no alleles are shared
            valid = tAlleles.isdisjoint(oAlleles)
            if (valid):
                target_alles.append((pos, tAlleles, oAlleles))

    fh.close()

    if (len(target_alles) == 0):
        msg = f"Did not find any fixed or private differences at the target locus"
        sys.exit(msg)
    
    print(f"Found a total of {len(target_alles)} fixed or private sites at the target locus")

    return 0

def write_sites() -> int:
    """check the eQTL and see if it is a fixed difference"""

    global chrom, target_alles

    ofh = open("Fixed.Private.Alleles.tsv", 'w')

    for alleles in target_alles:
        pos = alleles[0]
        tAL = ','.join(alleles[1]) # target alleles
        oAL = ','.join(alleles[2]) # other alleles
        ln  = f"{chrom}\t{pos}\t{tAL}\t{oAL}\n"
        ofh.write(ln) 
    ofh.close()

    return 0

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # find sites of interest
    get_sites()

    # write out the findings
    write_sites()

if __name__ == "__main__":
    main()