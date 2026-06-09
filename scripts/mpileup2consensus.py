#!/bin/python

import os
import argparse
import gzip
import textwrap


# set a global dict to get IUPAC

iupac = {'A'    : 'A',
         "AC"   : 'M',
         "ACG"  : 'V',
         "ACGT" : 'N',
         "ACT"  : 'H',
         "AG"   : 'R',
         "AGT"  : 'D',
         "AT"   : 'W',
         'C'    : 'C',
         "CG"   : 'S',
         "CGT"  : 'B',
         "CT"   : 'Y',
         'G'    : 'G',
         "GT"   : 'K',
         'T'    : 'T'
         }

def get_arguments():
    """get the arguments"""

    parser = argparse.ArgumentParser(
        description="generate a consensus sequence from the mpileup file")
    parser.add_argument(
        "-m", "--mpileup", help="mpileup file generated from samtools", required=True)
    parser.add_argument(
        "-f", "--fraction", help="fraction of sites required to consider a nucleotide", type=float, default=0.10)
    parser.add_argument(
        '-o', "--outfile", help="name of output file", type=str, default=''
    )

    args    = parser.parse_args()
    mpileup = args.mpileup
    fract   = args.fraction
    outfile = args.outfile

    assert os.path.isfile(mpileup), f"Could not location {mpileup}"
    assert (fract >= 0.0 and fract <= 1.0), f"{fract} was not a legal fraction"

    return mpileup, fract, outfile


def parse_mpileup(mpileup: str, fract: float, outfile: str):
    """a single function to go through the entire file"""

    """
    overall structure will look like this
    [[0, 0, 0, 0]]
    where the positions will be sorted by alpha
    [A, C, G, T]
    """

    characters = []
    fh         = gzip.open(mpileup, "rt") if mpileup.endswith(".gz") else open(mpileup, 'r')
    nuc_chars  = ['A', 'C', 'G', 'T']
    s, e       = float("inf"), -float("inf")
    ref        = ''

    for line in fh:
        fields = line.split('\t')
        sites  = fields[4]
        ref    = fields[0]
        count  = int(fields[3])
        nucs   = [0, 0, 0, 0]
        pos    = int(fields[1])
        s      = min(s, pos)
        e      = max(e, pos)

        for site in sites:
            if (site == 'A') or (site == 'a'):
                nucs[0] += 1
            elif (site == 'C') or (site == 'c'):
                nucs[1] += 1
            elif (site == 'G') or (site == 'g'):
                nucs[2] += 1
            elif (site == 'T') or (site == 't'):
                nucs[3] += 1
        
        num   = int(count * fract)
        chars = ''
        for i in range(4):
            if (nucs[i] >= num and nucs[i] > 0):
                chars += nuc_chars[i]
        if (chars == ''):
            characters.append('N')
        else:
            characters.append(iupac[chars])
    
        
    consensus = ''.join(characters)

    fh.close()

    if ref == '':
        ref = "Ref_Seq"

    if (outfile == ''):    
        outfile = f"{ref}_{s}-{e}_{fract}_consensue.fa"
    
    seq_header = f">{ref}_{s}-{e}\n"
    
    with open(outfile, 'w') as fh:
        fh.write(seq_header)
        fh.write(textwrap.fill(consensus, width=70))
        fh.write('\n')

def main():
    """get the inputs and create the consensus"""

    # get arguments
    mpileup, fract, outfile = get_arguments()

    # a single function to perform everything
    parse_mpileup(mpileup, fract, outfile)

if __name__ == "__main__":
    main()
