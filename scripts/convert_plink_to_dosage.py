#!/usr/bin/env python

import sys
import os
import subprocess
import gzip
from argparse import ArgumentParser
from collections import defaultdict

# from itertools import izip

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Converts plink file format to dosage format for prediXcan"
    )
    parser.add_argument(
        "-b", "--bfile", dest="bfile", required=True,
        help="prefix for binary ped, bim, and bam files."
    )
    parser.add_argument(
        "-o", "--out", dest="out", required=False,
        help="prefix for output files", default="chr"
    )
    parser.add_argument(
        "-p", "--plink-binary", dest="plink", default="plink_linux/plink",
        help="path to plink (1.9) binary"
    )

    if len(sys.argv) == 2:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    # First we get the minor allele dosages for *all* chromosomes:
    subprocess.call([
        args.plink, '--bfile', args.bfile,
        '--recode', 'A-transpose', '--out', args.out
    ])

    # Now calculate the minor allele frequency:
    subprocess.call([
        args.plink, '--bfile', args.bfile, '--freq', '--out', args.out
    ])

    # Thanks to Adam Whiteside (https://github.com/adamcw) for
    # cleaning up this code to be more efficient
    diter = iter(x.split() for x in open(args.out + ".traw"))
    fiter = iter(x.split() for x in open(args.out + ".frq"))
    biter = iter(x.split() for x in open(args.bfile + ".bim"))

    buff = defaultdict(list)

    # First skip the header row
    # diter.next()
    # fiter.next()

    # Then process the lines
    for i, (dcols, fcols, bcols) in enumerate(zip(diter, fiter, biter)):
        if i == 0:
            # Skip header row
            continue
        # Combine columns as per 'dosage' format.
        # First we add the information columns for each rsID
        nline = fcols[0:2] + [bcols[3]] + [fcols[3]] + [fcols[2]] + [fcols[4]]

        # Next add the (additive linear) dosage data for the samples.
        # Impute missing data as 2*MAF.
        MAF = float(fcols[4])
        nline = nline + [str(MAF * 2) if e == "NA" else e for e in dcols[6:]]

        # Add line to write buffer
        buff[fcols[0]].append(" ".join(nline) + "\n")

    # Write out the buffer
    # for curCHR, lines in buff.items():
    #     with gzip.open(args.out + curCHR + ".txt.gz", "wb") as ofile:
    #         ofile.writelines(lines)
    for curCHR, lines in buff.items():
        with gzip.open(args.out + ".txt.gz", "wb") as ofile:
            for line in lines:
                ofile.write(line.encode('utf-8'))

    # Remove left over files
    os.remove(args.out + ".frq")
    os.remove(args.out + ".traw")
    os.remove(args.out + ".log")
