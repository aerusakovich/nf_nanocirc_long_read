#!/usr/bin/env python3
"""
circnick_liftover_introns.py

Rejoins lifted intron coordinates from liftOver output back into
the original intronCov.bed file format.

Called by CIRCNICK_LIFTOVER module after liftOver runs on introns.

Usage:
    python3 circnick_liftover_introns.py \\
        --lifted_introns  intron_lifted.bed \\
        --orig_introns    combined_reads.introns...intronCov.bed \\
        --sample          SAMPLE_ID
"""

import os
import sys
import argparse


def parse_args(args=None):
    Description = "Rejoin liftOver results with original circnick intron coverage file."
    Epilog = "Example usage: circnick_liftover_introns.py --lifted_introns a.bed --orig_introns b.bed --sample S1"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--lifted_introns", required=True, help="BED4 output from liftOver for introns.")
    parser.add_argument("--orig_introns",   required=True, help="Original intronCov.bed file.")
    parser.add_argument("--sample",         required=True, help="Sample name for output filename.")
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input files -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input files -> {}\n{}: '{}'".format(
            error, context, context_str
        )
    print(error_str)
    sys.exit(1)


def main():
    args = parse_args()

    # load lifted coords keyed by synthetic name: circ_id_INTRON_linenum
    lifted = {}
    with open(args.lifted_introns) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            lifted[cols[3]] = (cols[0], cols[1], cols[2])

    out_path = "{}_lifted_intron_cov.bed".format(args.sample)
    written  = 0
    skipped  = 0

    with open(args.orig_introns) as fh, open(out_path, "w") as out:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 11:
                continue
            circ_id = cols[10].strip()
            key     = "{}_INTRON_{}".format(circ_id, lineno)

            if key not in lifted:
                sys.stderr.write("WARNING: intron {} not lifted — skipping\n".format(key))
                skipped += 1
                continue

            cols[0], cols[1], cols[2] = lifted[key]
            out.write("\t".join(cols) + "\n")
            written += 1

    print("[circnick_liftover_introns] {} introns written, {} skipped → {}".format(
        written, skipped, out_path
    ))


if __name__ == "__main__":
    main()
