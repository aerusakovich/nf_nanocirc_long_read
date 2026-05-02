#!/usr/bin/env python3
"""
circnick_liftover.py

After liftOver has converted circnick coordinates to a new genome build,
this script rejoins the lifted coordinates back into the original
circnick output files (annotated.txt and exon_usage.txt).

Called by the CIRCNICK_LIFTOVER Nextflow module after liftOver runs.

Usage:
    python3 circnick_liftover.py \\
        --lifted_annotated  annotated_lifted.bed \\
        --lifted_exons      exon_lifted.bed \\
        --orig_annotated    combined_reads.circRNA_candidates.annotated.txt \\
        --orig_exon_usage   combined_reads.circ_circRNA_exon_usage_length_of_exons.txt \\
        --sample            SAMPLE_ID
"""

import os
import sys
import argparse


def parse_args(args=None):
    Description = "Rejoin liftOver results with original circnick output files."
    Epilog = "Example usage: circnick_liftover.py --lifted_annotated a.bed --lifted_exons b.bed --orig_annotated c.txt --orig_exon_usage d.txt --sample S1"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument("--lifted_annotated",  required=True,
                        help="BED4 file output by liftOver for annotated.txt coords.")
    parser.add_argument("--lifted_exons",      required=True,
                        help="BED4 file output by liftOver for exon_usage coords.")
    parser.add_argument("--orig_annotated",    required=True,
                        help="Original circnick annotated.txt file.")
    parser.add_argument("--orig_exon_usage",   required=True,
                        help="Original circnick exon_usage file.")
    parser.add_argument("--sample",            required=True,
                        help="Sample name used in output filenames.")

    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input files -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input files -> {}\n{}: '{}'".format(
            error, context, context_str
        )
    print(error_str)
    sys.exit(1)


def read_lifted_bed(path):
    """
    Read a liftOver output BED4 file.
    Returns dict: name → (new_chrom, new_start, new_end)
    """
    lifted = {}

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            name = cols[3]
            lifted[name] = (cols[0], cols[1], cols[2])

    return lifted


def lift_annotated(orig_path, lifted_coords, failed_circs, out_path):
    """
    Rewrite annotated.txt replacing chr/start/end with lifted coordinates.
    Excludes any circRNA in failed_circs (had exon liftover failure).

    annotated.txt columns:
      0: internal_circRNA_name
      1: chr
      2: start
      3: end
      4: description
      5: BSJ_reads
      6: strandgene  (e.g. "+Vegfd")
      ...

    The liftOver key is the internal_circRNA_name (col 0).
    """
    written  = 0
    skipped  = 0

    with open(orig_path) as fh, open(out_path, "w") as out:
        for lineno, line in enumerate(fh):
            line = line.rstrip("\n")

            # keep header unchanged
            if lineno == 0:
                out.write(line + "\n")
                continue

            cols = line.split("\t")
            name = cols[0]

            # exclude if BSJ itself failed to lift
            if name not in lifted_coords:
                sys.stderr.write("WARNING: {} BSJ not lifted — skipping\n".format(name))
                skipped += 1
                continue

            # exclude if any exon failed to lift
            if name in failed_circs:
                skipped += 1
                continue

            new_chrom, new_start, new_end = lifted_coords[name]
            cols[1] = new_chrom
            cols[2] = new_start
            cols[3] = new_end

            out.write("\t".join(cols) + "\n")
            written += 1

    print("[lift_annotated] {} records written, {} skipped → {}".format(
        written, skipped, out_path
    ))


def lift_exon_usage(orig_path, lifted_coords, out_path):
    """
    Rewrite exon_usage.txt replacing exon start/end with lifted coordinates.

    exon_usage.txt columns:
      0: Internal circRNA IDs
      1: Exon used
      2: Exon covered by read
      3: Usage level
      4: exon name  (e.g. "Vegfd_chrX:164402006-164402650")
      5: start
      6: end
      7: length

    The liftOver key is a synthetic name we created:
      <circ_id>_EXON_<line_number>   (1-based, excluding header)

    Returns a set of circRNA IDs that had at least one exon fail to lift.
    These circRNAs are excluded entirely from the output.
    """
    written  = 0
    skipped  = 0
    exon_idx = 0   # counts data lines only (no header)

    # first pass — find which circRNAs have any failed exon
    failed_circs = set()
    with open(orig_path) as fh:
        for lineno, line in enumerate(fh):
            if lineno == 0:
                continue
            exon_idx += 1
            cols = line.rstrip("\n").split("\t")
            key  = "{}_EXON_{}".format(cols[0], exon_idx)
            if key not in lifted_coords:
                failed_circs.add(cols[0])
                sys.stderr.write(
                    "WARNING: exon {} not lifted — entire circRNA {} will be excluded\n".format(
                        key, cols[0]
                    )
                )

    if failed_circs:
        sys.stderr.write(
            "WARNING: {} circRNA(s) excluded due to partial liftover failure\n".format(
                len(failed_circs)
            )
        )

    # second pass — write only fully lifted circRNAs
    exon_idx = 0
    with open(orig_path) as fh, open(out_path, "w") as out:
        for lineno, line in enumerate(fh):
            line = line.rstrip("\n")

            if lineno == 0:
                out.write(line + "\n")
                continue

            exon_idx += 1
            cols    = line.split("\t")
            circ_id = cols[0]
            key     = "{}_EXON_{}".format(circ_id, exon_idx)

            # skip entire circRNA if any of its exons failed
            if circ_id in failed_circs:
                skipped += 1
                continue

            new_chrom, new_start, new_end = lifted_coords[key]
            cols[5] = new_start
            cols[6] = new_end
            try:
                cols[7] = str(int(new_end) - int(new_start))
            except (ValueError, IndexError):
                pass

            out.write("\t".join(cols) + "\n")
            written += 1

    print("[lift_exon_usage] {} exons written, {} skipped → {}".format(
        written, skipped, out_path
    ))

    return failed_circs


def main():
    args = parse_args()

    # ── load lifted coordinates ───────────────────────────────
    lifted_annotated = read_lifted_bed(args.lifted_annotated)
    lifted_exons     = read_lifted_bed(args.lifted_exons)

    print("[circnick_liftover] Lifted circRNAs: {}".format(len(lifted_annotated)))
    print("[circnick_liftover] Lifted exons:    {}".format(len(lifted_exons)))

    # ── exon usage first — identifies which circRNAs failed ──
    out_exon_usage = "{}_lifted_exon_usage.txt".format(args.sample)
    failed_circs   = lift_exon_usage(
        args.orig_exon_usage, lifted_exons, out_exon_usage
    )

    # ── annotated.txt — exclude any circRNA with failed exons ─
    out_annotated = "{}_lifted_annotated.txt".format(args.sample)
    lift_annotated(
        args.orig_annotated, lifted_annotated, failed_circs, out_annotated
    )

    # ── write failed list TSV ─────────────────────────────────
    failed_path = "{}_liftover_failed.tsv".format(args.sample)
    with open(failed_path, "w") as fh:
        fh.write("#circRNA_id\treason\n")
        for circ_id in sorted(failed_circs):
            fh.write("{}\tone_or_more_exons_not_lifted\n".format(circ_id))

    print("[circnick_liftover] {} circRNAs failed liftover → {}".format(
        len(failed_circs), failed_path
    ))


if __name__ == "__main__":
    main()
