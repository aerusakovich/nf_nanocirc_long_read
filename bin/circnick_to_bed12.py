#!/usr/bin/env python3
"""
circnick_to_bed12.py

Converts circnick-lrs output to BED12 format.
Requires three input files:
  - combined_reads.circRNA_candidates.annotated.txt  (one row per circRNA)
  - combined_reads.circ_circRNA_exon_usage_length_of_exons.txt (one row per exon)
  - combined_reads.introns.uniq.exon_remove.coverage.onlyCirc.novelExonMap.intronCov.bed
    (used to recover exon structure for circRNAs with unresolved exon_usage rows)

Usage:
    python3 circnick_to_bed12.py annotated.txt exon_usage.txt intron_cov.bed output.bed12
"""

import sys


def parse_args():
    if len(sys.argv) != 5:
        print("Usage: circnick_to_bed12.py <annotated.txt> <exon_usage.txt> <intron_cov.bed> <output.bed12>")
        sys.exit(1)
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input files -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input files -> {}\n{}: '{}'".format(
            error, context, context_str
        )
    print(error_str)
    sys.exit(1)


def isoform_to_blocks(exon_list, chrom_start, chrom_end):
    """
    Convert a list of (exon_start, exon_end) tuples to BED12 block columns.

    Exon coordinates are clamped to BSJ boundaries [chrom_start, chrom_end]
    to handle cases where circnick reports exons extending outside its own
    reported BSJ coordinates.

    Returns tuple of three strings (block_count, block_sizes, block_starts)
    or None if all exons are outside the BSJ span.
    """
    sizes  = []
    starts = []

    for exon_start, exon_end in exon_list:
        # clamp to BSJ boundaries
        exon_start = max(exon_start, chrom_start)
        exon_end   = min(exon_end,   chrom_end)

        if exon_start >= exon_end:
            continue   # exon completely outside BSJ — skip

        sizes.append(exon_end - exon_start)
        starts.append(exon_start - chrom_start)

    if not sizes:
        return None   # all exons were outside BSJ

    block_count  = str(len(sizes))
    block_sizes  = ",".join(str(s) for s in sizes)  + ","
    block_starts = ",".join(str(s) for s in starts) + ","

    return block_count, block_sizes, block_starts


def read_annotated(path):
    """
    Read annotated.txt into a dict keyed by circRNA ID.

    Returns:
        { 'mmu_circ_0002_Reg3a': {'chr': 'chr6', 'start': 78380708,
                                   'end': 78382351, 'strand': '+',
                                   'reads': '54'}, ... }
    """
    result = {}

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")

            if cols[0] == "internal_circRNA_name":
                continue

            if len(cols) < 7:
                print(f"WARNING annotated.txt line {lineno}: expected >=7 cols, got {len(cols)} — skipping")
                continue

            circ_id     = cols[0]
            chrom       = cols[1]
            start       = int(cols[2])
            end         = int(cols[3])
            reads       = cols[5]

            strand_gene = cols[6]
            strand      = strand_gene[0]

            if strand not in ("+", "-"):
                print(f"WARNING annotated.txt line {lineno}: unexpected strand '{strand}' — setting to '.'")
                strand = "."

            result[circ_id] = {
                "chr":    chrom,
                "start":  start,
                "end":    end,
                "strand": strand,
                "reads":  reads,
            }

    return result


def read_exon_usage(path):
    """
    Read exon_usage file into a dict keyed by circRNA ID.

    Returns:
        result  : { circ_id: [(start, end), ...] }
        failed  : set of circ_ids that had at least one unresolvable exon row
    """
    result = {}
    failed = set()

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")

            if cols[0] == "Internal circRNA IDs":
                continue

            if len(cols) < 7:
                print(f"WARNING exon_usage line {lineno}: expected >=7 cols, got {len(cols)} — skipping")
                continue

            circ_id = cols[0]

            if not cols[5].strip() or not cols[6].strip():
                failed.add(circ_id)
                continue

            try:
                exon_start = int(cols[5])
                exon_end   = int(cols[6])
            except ValueError:
                failed.add(circ_id)
                continue

            if circ_id not in result:
                result[circ_id] = []
            result[circ_id].append((exon_start, exon_end))

    for circ_id in result:
        result[circ_id].sort(key=lambda x: x[0])

    print(f"[read_exon_usage] {len(result)} circRNAs with resolved exons, "
          f"{len(failed)} circRNAs with unresolved exon rows")

    return result, failed


def read_intron_coverage(path):
    """
    Read intron coverage file into a dict keyed by circRNA ID.

    Returns:
        { 'mmu_circ_0004_Heatr1': [(12432779, 12433611), ...], ... }
    """
    result = {}

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")

            if len(cols) < 11:
                continue

            circ_id = cols[10].strip()
            if not circ_id or circ_id == ".":
                continue

            try:
                intron_start = int(cols[1])
                intron_end   = int(cols[2])
            except ValueError:
                print(f"WARNING intron_cov line {lineno}: non-numeric coords — skipping")
                continue

            if circ_id not in result:
                result[circ_id] = []
            result[circ_id].append((intron_start, intron_end))

    for circ_id in result:
        result[circ_id].sort(key=lambda x: x[0])

    return result


def derive_exons_from_introns(bsj_start, bsj_end, introns):
    """
    Derive exon coordinates by filling gaps between introns within the BSJ span.

    Example:
        bsj_start = 12432628, bsj_end = 12434051
        introns   = [(12432779, 12433611), (12433770, 12433916)]
        → exons: [(12432628, 12432779), (12433611, 12433770), (12433916, 12434051)]
    """
    exons   = []
    current = bsj_start

    for intron_start, intron_end in introns:
        if intron_start >= bsj_end or intron_end <= bsj_start:
            continue

        if current < intron_start:
            exons.append((current, intron_start))

        current = intron_end

    if current < bsj_end:
        exons.append((current, bsj_end))

    return exons


def convert(annotated_file, exon_file, intron_file, out_file):
    """
    Join all files and write BED12.
    """
    annotated          = read_annotated(annotated_file)
    exon_usage, failed = read_exon_usage(exon_file)
    intron_cov         = read_intron_coverage(intron_file)

    written         = 0
    recovered       = 0
    skipped_no_exon = 0
    skipped_other   = 0

    with open(out_file, "w") as out:
        for circ_id, info in annotated.items():

            chrom_start = info["start"]
            chrom_end   = info["end"]

            if circ_id in failed or circ_id not in exon_usage:
                # try intron-based recovery
                introns = intron_cov.get(circ_id, [])
                exons   = derive_exons_from_introns(chrom_start, chrom_end, introns)

                if not exons:
                    # fallback: single block spanning full BSJ
                    exons = [(chrom_start, chrom_end)]
                    print(f"WARNING: {circ_id} has no exon or intron data — writing as single block")
                    skipped_no_exon += 1
                else:
                    recovered += 1
            else:
                exons = exon_usage[circ_id]

            result = isoform_to_blocks(exons, chrom_start, chrom_end)
            if result is None:
                print(f"WARNING: {circ_id} all exons outside BSJ — skipping")
                skipped_other += 1
                continue

            block_count, block_sizes, block_starts = result

            bed12_line = "\t".join([
                info["chr"],
                str(chrom_start),
                str(chrom_end),
                circ_id,
                info["reads"],
                info["strand"],
                str(chrom_start),
                str(chrom_end),
                "0",
                block_count,
                block_sizes,
                block_starts,
            ])

            out.write(bed12_line + "\n")
            written += 1

    print(f"Done: {written} records written "
          f"({recovered} recovered from intron data), "
          f"{skipped_no_exon} single-block fallbacks, "
          f"{skipped_other} skipped (all exons outside BSJ)")


if __name__ == "__main__":
    annotated_file, exon_file, intron_file, out_file = parse_args()
    convert(annotated_file, exon_file, intron_file, out_file)