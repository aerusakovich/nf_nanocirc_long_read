#!/usr/bin/env python3
"""
cirilong_to_bed12.py

Converts CIRI-long.info output to BED12 format.

CIRI-long.info uses 1-based GFF-like coordinates.
BED format is 0-based, so we subtract 1 from all start coordinates.
End coordinates are unchanged (GFF inclusive end == BED exclusive end).

Some records contain multiple isoforms separated by '|' in the isoform
field. Each isoform is written as a separate BED12 record with a
'_iso1', '_iso2', ... suffix appended to the name.

Usage:
    python3 cirilong_to_bed12.py CIRI-long.info output.bed12
"""

import sys


def parse_args():
    if len(sys.argv) != 3:
        print("Usage: cirilong_to_bed12.py <CIRI-long.info> <output.bed12>")
        sys.exit(1)
    return sys.argv[1], sys.argv[2]


def parse_attributes(attr_string):
    """
    Parse the GFF-like attributes column (col 8).

    Input looks like:
        circ_id "chr1:13672126-13672764"; isoform "13672126-13672764"; gene_name "Xkr9";

    Returns a dict:
        { 'circ_id': 'chr1:13672126-13672764', 'isoform': '13672126-13672764', ... }
    """
    result = {}
    for chunk in attr_string.strip().split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(" ", 1)
        if len(parts) != 2:
            continue
        key   = parts[0].strip()
        value = parts[1].strip().strip('"')
        result[key] = value
    return result


def isoform_to_blocks(isoform_str, chrom_start):
    """
    Convert a single isoform string to BED12 block columns.

    isoform_str uses 1-based coordinates e.g. "17750265-17750325,17750434-17750864"
    chrom_start is already 0-based (we subtracted 1 before calling this).

    Each exon start is converted from 1-based to 0-based by subtracting 1.
    Each exon end is kept as-is (GFF inclusive end == BED exclusive end).

    Returns tuple of three strings:
        block_count   e.g. "2"
        block_sizes   e.g. "61,431,"
        block_starts  e.g. "0,169,"
    or None if the isoform string is empty or unparseable.
    """
    exons  = [e.strip() for e in isoform_str.split(",") if e.strip()]
    if not exons:
        return None

    sizes  = []
    starts = []

    for exon in exons:
        parts = exon.split("-")
        if len(parts) != 2:
            print(f"WARNING: unexpected exon format '{exon}' — skipping exon")
            continue
        exon_start = int(parts[0]) - 1   # 1-based → 0-based
        exon_end   = int(parts[1])        # end stays the same
        sizes.append(exon_end - exon_start)
        starts.append(exon_start - chrom_start)

    if not sizes:
        return None

    block_count  = str(len(exons))
    block_sizes  = ",".join(str(s) for s in sizes)  + ","
    block_starts = ",".join(str(s) for s in starts) + ","

    return block_count, block_sizes, block_starts


def convert(info_file, out_file):
    """
    Read CIRI-long.info line by line, write BED12.

    Handles multi-isoform records where the isoform field contains
    multiple isoforms separated by '|'. Each isoform is written as
    a separate BED12 record with '_iso1', '_iso2', ... suffix.
    """
    written = 0
    skipped = 0

    with open(info_file) as fh, open(out_file, "w") as out:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                print(f"WARNING line {lineno}: expected 9 columns, got {len(cols)} — skipping")
                skipped += 1
                continue

            chrom       = cols[0]
            chrom_start = int(cols[3]) - 1   # 1-based GFF → 0-based BED
            chrom_end   = int(cols[4])        # end unchanged
            score       = cols[5]             # read count
            strand      = cols[6] if cols[6] in ('+', '-', '.') else '.'
            attrs       = parse_attributes(cols[8])

            isoform_field = attrs.get("isoform")
            if not isoform_field:
                print(f"WARNING line {lineno}: no isoform attribute — skipping")
                skipped += 1
                continue

            base_name = attrs.get("circ_id", f"{chrom}:{chrom_start}-{chrom_end}")

            # split on | to get individual isoforms
            isoforms = [iso.strip() for iso in isoform_field.split("|") if iso.strip()]

            for iso_idx, isoform in enumerate(isoforms):
                # add suffix only when multiple isoforms
                if len(isoforms) > 1:
                    name = f"{base_name}_iso{iso_idx + 1}"
                else:
                    name = base_name

                result = isoform_to_blocks(isoform, chrom_start)
                if result is None:
                    print(f"WARNING line {lineno} iso{iso_idx + 1}: could not parse isoform '{isoform}' — skipping")
                    skipped += 1
                    continue

                block_count, block_sizes, block_starts = result

                bed12_line = "\t".join([
                    chrom,
                    str(chrom_start),
                    str(chrom_end),
                    name,
                    score,
                    strand,
                    str(chrom_start),   # thickStart = thickEnd = start → no CDS
                    str(chrom_start),   # thickEnd
                    "0",                # itemRgb
                    block_count,
                    block_sizes,
                    block_starts,
                ])

                out.write(bed12_line + "\n")
                written += 1

    print(f"Done: {written} records written, {skipped} skipped")


if __name__ == "__main__":
    info_file, out_file = parse_args()
    convert(info_file, out_file)