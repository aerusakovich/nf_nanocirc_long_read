#!/usr/bin/env python
"""
liftover_bed.py — BED4 liftover using the Python liftover library.
Argument order matches UCSC liftOver: input.bed chain output.bed unmapped.bed
"""
import sys
from liftover import get_lifter


def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: liftover_bed.py <input.bed> <chain.gz> <output.bed> <unmapped.bed>")

    in_bed, chain_path, out_bed, unmapped_bed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

    converter = get_lifter(chain_path)
    lifted = []
    unmapped = []

    with open(in_bed) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                unmapped.append(line)
                continue

            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]

            try:
                start_hits = converter[chrom][start]
                end_hits   = converter[chrom][end - 1]  # BED end is exclusive
            except Exception:
                unmapped.append(line)
                continue

            if start_hits and end_hits:
                s_chrom, s_pos, s_strand = start_hits[0]
                e_chrom, e_pos, e_strand = end_hits[0]
                if s_chrom == e_chrom and s_strand == e_strand:
                    if s_strand == '+':
                        lifted.append('\t'.join([s_chrom, str(s_pos), str(e_pos + 1), name]))
                    else:
                        lifted.append('\t'.join([s_chrom, str(e_pos), str(s_pos + 1), name]))
                else:
                    unmapped.append(line)
            else:
                unmapped.append(line)

    with open(out_bed, 'w') as fh:
        for r in lifted:
            fh.write(r + '\n')

    with open(unmapped_bed, 'w') as fh:
        for line in unmapped:
            fh.write(line + '\n')


if __name__ == '__main__':
    main()
