#!/usr/bin/env python3
"""
add_class_codes.py
------------------
Merges gffcompare class codes and reference annotation info into a
circRNA confidence TSV produced by merge_circrna.py / merge_exon_based.py.

AGAT's agat_convert_bed2gff.pl assigns sequential integers (1, 2, 3…) as
GFF feature IDs and stores the original BED name field as the GFF 'Name'
attribute.  gffcompare therefore reports those integers as qry_id / qry_gene_id
in the .tmap file — not the BSJ id.

This script resolves the mismatch by parsing the AGAT GFF to build an
integer_ID → Name (bsj_id) lookup before joining with the confidence TSV.

New columns appended to the TSV:
    class_code   – gffcompare class code (= c j o e i x u s p r .)
    ref_gene_id  – best-matching reference gene ID   (or '.' if intergenic)
    ref_id       – best-matching reference transcript (or '.' if none)
"""

import argparse
import csv
import re
import sys


def parse_gff_id_to_name(gff_path: str) -> dict:
    """
    Parse AGAT GFF output and return {ID_value: Name_value} for every
    feature that has both an ID and a Name attribute.

    AGAT writes:  ID=1;Name=chr1:1000-2000:+
    We build:     {'1': 'chr1:1000-2000:+'}
    """
    id_to_name = {}
    id_re   = re.compile(r'(?:^|;)ID=([^;]+)')
    name_re = re.compile(r'(?:^|;)Name=([^;]+)')

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            attrs = parts[8]
            id_m   = id_re.search(attrs)
            name_m = name_re.search(attrs)
            if id_m and name_m:
                id_to_name[id_m.group(1)] = name_m.group(1)

    return id_to_name


def parse_tmap(tmap_path: str, id_to_name: dict) -> dict:
    """
    Parse a gffcompare .tmap file and return a dict keyed by bsj_id.

    gffcompare .tmap columns (0-based):
        0  ref_gene_id
        1  ref_id
        2  class_code
        3  qry_gene_id   ← AGAT integer ID
        4  qry_id        ← AGAT integer ID (same here since no mRNA level)

    We translate qry_gene_id / qry_id → bsj_id via id_to_name.
    """
    lookup = {}  # bsj_id → info dict

    with open(tmap_path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 5:
                continue

            ref_gene_id = row[0]
            ref_id      = row[1]
            class_code  = row[2]
            qry_gene_id = row[3]
            qry_id      = row[4]

            info = {
                'class_code' : class_code,
                'ref_gene_id': ref_gene_id,
                'ref_id'     : ref_id,
            }

            # Translate AGAT integer IDs to bsj_ids via GFF Name attributes
            for agat_id in (qry_id, qry_gene_id):
                bsj_id = id_to_name.get(agat_id, agat_id)  # fallback: use as-is
                lookup[bsj_id] = info

    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tsv',    required=True, help='Confidence TSV from merge step')
    parser.add_argument('--tmap',   required=True, help='gffcompare .tmap output')
    parser.add_argument('--gff',    required=True, help='AGAT GFF output (used to resolve integer IDs to bsj_ids)')
    parser.add_argument('--output', required=True, help='Output annotated TSV path')
    args = parser.parse_args()

    id_to_name  = parse_gff_id_to_name(args.gff)
    tmap_lookup = parse_tmap(args.tmap, id_to_name)

    n_matched = 0
    n_total   = 0
    n_na      = 0

    with open(args.tsv) as in_fh, open(args.output, 'w', newline='') as out_fh:
        reader = csv.reader(in_fh, delimiter='\t')
        writer = csv.writer(out_fh, delimiter='\t')

        for row in reader:
            if not row:
                writer.writerow(row)
                continue

            if row[0].startswith('#') or row[0] == 'chrom':
                writer.writerow(row + ['class_code', 'ref_gene_id', 'ref_id'])
                continue

            n_total += 1

            # bsj_id is column index 4 in both BSJ-based and exon-based TSV formats
            bsj_id = row[4] if len(row) > 4 else ''

            info = tmap_lookup.get(bsj_id, {
                'class_code' : 'NA',
                'ref_gene_id': '.',
                'ref_id'     : '.',
            })

            if info['class_code'] == 'NA':
                n_na += 1
            else:
                n_matched += 1

            writer.writerow(row + [info['class_code'], info['ref_gene_id'], info['ref_id']])

    if n_total > 0:
        pct = 100.0 * n_matched / n_total
        print(
            f"[add_class_codes] Annotated {n_matched}/{n_total} circRNAs "
            f"({pct:.1f}%) with gffcompare class codes. Unmatched: {n_na}",
            file=sys.stderr
        )


if __name__ == '__main__':
    main()
