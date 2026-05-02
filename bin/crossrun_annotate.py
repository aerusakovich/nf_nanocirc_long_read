#!/usr/bin/env python3
"""
crossrun_annotate.py

Reads the raw cross-run hybrid confidence TSV (output of add_isoform_confidence)
and writes two output files:

  <prefix>_confidence.tsv  — full intermediate: all columns from smart_merge,
                              plus n_samples, per-sample bsj/isoform consensus,
                              and circRNA type.

  <prefix>_clean.tsv       — wet-lab-friendly: only the essential columns,
                              rows filtered to bsj_confidence >= min_count.

  <prefix>.bed12           — BED12 filtered to match <prefix>_clean.tsv.

Usage:
    crossrun_annotate.py \\
        --input_tsv  GROUP_smart_consensus_hybrid_confidence.tsv \\
        --input_bed  GROUP_smart_consensus_hybrid.bed12 \\
        --sample_names sample1 sample2 sample3 \\
        --sample_tsvs  s1_confidence.tsv s2_confidence.tsv s3_confidence.tsv \\
        --gene_bed   genes.bed \\
        --exon_bed   exons.bed \\
        --bsj_tol 5 \\
        --min_count 2 \\
        --prefix GROUP_tier
"""

import argparse
import csv
import os
import subprocess
from collections import defaultdict

CLEAN_COLUMNS = [
    '#chrom', 'start', 'end', 'strand',
    'sel_block_count', 'sel_block_sizes', 'sel_block_starts',
    'bsj_id', 'bsj_confidence', 'isoform_confidence', 'n_samples',
    'type',
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input_tsv',    required=True, help='Raw hybrid confidence TSV (read-only)')
    p.add_argument('--input_bed',    required=True, help='Raw hybrid BED12 (read-only)')
    p.add_argument('--sample_names', nargs='+', required=True)
    p.add_argument('--sample_tsvs',  nargs='+', required=True)
    p.add_argument('--gene_bed',     required=True, help='Gene-level BED (6-col, from GTF)')
    p.add_argument('--exon_bed',     required=True, help='Exon-level BED (6-col, from GTF)')
    p.add_argument('--bsj_tol',      type=int, default=5)
    p.add_argument('--min_count',    type=int, required=True)
    p.add_argument('--prefix',       required=True, help='Output file prefix')
    return p.parse_args()


def parse_bsj_id(bsj_id):
    try:
        parts  = bsj_id.split(':')
        strand = parts[-1]
        chrom  = ':'.join(parts[:-2])
        start, end = map(int, parts[-2].split('-'))
        return chrom, start, end, strand
    except Exception:
        return None


def load_sample_index(tsv_path):
    """Return {(chrom, strand): [(start, end, bsj_consensus, isoform_consensus)]}."""
    index = defaultdict(list)
    try:
        with open(tsv_path) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                parsed = parse_bsj_id(row.get('bsj_id', '').split('|')[0])
                if parsed is None:
                    continue
                chrom, start, end, strand = parsed
                index[(chrom, strand)].append((
                    start, end,
                    row.get('bsj_consensus', ''),
                    row.get('isoform_consensus', ''),
                ))
    except FileNotFoundError:
        pass
    return index


def lookup(index, bsj_id, bsj_tol):
    parsed = parse_bsj_id(bsj_id)
    if parsed is None:
        return '', ''
    chrom, cs, ce, strand = parsed
    for ss, se, bc, ic in index.get((chrom, strand), []):
        if abs(ss - cs) <= bsj_tol and abs(se - ce) <= bsj_tol:
            return bc, ic
    return '', ''


# ── circRNA type classification ────────────────────────────────────────────────

def _spliced_length(block_sizes_str):
    return sum(int(x) for x in block_sizes_str.split(',') if x.strip())


def classify_types(rows, gene_bed, exon_bed):
    """
    Classify each circRNA as: eciRNA | EIciRNA | ciRNA | antisense | intergenic.

    Uses bedtools intersect -split -wo on BED12 circRNA records; overlap fractions
    are computed in Python from the returned overlap length divided by spliced length,
    avoiding the bedtools -f/-F fraction bug with -split.

      eciRNA     — same-strand gene; exon overlap covers 100% of spliced length
      EIciRNA    — same-strand gene; partial exon overlap (retains intronic content)
      ciRNA      — same-strand gene; no exon overlap (purely intronic)
      antisense  — opposite-strand gene overlap only
      intergenic — no gene overlap on either strand

    Returns {bsj_id: type_string}.
    """
    if not rows:
        return {}

    tmp_bed = '.type_class_tmp.bed'
    spliced = {}
    with open(tmp_bed, 'w') as fh:
        for row in rows:
            chrom  = row.get('#chrom', '')
            start  = row.get('start', '')
            end    = row.get('end', '')
            strand = row.get('strand', '.')
            bsj_id = row.get('bsj_id', '').split('|')[0]
            # cross-run rows carry BED12 block columns from smart_merge output
            block_count  = row.get('sel_block_count',  '1')
            block_sizes  = row.get('sel_block_sizes',  '') or f"{int(end)-int(start)},"
            block_starts = row.get('sel_block_starts', '0')
            if chrom and start and end:
                fh.write(
                    f"{chrom}\t{start}\t{end}\t{bsj_id}\t0\t{strand}"
                    f"\t{start}\t{end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n"
                )
                spliced[bsj_id] = _spliced_length(block_sizes)

    def run_intersect(b_file, strand_flag):
        """Run bedtools -split -wo; return {bsj_id: total_overlap_bp}."""
        r = subprocess.run(
            f"bedtools intersect -a {tmp_bed} -b {b_file} {strand_flag} -split -wo",
            shell=True, capture_output=True, text=True
        )
        totals = defaultdict(int)
        for line in r.stdout.strip().split('\n'):
            p = line.split('\t')
            # A=12 cols, B=6 cols, overlap=1 → 19 cols total
            if len(p) >= 19 and p[3]:
                try:
                    totals[p[3]] += int(p[18])
                except ValueError:
                    pass
        return totals

    gene_sense_ovlp     = run_intersect(gene_bed, '-s')
    gene_antisense_ovlp = run_intersect(gene_bed, '-S')
    exon_ovlp           = run_intersect(exon_bed,  '-s')

    types = {}
    for bsj_id, slen in spliced.items():
        in_sense     = bsj_id in gene_sense_ovlp
        in_antisense = bsj_id in gene_antisense_ovlp

        if in_sense:
            exon_bp   = exon_ovlp.get(bsj_id, 0)
            exon_frac = exon_bp / slen if slen > 0 else 0.0
            if exon_frac >= 1.0:
                typ = 'eciRNA'
            elif exon_frac > 0:
                typ = 'EIciRNA'
            else:
                typ = 'ciRNA'
        elif in_antisense:
            typ = 'antisense'
        else:
            typ = 'intergenic'

        types[bsj_id] = typ

    try:
        os.remove(tmp_bed)
    except OSError:
        pass

    return types


def main():
    args = parse_args()

    sample_indexes = {name: load_sample_index(path)
                      for name, path in zip(args.sample_names, args.sample_tsvs)}

    # Read raw input TSV once
    with open(args.input_tsv) as fh:
        reader     = csv.DictReader(fh, delimiter='\t')
        in_fields  = list(reader.fieldnames)
        rows       = list(reader)

    # Classify types via bedtools
    types = classify_types(rows, args.gene_bed, args.exon_bed)

    # Per-sample columns to append
    per_sample_cols = []
    for name in args.sample_names:
        per_sample_cols += [f'{name}_bsj_consensus', f'{name}_isoform_consensus']

    full_fields  = in_fields + ['n_samples', 'type'] + per_sample_cols
    clean_fields = CLEAN_COLUMNS + per_sample_cols

    filtered_ids = set()

    with open(f'{args.prefix}_confidence.tsv', 'w', newline='') as fh_full, \
         open(f'{args.prefix}_clean.tsv',      'w', newline='') as fh_clean:

        full_w  = csv.DictWriter(fh_full,  fieldnames=full_fields,  delimiter='\t', extrasaction='ignore')
        clean_w = csv.DictWriter(fh_clean, fieldnames=clean_fields, delimiter='\t', extrasaction='ignore')
        full_w.writeheader()
        clean_w.writeheader()

        for row in rows:
            bsj_id = row.get('bsj_id', '').split('|')[0]
            row['n_samples'] = row.get('bsj_confidence', '')
            row['type']      = types.get(bsj_id, '')
            for name in args.sample_names:
                bc, ic = lookup(sample_indexes[name], bsj_id, args.bsj_tol)
                row[f'{name}_bsj_consensus']     = bc
                row[f'{name}_isoform_consensus'] = ic

            full_w.writerow(row)

            try:
                count = int(row.get('bsj_confidence', 0))
            except ValueError:
                count = 0
            if count >= args.min_count:
                clean_w.writerow(row)
                filtered_ids.add(bsj_id)

    # Write filtered BED12
    with open(args.input_bed) as fi, open(f'{args.prefix}.bed12', 'w') as fo:
        for line in fi:
            cols = line.split('\t')
            if len(cols) >= 4 and cols[3].split('|')[0] in filtered_ids:
                fo.write(line)


if __name__ == '__main__':
    main()
