#!/usr/bin/env python3
"""
circrna_clean.py

Produces a clean, wet-lab-friendly TSV from an annotated circRNA confidence TSV
(output of add_class_codes.py).

Adds two columns that are not in the annotated TSV:
  type        — eciRNA | EIciRNA | ciRNA | antisense | intergenic
                (bedtools intersect against gene/exon BED files derived from GTF)
  expression  — read count, priority: isocirc > ciri-long > circfl > circnick

Clean TSV columns:
  #chrom, start, end, strand,
  sel_block_count, sel_block_sizes, sel_block_starts,
  bsj_id, bsj_confidence, isoform_confidence,
  class_code, ref_gene_id,
  type, expression

Isoform rows (bsj_id containing '|iso') are excluded — clean TSV shows main
circRNAs only.  The full annotated TSV (with isoforms) remains unmodified.

Usage:
    circrna_clean.py \\
        --annotated_tsv  sample_discovery.annotated.tsv \\
        --gene_bed       genes.bed \\
        --exon_bed       exons.bed \\
        --prefix         sample_discovery \\
        [--bsj_tol 5] \\
        [--iso_expr    isocirc.out] \\
        [--ciri_expr   CIRI-long.expression] \\
        [--circfl_expr circfl_pass.txt] \\
        [--nick_expr   circRNA_candidates.annotated.txt]
"""

import argparse
import csv
import os
import subprocess
from collections import defaultdict

CLEAN_COLUMNS = [
    '#chrom', 'start', 'end', 'strand',
    'sel_block_count', 'sel_block_sizes', 'sel_block_starts',
    'bsj_id', 'bsj_confidence', 'isoform_confidence',
    'class_code', 'ref_gene_id',
    'type',
    'supporting_reads',
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotated_tsv', required=True)
    p.add_argument('--gene_bed',      required=True)
    p.add_argument('--exon_bed',      required=True)
    p.add_argument('--prefix',        required=True)
    p.add_argument('--bsj_tol',       type=int, default=5)
    p.add_argument('--iso_expr',    default='', help='isocirc.out')
    p.add_argument('--ciri_expr',   default='', help='CIRI-long.expression')
    p.add_argument('--circfl_expr', default='', help='circfl_pass.txt')
    p.add_argument('--nick_expr',   default='', help='circRNA_candidates.annotated.txt')
    return p.parse_args()


# ── Expression loading ─────────────────────────────────────────────────────────

def _build_index(entries):
    """entries: [(chrom, start, end, strand, count)]. Returns {(chrom,strand):[(s,e,count)]}."""
    idx = defaultdict(list)
    for chrom, start, end, strand, count in entries:
        idx[(chrom, strand)].append((start, end, count))
    return idx


def _lookup(idx, chrom, start, end, strand, tol):
    for s, e, count in idx.get((chrom, strand), []):
        if abs(s - start) <= tol and abs(e - end) <= tol:
            return count
    # fallback: strand-agnostic entries (stored under '.')
    for s, e, count in idx.get((chrom, '.'), []):
        if abs(s - start) <= tol and abs(e - end) <= tol:
            return count
    return None


def load_iso_expr(path):
    """isocirc.out — no header; chrom=col2, start=col3 (0-based BED), end=col4, expr=NF-1."""
    entries = []
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 4:
                continue
            try:
                chrom  = c[1]
                start  = int(c[2])
                end    = int(c[3])
                expr   = int(c[-2]) if len(c) > 4 and c[-2] else 0
                strand = c[5] if len(c) > 5 and c[5] in ('+', '-') else '.'
                entries.append((chrom, start, end, strand, expr))
            except (ValueError, IndexError):
                pass
    return _build_index(entries)


def load_ciri_expr(path):
    """CIRI-long.expression — header row 1; col1=coord (chrom:start-end, 1-based), col2=expr."""
    entries = []
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 2:
                continue
            try:
                coord = c[0]
                rpos  = coord.rfind(':')
                chrom = coord[:rpos]
                se    = coord[rpos+1:].split('-')
                start = int(se[0]) - 1      # CIRI-long 1-based → BED 0-based
                end   = int(se[1])
                expr  = int(c[1]) if c[1] else 0
                entries.append((chrom, start, end, '.', expr))
            except (ValueError, IndexError):
                pass
    return _build_index(entries)


def load_circfl_expr(path):
    """circfl_pass.txt — header row 1; chr=col3, start=col4 (1-based), end=col5, readCount=col17."""
    entries = []
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 17:
                continue
            try:
                chrom = c[2]
                start = int(c[3]) - 1       # circFL 1-based → BED 0-based
                end   = int(c[4])
                expr  = int(c[16]) if c[16] else 0
                entries.append((chrom, start, end, '.', expr))
            except (ValueError, IndexError):
                pass
    return _build_index(entries)


def load_nick_expr(path):
    """circRNA_candidates.annotated.txt — header row 1; chrom=col2, start=col3 (0-based), end=col4, count=col6."""
    entries = []
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 6:
                continue
            try:
                chrom  = c[1]
                start  = int(c[2])
                end    = int(c[3])
                expr   = int(c[5]) if c[5] else 0
                strand = c[4] if len(c) > 4 and c[4] in ('+', '-') else '.'
                entries.append((chrom, start, end, strand, expr))
            except (ValueError, IndexError):
                pass
    return _build_index(entries)


def get_expression(chrom, start, end, strand, tol, *indexes):
    """Return first non-None expression value across indexes (priority order)."""
    for idx in indexes:
        val = _lookup(idx, chrom, start, end, strand, tol)
        if val is not None:
            return val
    return ''


# ── circRNA type classification ────────────────────────────────────────────────

def classify_types(rows, gene_bed, exon_bed):
    """
    Classify each main circRNA as: eciRNA, EIciRNA, ciRNA, antisense, intergenic.
    Returns {bsj_id: type_string}.
    """
    if not rows:
        return {}

    tmp_bed = '.type_class_tmp.bed'
    id_to_coord = {}
    with open(tmp_bed, 'w') as fh:
        for row in rows:
            chrom  = row.get('#chrom', '')
            start  = row.get('start', '')
            end    = row.get('end', '')
            strand = row.get('strand', '.')
            bsj_id = row.get('bsj_id', '').split('|')[0]
            if chrom and start and end:
                fh.write(f"{chrom}\t{start}\t{end}\t{bsj_id}\t0\t{strand}\n")
                id_to_coord[bsj_id] = (chrom, start, end, strand)

    def intersect_ids(flags):
        r = subprocess.run(
            f"bedtools intersect -a {tmp_bed} -b {flags} -wa",
            shell=True, capture_output=True, text=True
        )
        ids = set()
        for line in r.stdout.strip().split('\n'):
            p = line.split('\t')
            if len(p) >= 4 and p[3]:
                ids.add(p[3])
        return ids

    gene_sense     = intersect_ids(f"{gene_bed} -s")
    exon_full      = intersect_ids(f"{exon_bed} -s -f 1.0")
    exon_partial   = intersect_ids(f"{exon_bed} -s")
    gene_antisense = intersect_ids(f"{gene_bed} -S")

    types = {}
    for bsj_id in id_to_coord:
        if bsj_id in gene_sense:
            if bsj_id in exon_full:
                typ = 'eciRNA'
            elif bsj_id in exon_partial:
                typ = 'EIciRNA'
            else:
                typ = 'ciRNA'
        elif bsj_id in gene_antisense:
            typ = 'antisense'
        else:
            typ = 'intergenic'
        types[bsj_id] = typ

    try:
        os.remove(tmp_bed)
    except OSError:
        pass

    return types


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    iso_idx    = load_iso_expr(args.iso_expr)
    ciri_idx   = load_ciri_expr(args.ciri_expr)
    circfl_idx = load_circfl_expr(args.circfl_expr)
    nick_idx   = load_nick_expr(args.nick_expr)

    # Read annotated TSV; skip isoform rows
    with open(args.annotated_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = [row for row in reader if '|iso' not in row.get('bsj_id', '')]

    # Classify types (batch bedtools call)
    types = classify_types(rows, args.gene_bed, args.exon_bed)

    with open(f'{args.prefix}_clean.tsv', 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=CLEAN_COLUMNS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            bsj_id = row.get('bsj_id', '').split('|')[0]
            try:
                start  = int(row.get('start', 0))
                end    = int(row.get('end', 0))
            except ValueError:
                start, end = 0, 0
            chrom  = row.get('#chrom', '')
            strand = row.get('strand', '.')

            row['type']             = types.get(bsj_id, '')
            row['supporting_reads'] = get_expression(
                chrom, start, end, strand, args.bsj_tol,
                iso_idx, ciri_idx, circfl_idx, nick_idx
            )
            writer.writerow(row)


if __name__ == '__main__':
    main()
