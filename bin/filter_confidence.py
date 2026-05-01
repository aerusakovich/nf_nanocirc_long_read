#!/usr/bin/env python3
"""
filter_confidence.py
--------------------
Filters a scored confidence TSV + matching BED12 using the two independent
confidence axes produced by add_isoform_confidence.py:

    bsj_consensus     — quality of BSJ detection across tools (Low/Medium/High)
    isoform_consensus — quality of exon-structure agreement across tools

Filter modes:
  trusted_only  — Balanced mode.
                  Removes Low entries on each axis UNLESS the call comes from
                  a trusted tool (CIRI-long or IsoCirc), which are kept even
                  when confidence is Low.
                  - bsj_consensus Low → keep if bsj_source in {cirilong, isocirc}
                  - isoform_consensus Low → keep if any tool in isoform_tools
                    is in {cirilong, isocirc}
                  CircFL and CircNick-LRS Low calls are removed on both axes.

  no_low        — High-confidence mode.
                  Removes ALL entries with Low on EITHER axis, no exceptions.

  high_only     — Strictest filter.
                  Keeps ONLY entries where BOTH axes are 'High'. Medium is removed.

Usage:
    filter_confidence.py \\
        --bed    sample_smart_consensus.bed12 \\
        --tsv    sample_smart_consensus_confidence.tsv \\
        --mode   trusted_only \\
        --prefix sample_smart_consensus_filtered
"""

import argparse
import sys
import os

TRUSTED_TOOLS = {'cirilong', 'isocirc'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bed',    required=True, help='BED12 input file')
    p.add_argument('--tsv',    required=True, help='Scored confidence TSV input file')
    p.add_argument('--mode',   required=True, choices=['trusted_only', 'no_low', 'high_only'],
                   help='Filter mode')
    p.add_argument('--prefix', required=True,
                   help='Output file prefix (outputs: prefix.bed12, prefix_confidence.tsv)')
    return p.parse_args()


def _bsj_low_passes(cols, idx, mode):
    """
    Decide whether a Low bsj_consensus record is kept.
    no_low → always False.
    trusted_only → True only if the BSJ source is a trusted tool.
    """
    if mode == 'no_low':
        return False
    # trusted_only: check bsj_source (smart TSVs) or any trusted tool flag
    if 'bsj_source' in idx:
        return cols[idx['bsj_source']] in TRUSTED_TOOLS
    return any(t in idx and cols[idx[t]] == '1' for t in TRUSTED_TOOLS)


def _isoform_low_passes(cols, idx, mode):
    """
    Decide whether a Low isoform_consensus record is kept.
    no_low → always False.
    trusted_only → True only if isoform_tools contains a trusted tool.
    """
    if mode == 'no_low':
        return False
    # trusted_only: check isoform_tools column (smart TSVs)
    if 'isoform_tools' in idx:
        tools = set(cols[idx['isoform_tools']].split(','))
        return bool(tools & TRUSTED_TOOLS)
    # fallback for non-smart TSVs: any trusted tool flag
    return any(t in idx and cols[idx[t]] == '1' for t in TRUSTED_TOOLS)


def passes_filter(cols, idx, mode):
    """
    Return True if this TSV row should be kept.
    idx: dict mapping column name → column index.

    Both axes must independently pass:
      bsj_consensus      != Low  (or Low + trusted tool exception)
      isoform_consensus  != Low  (or Low + trusted tool exception)

    high_only: both axes must be exactly 'High' — Medium is also removed.
    """
    bsj_cons = cols[idx['bsj_consensus']]    if 'bsj_consensus'     in idx else 'NA'
    iso_cons = cols[idx['isoform_consensus']] if 'isoform_consensus' in idx else 'NA'

    if mode == 'high_only':
        return bsj_cons == 'High' and iso_cons == 'High'

    bsj_ok = (bsj_cons != 'Low') or _bsj_low_passes(cols, idx, mode)
    iso_ok = (iso_cons != 'Low') or _isoform_low_passes(cols, idx, mode)

    return bsj_ok and iso_ok


def main():
    args = parse_args()

    kept_ids     = set()
    tsv_out      = []
    header_idx   = {}

    # ── Pass 1: filter TSV ────────────────────────────────────────────────────
    with open(args.tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue

            # Header detection: starts with '#' or is the first non-empty line
            if line.startswith('#') or not header_idx:
                header = line.lstrip('#').split('\t')
                # Re-attach '#' to first column for output fidelity
                raw_header = line
                header_idx = {col: i for i, col in enumerate(header)}
                # Also index the raw '#chrom' as 'chrom'
                if line.startswith('#chrom'):
                    header_idx['chrom'] = 0
                tsv_out.append(raw_header)
                continue

            if 'bsj_consensus' not in header_idx and 'isoform_consensus' not in header_idx:
                print(
                    '[filter_confidence] WARNING: bsj_consensus / isoform_consensus '
                    'columns not found — TSV may not have been scored by '
                    'add_isoform_confidence.py. Passing all rows through.',
                    file=sys.stderr
                )
                tsv_out.append(line)
                continue

            cols   = line.split('\t')
            bsj_id = cols[4]   # bsj_id is always column 4 in all confidence TSVs

            if passes_filter(cols, header_idx, args.mode):
                kept_ids.add(bsj_id)
                tsv_out.append(line)

    # ── Pass 2: filter BED12 by kept bsj_ids ─────────────────────────────────
    bed_out = []
    with open(args.bed) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.split('\t')
            if cols[3] in kept_ids:   # BED12 name field (col 4, 0-indexed 3) = bsj_id
                bed_out.append(line)

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(f'{args.prefix}.bed12', 'w') as fh:
        if bed_out:
            fh.write('\n'.join(bed_out) + '\n')

    with open(f'{args.prefix}_confidence.tsv', 'w') as fh:
        if tsv_out:
            fh.write('\n'.join(tsv_out) + '\n')

    print(
        f'[filter_confidence] {args.mode}: '
        f'{len(kept_ids)} entries kept → {args.prefix}.bed12',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
