#!/usr/bin/env python3
"""
add_isoform_confidence.py

Adds three confidence components and two independent consensus labels to the
confidence TSV produced by merge_circrna.py / smart_merge.py.

Components:
    bsj_score       : percentage of active tools detecting the BSJ → 1-4
                        ≤ 25%  → 1
                        ≤ 50%  → 2
                        ≤ 75%  → 3
                        > 75%  → 4
    isoform_score   : percentage of active tools supporting this exact exon
                      structure → 1-4  (same binning, minimum 1)
    overlap_score   : average pairwise spliced-length overlap fraction → 1-4
                        < 25%  → 1
                        25-50% → 2
                        50-75% → 3
                        ≥ 75%  → 4
                      1 if no pairs (single tool)

Independent consensus labels (score 1→Low, 2-3→Medium, 4→High):
    bsj_consensus     : quality of BSJ detection across tools
    isoform_consensus : quality of exon-structure agreement across tools

These two metrics are intentionally separate — a circRNA can have a
well-supported BSJ but an uncertain exon structure (or vice versa).
filter_confidence.py uses them independently:
  no_low       → removes records with Low on EITHER bsj_consensus OR
                 isoform_consensus
  trusted_only → removes Low bsj_consensus unless from a trusted tool
                 (CIRI-long / IsoCirc); also removes Low isoform_consensus

NOTE: scores reflect agreement among the tools that were run.
      Running fewer than 4 tools reduces scoring resolution.

Usage:
    python3 add_isoform_confidence.py \\
        --confidence  sample_smart_consensus_confidence.tsv \\
        --pairs       isocirc_vs_circfl.txt ... \\
        --output      sample_smart_consensus_confidence_scored.tsv \\
        --min_overlap 0.95 \\
        --n_active    4
"""

import os
import sys
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--confidence",   required=True)
    parser.add_argument("--pairs",        required=True, nargs="+")
    parser.add_argument("--output",       required=True)
    parser.add_argument("--min_overlap",  type=float, default=0.95)
    parser.add_argument("--n_active",     type=int,   default=4,
                        help="Number of active detection tools (used for percentage-based scoring)")
    parser.add_argument("--strip_isoform_suffix", action="store_true", default=False,
                        help="Strip '|iso*' suffix from bsj_id before pair lookups "
                             "(required for smart_merge TSVs where isoforms are labelled "
                             "as chr:start-end:strand|iso1)")
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input files -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input files -> {}\n{}: '{}'".format(
            error, context, context_str)
    print(error_str)
    sys.exit(1)


def parse_tool_names(filename):
    basename = os.path.basename(filename)
    stem     = basename.rsplit(".", 1)[0]
    if "_vs_" not in stem:
        print_error("Pairwise filename must contain '_vs_'.",
                    context="Filename", context_str=basename)
    return stem.split("_vs_", 1)


def make_bsj_key(chrom, start, end, strand):
    return "{}:{}-{}:{}".format(chrom, start, end, strand)


def spliced_length(block_sizes_str):
    return sum(int(x) for x in block_sizes_str.split(",") if x.strip())


def count_to_score(count, n_total):
    """Convert count/n_total ratio to 1-4 score using percentage bins."""
    if n_total == 0:
        return 1
    pct = count / n_total
    if pct <= 0.25:
        return 1
    elif pct <= 0.50:
        return 2
    elif pct <= 0.75:
        return 3
    else:
        return 4


def frac_to_overlap_score(avg_frac):
    """Convert average pairwise fraction to 1-4 score."""
    if avg_frac < 0.25:
        return 1
    elif avg_frac < 0.50:
        return 2
    elif avg_frac < 0.75:
        return 3
    else:
        return 4


def score_to_cat(score):
    """Convert a 1-4 component score to Low / Medium / High."""
    if score == 1:
        return "Low"
    elif score <= 3:
        return "Medium"
    else:
        return "High"


def read_overlapping_pairs(path, tool_a, tool_b, min_overlap):
    """
    Read one bedtools intersect -split -wo output file.

    Returns:
        overlapping_a  : set of bsj_keys from tool_a passing threshold
        overlapping_b  : set of bsj_keys from tool_b passing threshold
        fractions      : { (key_a, key_b): (frac_a, frac_b) } for ALL pairs
    """
    overlapping_a = set()
    overlapping_b = set()
    fractions     = {}

    if not os.path.isfile(path):
        print_error("Pairwise file not found.", context="Path", context_str=path)

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 25:
                print_error(
                    "Expected 25 columns, got {}".format(len(cols)),
                    context="{}_vs_{} line".format(tool_a, tool_b),
                    context_str=str(lineno))

            key_a = make_bsj_key(cols[0],  cols[1],  cols[2],  cols[5])
            key_b = make_bsj_key(cols[12], cols[13], cols[14], cols[17])

            try:
                s_a   = spliced_length(cols[10])
                s_b   = spliced_length(cols[22])
                ovlp  = int(cols[24])
            except (ValueError, IndexError) as e:
                sys.stderr.write(f"WARNING line {lineno}: {e}\n")
                continue

            if s_a == 0 or s_b == 0:
                continue

            frac_a = ovlp / s_a
            frac_b = ovlp / s_b

            fractions[(key_a, key_b)] = (frac_a, frac_b)

            if frac_a >= min_overlap and frac_b >= min_overlap:
                overlapping_a.add(key_a)
                overlapping_b.add(key_b)

    return overlapping_a, overlapping_b, fractions


def compute_scores(bsj_id, bsj_confidence, tool_flags, all_pair_results, n_active,
                   isoform_tools_set=None):
    """
    Compute all three score components for one circRNA.

    isoform_tools_set : optional set of tool names that agree on this specific
                        isoform's exon structure (from smart_merge 'isoform_tools'
                        column).  When provided:
                          - isoform_confidence = len(isoform_tools_set)
                          - overlap fracs are restricted to pairs where BOTH tools
                            are in isoform_tools_set, so minority isoforms (single
                            tool) correctly get overlap_score=1 rather than
                            inheriting the locus-level pairwise fractions.
                        When absent (legacy merge TSVs), the original BSJ-key
                        lookup behaviour is preserved.

    Returns:
        isoform_confidence : int  — number of tools supporting this specific structure
        bsj_score          : int  — BSJ agreement score (1-4)
        bsj_consensus      : str  — Low / Medium / High for BSJ detection
        isoform_score      : int  — isoform structure agreement score (1-4, min 1)
        isoform_consensus  : str  — Low / Medium / High for exon structure
        overlap_score      : int  — 1-4 based on avg pairwise fraction, 1 if no pairs
        pair_fracs         : dict — { (tool_a, tool_b): (frac_a, frac_b) or (None, None) }
    """
    detected = {tool for tool, flag in tool_flags.items() if flag == "1"}

    pair_fracs     = {}
    all_pair_fracs = []

    for (tool_a, tool_b), (overlapping_a, overlapping_b, fractions) in all_pair_results.items():
        if tool_a not in detected or tool_b not in detected:
            pair_fracs[(tool_a, tool_b)] = (None, None)
            continue

        # find fractions for this bsj_id in this pair
        frac_a = frac_b = None
        if (bsj_id, bsj_id) in fractions:
            frac_a, frac_b = fractions[(bsj_id, bsj_id)]
        else:
            for (ka, kb), (fa, fb) in fractions.items():
                if ka == bsj_id or kb == bsj_id:
                    frac_a, frac_b = fa, fb
                    break

        pair_fracs[(tool_a, tool_b)] = (frac_a, frac_b)

        # When isoform_tools_set is given, only count overlap fracs for pairs
        # where both tools support THIS specific isoform structure.
        if isoform_tools_set is not None:
            if tool_a in isoform_tools_set and tool_b in isoform_tools_set:
                if frac_a is not None:
                    all_pair_fracs.extend([frac_a, frac_b])
        else:
            # legacy path: accumulate all fracs as before
            if frac_a is not None:
                all_pair_fracs.extend([frac_a, frac_b])

    # ── isoform_confidence ────────────────────────────────────────────────
    if isoform_tools_set is not None:
        # Use the isoform_tools column directly: it records which tools produced
        # this exact exon structure, so it IS the per-isoform agreement count.
        isoform_confidence = len(isoform_tools_set)
    else:
        # Legacy path: derive from pairwise BSJ-key overlaps
        tools_with_isoform = set()
        for (tool_a, tool_b), (overlapping_a, overlapping_b, _) in all_pair_results.items():
            if bsj_id in overlapping_a or bsj_id in overlapping_b:
                tools_with_isoform.add(tool_a)
                tools_with_isoform.add(tool_b)
        isoform_confidence = len(tools_with_isoform)

    # ── component 1: bsj_score / bsj_consensus ───────────────────────────
    bsj_score     = count_to_score(bsj_confidence, n_active)
    bsj_consensus = score_to_cat(bsj_score)

    # ── component 2: isoform_score / isoform_consensus ───────────────────
    isoform_score     = max(1, count_to_score(isoform_confidence, n_active))
    isoform_consensus = score_to_cat(isoform_score)

    # ── component 3: overlap_score ────────────────────────────────────────
    if all_pair_fracs:
        avg_frac      = sum(all_pair_fracs) / len(all_pair_fracs)
        overlap_score = frac_to_overlap_score(avg_frac)
    else:
        overlap_score = 1

    return (isoform_confidence, bsj_score, bsj_consensus,
            isoform_score, isoform_consensus, overlap_score, pair_fracs)


def main():
    args = parse_args()

    all_pair_results = {}
    pair_order       = []

    for pair_file in args.pairs:
        tool_a, tool_b = parse_tool_names(pair_file)
        overlapping_a, overlapping_b, fractions = read_overlapping_pairs(
            pair_file, tool_a, tool_b, args.min_overlap)
        all_pair_results[(tool_a, tool_b)] = (overlapping_a, overlapping_b, fractions)
        pair_order.append((tool_a, tool_b))
        print("Loaded {}_vs_{}: {} passing keys in {}, {} in {} (threshold={})".format(
            tool_a, tool_b, len(overlapping_a), tool_a,
            len(overlapping_b), tool_b, args.min_overlap))

    out_lines = []

    with open(args.confidence) as fh:
        header           = None
        tool_cols        = []
        iso_tools_idx    = None   # index of 'isoform_tools' column, if present

        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                header       = line.split("\t")
                bsj_conf_idx = header.index("bsj_confidence")
                iso_conf_idx = header.index("isoform_confidence")
                tool_cols    = [c for c in header[bsj_conf_idx + 1 : iso_conf_idx]
                                if not c.endswith("_block_sizes")
                                and not c.endswith("_block_starts")]

                # 'isoform_tools' is present in smart_merge TSVs; use it for
                # per-isoform structure agreement scoring.
                iso_tools_idx = header.index("isoform_tools") if "isoform_tools" in header else None

                pair_frac_cols = []
                for tool_a, tool_b in pair_order:
                    pair_frac_cols.append("{}_vs_{}_frac_{}".format(tool_a, tool_b, tool_a))
                    pair_frac_cols.append("{}_vs_{}_frac_{}".format(tool_a, tool_b, tool_b))

                new_header = (line + "\t" +
                    "\t".join(pair_frac_cols +
                              ["bsj_score", "bsj_consensus",
                               "isoform_score", "isoform_consensus",
                               "overlap_score"]))
                out_lines.append(new_header)
                continue

            if header is None:
                print_error("No header line found.", context="Line", context_str=str(lineno))

            cols           = line.split("\t")
            bsj_id         = cols[4]
            # For smart_merge isoform entries (e.g. chr:100-200:+|iso1), strip
            # the suffix before looking up in pairwise overlap files, which are
            # keyed on plain BSJ coordinates.
            lookup_id      = bsj_id.split('|')[0] if args.strip_isoform_suffix else bsj_id
            bsj_confidence = int(cols[bsj_conf_idx])
            tool_flags     = {tool: cols[header.index(tool)] for tool in tool_cols}

            # When 'isoform_tools' is present, use it for per-isoform scoring so
            # that minority isoforms (fewer supporting tools) score lower than the
            # consensus isoform at the same BSJ.
            isoform_tools_set = None
            if iso_tools_idx is not None:
                raw = cols[iso_tools_idx].strip()
                if raw and raw != '.':
                    isoform_tools_set = set(raw.split(','))

            (isoform_confidence, bsj_score, bsj_consensus,
             isoform_score, isoform_consensus, overlap_score, pair_fracs) = compute_scores(
                lookup_id, bsj_confidence, tool_flags, all_pair_results, args.n_active,
                isoform_tools_set=isoform_tools_set)

            cols[iso_conf_idx] = str(isoform_confidence)

            frac_values = []
            for tool_a, tool_b in pair_order:
                frac_a, frac_b = pair_fracs.get((tool_a, tool_b), (None, None))
                frac_values.append("{:.4f}".format(frac_a) if frac_a is not None else "")
                frac_values.append("{:.4f}".format(frac_b) if frac_b is not None else "")

            out_lines.append(
                "\t".join(cols) + "\t" +
                "\t".join(frac_values + [
                    str(bsj_score), bsj_consensus,
                    str(isoform_score), isoform_consensus,
                    str(overlap_score)
                ])
            )

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")

    print("Done: wrote {}".format(args.output))


if __name__ == "__main__":
    main()