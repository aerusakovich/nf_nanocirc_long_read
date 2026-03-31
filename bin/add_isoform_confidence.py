#!/usr/bin/env python3
"""
add_isoform_confidence.py

Adds three confidence components and a final tool_consensus score to the
confidence TSV produced by merge_circrna.py.

Components:
    bsj_score       : percentage of active tools detecting the BSJ → 1-4
                        ≤ 25%  → 1
                        ≤ 50%  → 2
                        ≤ 75%  → 3
                        > 75%  → 4
    isoform_score   : percentage of active tools with isoform support → 1-4
                      (same binning, minimum 1 — a single tool agrees with itself)
    overlap_score   : average pairwise spliced-length overlap fraction → 1-4
                        < 25%  → 1
                        25-50% → 2
                        50-75% → 3
                        ≥ 75%  → 4
                      1 if no pairs (single tool)

Final score = bsj_score + isoform_score + overlap_score  (max 12, min 3)
    3-4   → Low
    5-8   → Medium
    9-12  → High

NOTE: tool_consensus reflects agreement among the tools that were run.
      Running fewer than 4 tools reduces scoring resolution.

Usage:
    python3 add_isoform_confidence.py \\
        --confidence  sample_strict_confidence.tsv \\
        --pairs       isocirc_vs_circfl.txt ... \\
        --output      sample_strict_confidence_final.tsv \\
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


def final_score_to_cat(score):
    """Convert final score (3-12) to tool_consensus category."""
    if score <= 4:
        return "Low"
    elif score <= 8:
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


def compute_scores(bsj_id, bsj_confidence, tool_flags, all_pair_results, n_active):
    """
    Compute all three score components for one circRNA.

    Returns:
        isoform_confidence : int  — number of tools with isoform support
        bsj_score          : int  — percentage-based score (1-4)
        isoform_score      : int  — percentage-based, min 1 (1 tool always agrees with itself)
        overlap_score      : int  — 1-4 based on avg pairwise fraction, 1 if no pairs
        final_score        : int  — sum of three components (3-12)
        tool_consensus     : str  — Low / Medium / High
        pair_fracs         : dict — { (tool_a, tool_b): (frac_a, frac_b) or (None, None) }
    """
    detected = {tool for tool, flag in tool_flags.items() if flag == "1"}

    # ── isoform_confidence ─────────────────────────────────────────────────
    tools_with_isoform = set()
    all_pair_fracs     = []   # collect all individual fracs for avg

    pair_fracs = {}
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

        if bsj_id in overlapping_a or bsj_id in overlapping_b:
            tools_with_isoform.add(tool_a)
            tools_with_isoform.add(tool_b)

        if frac_a is not None:
            all_pair_fracs.extend([frac_a, frac_b])

    isoform_confidence = len(tools_with_isoform)

    # ── component 1: bsj_score ────────────────────────────────────────────
    # percentage of active tools detecting this BSJ, binned 1-4
    bsj_score = count_to_score(bsj_confidence, n_active)

    # ── component 2: isoform_score ────────────────────────────────────────
    # percentage of active tools with isoform support, binned 1-4, minimum 1
    isoform_score = max(1, count_to_score(isoform_confidence, n_active))

    # ── component 3: overlap_score ────────────────────────────────────────
    # 1 if no pairs (single tool or no overlaps found)
    if all_pair_fracs:
        avg_frac      = sum(all_pair_fracs) / len(all_pair_fracs)
        overlap_score = frac_to_overlap_score(avg_frac)
    else:
        overlap_score = 1

    final_score    = bsj_score + isoform_score + overlap_score
    tool_consensus = final_score_to_cat(final_score)

    return (isoform_confidence, bsj_score, isoform_score,
            overlap_score, final_score, tool_consensus, pair_fracs)


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
        header    = None
        tool_cols = []

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

                pair_frac_cols = []
                for tool_a, tool_b in pair_order:
                    pair_frac_cols.append("{}_vs_{}_frac_{}".format(tool_a, tool_b, tool_a))
                    pair_frac_cols.append("{}_vs_{}_frac_{}".format(tool_a, tool_b, tool_b))

                new_header = (line + "\t" +
                    "\t".join(pair_frac_cols +
                              ["bsj_score", "isoform_score", "overlap_score",
                               "final_score", "tool_consensus"]))
                out_lines.append(new_header)
                continue

            if header is None:
                print_error("No header line found.", context="Line", context_str=str(lineno))

            cols           = line.split("\t")
            bsj_id         = cols[4]
            bsj_confidence = int(cols[bsj_conf_idx])
            tool_flags     = {tool: cols[header.index(tool)] for tool in tool_cols}

            (isoform_confidence, bsj_score, isoform_score,
             overlap_score, final_score, tool_consensus, pair_fracs) = compute_scores(
                bsj_id, bsj_confidence, tool_flags, all_pair_results, args.n_active)

            cols[iso_conf_idx] = str(isoform_confidence)

            frac_values = []
            for tool_a, tool_b in pair_order:
                frac_a, frac_b = pair_fracs.get((tool_a, tool_b), (None, None))
                frac_values.append("{:.4f}".format(frac_a) if frac_a is not None else "")
                frac_values.append("{:.4f}".format(frac_b) if frac_b is not None else "")

            out_lines.append(
                "\t".join(cols) + "\t" +
                "\t".join(frac_values + [
                    str(bsj_score), str(isoform_score),
                    str(overlap_score), str(final_score), tool_consensus
                ])
            )

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")

    print("Done: wrote {}".format(args.output))


if __name__ == "__main__":
    main()