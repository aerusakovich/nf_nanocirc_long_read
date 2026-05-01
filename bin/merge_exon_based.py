#!/usr/bin/env python3
"""
merge_exon_based.py

Groups circRNAs by exon structure similarity using pairwise bedtools
intersect output (run with -split -wo).

Confidence scoring (same system as add_isoform_confidence.py):
    bsj_score       : percentage of active tools in group → 1-4
                        ≤ 25%  → 1 | ≤ 50%  → 2 | ≤ 75%  → 3 | > 75%  → 4
    isoform_score   : percentage of active tools with isoform overlap → 1-4, min 1
    overlap_score   : avg pairwise fraction → 1-4 (1 if no pairs)

Independent consensus labels (score 1→Low, 2-3→Medium, 4→High):
    bsj_consensus     : quality of BSJ detection across tools
    isoform_consensus : quality of exon-structure agreement across tools

NOTE: scores reflect agreement among the tools that were run.
      Running fewer than 4 tools reduces scoring resolution.

Usage:
    python3 merge_exon_based.py \\
        --pairs   isocirc_vs_circfl.txt ... \\
        --sample  SAMPLE_ID \\
        --outdir  results/ \\
        --min_overlap 0.95 \\
        --n_active    4
"""

import os
import sys
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs",        required=True, nargs="+")
    parser.add_argument("--sample",       required=True)
    parser.add_argument("--outdir",       default=".")
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


def find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(parent, x, y):
    root_x = find(parent, x)
    root_y = find(parent, y)
    if root_x != root_y:
        parent[root_y] = root_x


def make_bsj_key(chrom, start, end, strand):
    return "{}:{}-{}:{}".format(chrom, start, end, strand)


def spliced_length(block_sizes_str):
    return sum(int(x) for x in block_sizes_str.split(",") if x.strip())


def count_to_score(count, n_total):
    """Convert count/n_total ratio to 1-4 score using percentage bins."""
    if n_total == 0:
        return 1
    pct = count / n_total
    if pct <= 0.25: return 1
    elif pct <= 0.50: return 2
    elif pct <= 0.75: return 3
    else: return 4


def frac_to_overlap_score(avg_frac):
    if avg_frac < 0.25: return 1
    elif avg_frac < 0.50: return 2
    elif avg_frac < 0.75: return 3
    else: return 4


def score_to_cat(score):
    """Convert a 1-4 component score to Low / Medium / High."""
    if score == 1: return "Low"
    elif score <= 3: return "Medium"
    else: return "High"


def read_all_records(pair_files, min_overlap):
    all_nodes         = set()
    overlapping_pairs = []
    active_tools      = []
    # store fracs per node pair for scoring
    node_pair_fracs   = {}  # (node_a, node_b): (frac_a, frac_b)
    nodes_with_isoform = set()

    for pair_file in pair_files:
        tool_a, tool_b = parse_tool_names(pair_file)
        if tool_a not in active_tools: active_tools.append(tool_a)
        if tool_b not in active_tools: active_tools.append(tool_b)

        with open(pair_file) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line: continue

                cols = line.split("\t")
                if len(cols) < 25:
                    print_error("Expected 25 columns, got {}".format(len(cols)),
                                context="{}_vs_{} line".format(tool_a, tool_b),
                                context_str=str(lineno))

                chrom_a  = cols[0];  start_a = cols[1];  end_a = cols[2];  strand_a = cols[5]
                chrom_b  = cols[12]; start_b = cols[13]; end_b = cols[14]; strand_b = cols[17]

                try:
                    s_a  = spliced_length(cols[10])
                    s_b  = spliced_length(cols[22])
                    ovlp = int(cols[24])
                except (ValueError, IndexError) as e:
                    sys.stderr.write(f"WARNING line {lineno}: {e}\n")
                    continue

                if s_a == 0 or s_b == 0: continue

                frac_a = ovlp / s_a
                frac_b = ovlp / s_b

                node_a = (tool_a, chrom_a, start_a, end_a, strand_a)
                node_b = (tool_b, chrom_b, start_b, end_b, strand_b)

                all_nodes.add(node_a)
                all_nodes.add(node_b)
                node_pair_fracs[(node_a, node_b)] = (frac_a, frac_b)

                if frac_a >= min_overlap and frac_b >= min_overlap:
                    overlapping_pairs.append((node_a, node_b))
                    nodes_with_isoform.add(node_a)
                    nodes_with_isoform.add(node_b)

    return all_nodes, overlapping_pairs, active_tools, nodes_with_isoform, node_pair_fracs


def build_groups(all_nodes, overlapping_pairs):
    parent = {node: node for node in all_nodes}
    for node_a, node_b in overlapping_pairs:
        union(parent, node_a, node_b)

    root_to_nodes = {}
    for node in all_nodes:
        root = find(parent, node)
        if root not in root_to_nodes: root_to_nodes[root] = []
        root_to_nodes[root].append(node)

    groups = []
    for root, nodes in root_to_nodes.items():
        group = {}
        for tool_name, chrom, start, end, strand in nodes:
            if tool_name not in group: group[tool_name] = []
            group[tool_name].append((chrom, start, end, strand))
        groups.append(group)
    return groups


def compute_group_scores(group, nodes_with_isoform, node_pair_fracs, n_active):
    """
    Compute all score components for a group.
    """
    tools_in_group = set(group.keys())
    # component 1: percentage of active tools detecting this group, binned 1-4
    bsj_score = count_to_score(len(tools_in_group), n_active)

    # component 2: percentage of active tools with isoform support — min 1
    tools_with_isoform = set()
    for tool_name, coord_list in group.items():
        for coords in coord_list:
            node = (tool_name,) + coords
            if node in nodes_with_isoform:
                tools_with_isoform.add(tool_name)
    isoform_confidence = len(tools_with_isoform)
    isoform_score      = max(1, count_to_score(isoform_confidence, n_active))

    # component 3: avg pairwise fraction — 1 if no pairs
    all_fracs = []
    for tool_name, coord_list in group.items():
        for coords in coord_list:
            node_a = (tool_name,) + coords
            for other_tool, other_coords_list in group.items():
                if other_tool == tool_name: continue
                for other_coords in other_coords_list:
                    node_b = (other_tool,) + other_coords
                    if (node_a, node_b) in node_pair_fracs:
                        fa, fb = node_pair_fracs[(node_a, node_b)]
                        all_fracs.extend([fa, fb])
                    elif (node_b, node_a) in node_pair_fracs:
                        fa, fb = node_pair_fracs[(node_b, node_a)]
                        all_fracs.extend([fa, fb])

    if all_fracs:
        overlap_score = frac_to_overlap_score(sum(all_fracs) / len(all_fracs))
    else:
        overlap_score = 1

    bsj_consensus     = score_to_cat(bsj_score)
    isoform_consensus = score_to_cat(isoform_score)

    return isoform_confidence, bsj_score, bsj_consensus, isoform_score, isoform_consensus, overlap_score


def write_outputs(groups, active_tools, nodes_with_isoform,
                  node_pair_fracs, sample, outdir, n_active):
    union_lines      = []
    inter_lines      = []
    union_conf_lines = []
    inter_conf_lines = []

    header_cols = (
        ["#chrom", "start", "end", "strand", "bsj_id"]
        + active_tools
        + ["bsj_score", "bsj_consensus",
           "isoform_confidence", "isoform_score", "isoform_consensus",
           "overlap_score"]
    )
    header_line = "\t".join(header_cols)
    union_conf_lines.append(header_line)
    inter_conf_lines.append(header_line)

    for group in sorted(groups, key=lambda g: sorted(
        (tool, coords[0])
        for tool, coord_list in g.items()
        for coords in coord_list
    )):
        tools_in_group = set(group.keys())

        rep_coords = None
        for tool in active_tools:
            if tool in group:
                rep_coords = group[tool][0]
                break
        if rep_coords is None: continue

        chrom, start, end, strand = rep_coords
        try:
            start_int = int(start)
            end_int   = int(end)
        except ValueError:
            print_error("Non-numeric coordinates.", context="coords",
                        context_str=str(rep_coords))

        (isoform_confidence, bsj_score, bsj_consensus,
         isoform_score, isoform_consensus, overlap_score) = compute_group_scores(
            group, nodes_with_isoform, node_pair_fracs, n_active)

        bsj_id     = make_bsj_key(chrom, start, end, strand)
        block_size = end_int - start_int

        bed12_line = "\t".join([
            chrom, str(start_int), str(end_int),
            bsj_id, str(bsj_score), strand,
            str(start_int), str(start_int), "0",
            "1", str(block_size) + ",", "0,",
        ])

        tool_flags = ["1" if t in tools_in_group else "0" for t in active_tools]
        conf_row   = "\t".join(
            [chrom, str(start_int), str(end_int), strand, bsj_id]
            + tool_flags
            + [str(bsj_score), bsj_consensus,
               str(isoform_confidence), str(isoform_score), isoform_consensus,
               str(overlap_score)]
        )

        union_lines.append(bed12_line)
        union_conf_lines.append(conf_row)

        if tools_in_group == set(active_tools):
            inter_lines.append(bed12_line)
            inter_conf_lines.append(conf_row)

    def _write(lines, path):
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    _write(union_lines,      os.path.join(outdir, "{}_exon_union.bed12".format(sample)))
    _write(union_conf_lines, os.path.join(outdir, "{}_exon_union_confidence.tsv".format(sample)))
    _write(inter_lines,      os.path.join(outdir, "{}_exon_intersection.bed12".format(sample)))
    _write(inter_conf_lines, os.path.join(outdir, "{}_exon_intersection_confidence.tsv".format(sample)))

    print("[{}] exon mode: {} groups in union, {} in intersection".format(
        sample, len(union_lines), len(inter_lines)))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    all_nodes, overlapping_pairs, active_tools, nodes_with_isoform, node_pair_fracs = \
        read_all_records(args.pairs, args.min_overlap)

    print("[{}] Tools seen: {}".format(args.sample, active_tools))
    print("[{}] Total nodes: {}, pairs passing threshold: {}".format(
        args.sample, len(all_nodes), len(overlapping_pairs)))

    groups = build_groups(all_nodes, overlapping_pairs)
    print("[{}] Union-find produced {} groups".format(args.sample, len(groups)))

    write_outputs(groups, active_tools, nodes_with_isoform,
                  node_pair_fracs, args.sample, args.outdir, args.n_active)


if __name__ == "__main__":
    main()