#!/usr/bin/env python3
"""
merge_circrna.py

Merges BED12 outputs from up to 4 circRNA detection tools.
Produces union, intersection, and confidence tables at two
stringency levels: strict (exact BSJ match) and relaxed (±N bp).

Usage:
    python3 merge_circrna.py \\
        --sample SAMPLE_ID \\
        --isocirc  isocirc.bed \\
        --circfl   circFL_final.bed \\
        --cirilong cirilong.bed12 \\
        --circnick circnick.bed12 \\
        --outdir   results/
"""

import os
import sys
import argparse


def parse_args(args=None):
    Description = "Merge circRNA tool BED12 outputs into union/intersection/confidence tables."
    Epilog = "Example usage: merge_circrna.py --sample S1 --isocirc a.bed --circfl b.bed --outdir out/"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument("--sample",   required=True,  help="Sample name used in output filenames.")
    parser.add_argument("--isocirc",  default=None,   help="BED12 from isocirc (optional).")
    parser.add_argument("--circfl",   default=None,   help="BED12 from circFL-seq (optional).")
    parser.add_argument("--cirilong", default=None,   help="BED12 from CIRI-long (optional).")
    parser.add_argument("--circnick", default=None,   help="BED12 from circnick-lrs (optional).")
    parser.add_argument("--outdir",   default=".",    help="Output directory (default: current dir).")
    parser.add_argument("--tolerance", type=int, default=5,
                        help="BSJ coordinate tolerance in bp for relaxed mode (default: 5).")

    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input files -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input files -> {}\n{}: '{}'".format(
            error, context, context_str
        )
    print(error_str)
    sys.exit(1)


def read_bed12(path, tool_name):
    """
    Read a BED12 file and return a list of dicts.
    """
    records = []

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")

            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            cols = line.split("\t")
            if len(cols) < 12:
                print_error(
                    "Expected 12 columns, got {}".format(len(cols)),
                    context="{} line".format(tool_name),
                    context_str=str(lineno)
                )

            records.append({
                "chrom":        cols[0],
                "start":        int(cols[1]),
                "end":          int(cols[2]),
                "name":         cols[3],
                "score":        cols[4],
                "strand":       cols[5],
                "thick_start":  cols[6],
                "thick_end":    cols[7],
                "rgb":          cols[8],
                "block_count":  cols[9],
                "block_sizes":  cols[10],
                "block_starts": cols[11],
                "tool":         tool_name,
            })

    return records


def group_strict(records):
    """
    Group records by exact BSJ: (chrom, start, end, strand).
    """
    groups = {}

    for rec in records:
        key = (rec["chrom"], rec["start"], rec["end"], rec["strand"])
        if key not in groups:
            groups[key] = {}
        tool = rec["tool"]
        if tool not in groups[key]:
            groups[key][tool] = []
        groups[key][tool].append(rec)

    return groups


def group_relaxed(records, tolerance):
    """
    Group records by BSJ with coordinate tolerance.
    """
    canonical_keys = []
    groups = {}

    for rec in records:
        chrom  = rec["chrom"]
        start  = rec["start"]
        end    = rec["end"]
        strand = rec["strand"]
        tool   = rec["tool"]

        matched_key = None
        for key in canonical_keys:
            k_chrom, k_start, k_end, k_strand = key
            if (k_chrom  == chrom
                and k_strand == strand
                and abs(k_start - start) <= tolerance
                and abs(k_end   - end)   <= tolerance):
                matched_key = key
                break

        if matched_key is None:
            matched_key = (chrom, start, end, strand)
            canonical_keys.append(matched_key)
            groups[matched_key] = {}

        if tool not in groups[matched_key]:
            groups[matched_key][tool] = []
        groups[matched_key][tool].append(rec)

    return groups


def best_record(records):
    """
    From a list of records (same tool, same BSJ), return the one
    with the highest numeric score (read count).
    """
    def _score(rec):
        try:
            return int(rec["score"])
        except ValueError:
            return 0

    return max(records, key=_score)


def write_outputs(groups, active_tools, sample, mode, outdir):
    """
    Write union BED12, intersection BED12, and two confidence TSVs
    (one for union BSJs, one for intersection BSJs only).

    Confidence TSV includes per-tool block_sizes and block_starts columns
    so isoform structure can be compared directly.
    """
    union_lines      = []
    inter_lines      = []
    union_conf_lines = []
    inter_conf_lines = []

    # TSV header — two block columns per tool, then isoform_confidence
    tool_flag_cols  = active_tools
    tool_block_cols = []
    for t in active_tools:
        tool_block_cols.append("{}_block_sizes".format(t))
        tool_block_cols.append("{}_block_starts".format(t))

    header_cols = (
        ["#chrom", "start", "end", "strand", "bsj_id", "bsj_confidence"]
        + tool_flag_cols
        + tool_block_cols
        + ["isoform_confidence"]
    )
    header_line = "\t".join(header_cols)
    union_conf_lines.append(header_line)
    inter_conf_lines.append(header_line)

    for bsj_key, tool_map in sorted(groups.items()):
        chrom, start, end, strand = bsj_key
        bsj_id         = "{}:{}-{}:{}".format(chrom, start, end, strand)
        bsj_confidence = len(tool_map)

        # ── pick best representative record across all tools ──
        all_recs_for_bsj = [
            best_record(recs)
            for recs in tool_map.values()
        ]
        rep = max(all_recs_for_bsj, key=lambda r: int(r["score"]) if r["score"].isdigit() else 0)

        bed12_line = "\t".join([
            chrom,
            str(start),
            str(end),
            bsj_id,
            str(bsj_confidence),
            strand,
            str(start),         # thickStart = thickEnd = start → no CDS
            str(start),         # thickEnd
            "0",                # rgb
            rep["block_count"],
            rep["block_sizes"],
            rep["block_starts"],
        ])

        # ── per-tool presence flags ───────────────────────────
        tool_flags = ["1" if t in tool_map else "0" for t in active_tools]

        # ── per-tool block columns ────────────────────────────
        tool_blocks = []
        for t in active_tools:
            if t in tool_map:
                rec = best_record(tool_map[t])
                tool_blocks.append(rec["block_sizes"])
                tool_blocks.append(rec["block_starts"])
            else:
                tool_blocks.append("")   # block_sizes empty
                tool_blocks.append("")   # block_starts empty

        conf_row = "\t".join(
            [chrom, str(start), str(end), strand, bsj_id, str(bsj_confidence)]
            + tool_flags
            + tool_blocks
            + [""]   # isoform_confidence placeholder
        )

        # ── union: every BSJ ──────────────────────────────────
        union_lines.append(bed12_line)
        union_conf_lines.append(conf_row)

        # ── intersection: only BSJs found by ALL active tools ─
        if set(tool_map.keys()) == set(active_tools):
            inter_lines.append(bed12_line)
            inter_conf_lines.append(conf_row)

    # ── write files ───────────────────────────────────────────
    def _write(lines, path):
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    union_bed_path  = os.path.join(outdir, "{}_{}_union.bed12".format(sample, mode))
    union_conf_path = os.path.join(outdir, "{}_{}_union_confidence.tsv".format(sample, mode))
    inter_bed_path  = os.path.join(outdir, "{}_{}_intersection.bed12".format(sample, mode))
    inter_conf_path = os.path.join(outdir, "{}_{}_intersection_confidence.tsv".format(sample, mode))

    _write(union_lines,      union_bed_path)
    _write(union_conf_lines, union_conf_path)
    _write(inter_lines,      inter_bed_path)
    _write(inter_conf_lines, inter_conf_path)

    print("[{}] {} mode: {} BSJs in union, {} in intersection".format(
        sample, mode, len(union_lines), len(inter_lines)
    ))


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tool_files = {
        "isocirc":  args.isocirc,
        "circfl":   args.circfl,
        "cirilong": args.cirilong,
        "circnick": args.circnick,
    }
    active_tools = {
        name: path
        for name, path in tool_files.items()
        if path is not None and os.path.isfile(path)
    }

    if len(active_tools) < 1:
        print_error("At least 1 tool BED12 file must be provided.")

    print("[{}] Tools provided: {}".format(args.sample, list(active_tools.keys())))

    all_records = []
    for tool_name, path in active_tools.items():
        recs = read_bed12(path, tool_name)
        print("[{}] {}: {} records loaded".format(args.sample, tool_name, len(recs)))
        all_records.extend(recs)

    strict_groups = group_strict(all_records)
    write_outputs(strict_groups, list(active_tools.keys()), args.sample, "strict", args.outdir)

    relaxed_groups = group_relaxed(all_records, args.tolerance)
    write_outputs(relaxed_groups, list(active_tools.keys()), args.sample, "relaxed", args.outdir)


if __name__ == "__main__":
    main()