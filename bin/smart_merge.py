#!/usr/bin/env python3
"""
smart_merge.py
--------------
Produces consensus-based merged BED12 outputs from individual tool
detections, grouped by relaxed BSJ (±tolerance bp, strand-aware).

Four selection modes:

  smart_consensus       — BSJ by majority vote (ties → BSJ priority),
                          structure by majority vote among tools sharing
                          the exact winning BSJ (ties → struct priority)

  smart_consensus_xstruct — same BSJ vote as smart_consensus, but structure
                          voting uses ALL tools in the relaxed group,
                          comparing exon positions in absolute genomic
                          coordinates (±struct_tolerance bp per boundary).
                          Tools with a slightly different BSJ can therefore
                          cast a vote for or against the winning BSJ's
                          structure.  When the cross-BSJ winning structure
                          has no representative at the winning BSJ coords,
                          it is rebased (block_starts adjusted by the small
                          BSJ offset, ≤ struct_tolerance bp).

  smart_consensus_hybrid — same BSJ vote + absolute-coordinate structure vote
                          restricted to exact-BSJ tools only (no rebasing).
                          Minority-BSJ tools emit isoforms at their own coords.
                          Produces the pipeline's 'discovery' output.

  smart_priority        — BSJ unconditionally from BSJ priority tool,
                          structure unconditionally from struct priority tool

Within each relaxed-BSJ group, ALL distinct (BSJ / absolute-structure)
combinations are emitted: the winning combination is labelled 'main';
minority structures are emitted as separate isoform rows ('iso1', 'iso2', …).
The BED12 name field encodes the isoform: 'chr:start-end:strand' for main,
'chr:start-end:strand|iso1' for minority isoforms.

BSJ priority order:
    CIRI-long > circFL > IsoCirc > CircNick-LRS

Exon structure priority order:
    IsoCirc > CIRI-long > CircNick-LRS > circFL

Usage:
    smart_merge.py \\
        --sample   SAMPLE_ID \\
        --tool_names  cirilong isocirc circfl \\
        --bed_files   a.bed    b.bed   c.bed  \\
        --tolerance   5 \\
        --n_active    3 \\
        --outdir      results/
"""

import os
import sys
import argparse
from collections import defaultdict


# ── Priority tables ────────────────────────────────────────────────────────────
BSJ_PRIORITY    = ['cirilong', 'circfl', 'isocirc', 'circnick']
STRUCT_PRIORITY = ['isocirc', 'cirilong', 'circnick', 'circfl']
TRUSTED_TOOLS   = {'cirilong', 'isocirc'}


def parse_args(args=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',     required=True)
    p.add_argument('--tool_names', required=True, nargs='+',
                   help='Tool names in the same order as --bed_files')
    p.add_argument('--bed_files',  required=True, nargs='+',
                   help='BED12 paths in the same order as --tool_names')
    p.add_argument('--tolerance',       type=int, default=5,
                   help='BSJ coordinate tolerance in bp (default: 5)')
    p.add_argument('--struct_tolerance', type=int, default=None,
                   help='Per-exon-boundary tolerance for cross-BSJ structure '
                        'voting in smart_consensus_xstruct (bp). '
                        'Defaults to --tolerance if not set.')
    p.add_argument('--n_active',   type=int, default=None,
                   help='Total active tools (used only for confidence info)')
    p.add_argument('--outdir',     default='.',
                   help='Output directory (default: current dir)')
    return p.parse_args(args)


def print_error(msg):
    print(f'ERROR: {msg}', file=sys.stderr)
    sys.exit(1)


# ── I/O ────────────────────────────────────────────────────────────────────────

def read_bed12(path, tool_name):
    """Read BED12 and return list of record dicts."""
    records = []
    if not path or not os.path.exists(path):
        return records
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.split('\t')
            if len(cols) < 12:
                print_error(f'{tool_name} line {lineno}: expected 12 columns, got {len(cols)}')
            records.append({
                'chrom':        cols[0],
                'start':        int(cols[1]),
                'end':          int(cols[2]),
                'name':         cols[3],
                'score':        cols[4],
                'strand':       cols[5] if cols[5] in ('+', '-', '.') else '.',
                'thick_start':  cols[6],
                'thick_end':    cols[7],
                'rgb':          cols[8],
                'block_count':  cols[9],
                'block_sizes':  cols[10],
                'block_starts': cols[11],
                'tool':         tool_name,
            })
    return records


# ── Grouping ───────────────────────────────────────────────────────────────────

def group_relaxed(records, tolerance):
    """
    Group records by relaxed BSJ: (chrom, strand) + |Δstart| ≤ tol + |Δend| ≤ tol.
    Canonical key is the first-seen BSJ within each tolerance window.
    Returns {(chrom,start,end,strand): {tool: [record, ...]}}
    """
    canonical_keys = []
    groups = {}
    for rec in records:
        chrom, start, end, strand, tool = (
            rec['chrom'], rec['start'], rec['end'], rec['strand'], rec['tool']
        )
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
        groups[matched_key].setdefault(tool, []).append(rec)
    return groups


# ── Helpers ────────────────────────────────────────────────────────────────────

def best_record(recs):
    """Highest-score record from a list (read count proxy)."""
    def _score(r):
        try:
            return int(r['score'])
        except ValueError:
            return 0
    return max(recs, key=_score)


def bsj_key(rec):
    return (rec['start'], rec['end'])


def struct_key(rec):
    return (rec['block_count'], rec['block_sizes'], rec['block_starts'])


def make_bsj_id(chrom, start, end, strand, suffix=None):
    base = f'{chrom}:{start}-{end}:{strand}'
    return f'{base}|{suffix}' if suffix else base


def max_score_of(recs_dict):
    """Largest integer score across a {tool: record} dict."""
    best = 0
    for rec in recs_dict.values():
        try:
            best = max(best, int(rec['score']))
        except ValueError:
            pass
    return best


def _priority_tool(tools_present, priority_list):
    """Return highest-priority tool from a collection, or first if none listed."""
    for t in priority_list:
        if t in tools_present:
            return t
    return next(iter(tools_present))


# ── Absolute-coordinate structure helpers ──────────────────────────────────────

def absolute_exon_coords(rec):
    """
    Convert BED12 (relative) block_starts to absolute (genomic_start, size)
    tuples for each exon.
    """
    base   = rec['start']
    sizes  = [int(x) for x in rec['block_sizes'].split(',')  if x.strip()]
    starts = [int(x) for x in rec['block_starts'].split(',') if x.strip()]
    return tuple((base + s, sz) for s, sz in zip(starts, sizes))


def abs_struct_similar(a, b, tolerance):
    """
    Two absolute structures are similar if they have the same exon count,
    each corresponding exon START is within tolerance bp, and exon SIZES
    are identical.
    """
    if len(a) != len(b):
        return False
    return all(abs(ea[0] - eb[0]) <= tolerance and ea[1] == eb[1]
               for ea, eb in zip(a, b))


def group_by_abs_struct(tool_best, tolerance):
    """
    Group all tools by absolute exon coordinate similarity.
    Returns list of [abs_struct_representative, [tool, ...]].
    First tool seen becomes the representative for each group.
    """
    groups = []  # [ [abs_struct, [tools]] ]
    for tool, rec in tool_best.items():
        abs_struct = absolute_exon_coords(rec)
        matched = False
        for grp in groups:
            if abs_struct_similar(abs_struct, grp[0], tolerance):
                grp[1].append(tool)
                matched = True
                break
        if not matched:
            groups.append([abs_struct, [tool]])
    return groups


def vote_struct_groups(groups):
    """
    Vote over a list of [abs_struct_rep, [tools]].
    Plurality wins; ties broken by STRUCT_PRIORITY.
    Returns (winning_index, agree_count).
    """
    max_count = max(len(tools) for _, tools in groups)
    winners   = [i for i, (_, tools) in enumerate(groups)
                 if len(tools) == max_count]

    if len(winners) == 1:
        return winners[0], max_count

    for prio_tool in STRUCT_PRIORITY:
        for idx in winners:
            if prio_tool in groups[idx][1]:
                return idx, max_count

    return winners[0], max_count


def rebase_struct(rec, new_start, new_end):
    """
    Return a copy of rec with block_starts rebased so the record's chromStart
    becomes new_start.  Absolute exon positions are preserved; only the
    relative offsets change.  Used when the structure-vote winner has a
    slightly different BSJ than the winning BSJ (offset ≤ struct_tolerance).
    """
    offset     = rec['start'] - new_start        # positive if rec is to the right
    old_starts = [int(x) for x in rec['block_starts'].split(',') if x.strip()]
    new_starts = [s + offset for s in old_starts]
    new_rec    = dict(rec)
    new_rec['start']        = new_start
    new_rec['end']          = new_end
    new_rec['block_starts'] = ','.join(str(s) for s in new_starts)
    return new_rec


# ── Voting ─────────────────────────────────────────────────────────────────────

def vote_majority(votes_dict, priority_list):
    """
    Given {key: [tool, ...]} vote tallies, return (winning_key, agree_count).
    Strict majority (>N/2) wins; ties or all-disagree fall to priority_list.
    """
    n = sum(len(tools) for tools in votes_dict.values())
    max_count = max(len(tools) for tools in votes_dict.values())

    if max_count > n / 2:
        winners = [k for k, tools in votes_dict.items() if len(tools) == max_count]
        if len(winners) == 1:
            return winners[0], max_count
        # Tie among top-count keys — use priority tool
        for tool in priority_list:
            for k in winners:
                if tool in votes_dict[k]:
                    return k, max_count

    # No majority — fall to priority list
    all_tools_present = {t for tools in votes_dict.values() for t in tools}
    for tool in priority_list:
        if tool in all_tools_present:
            for k, tools in votes_dict.items():
                if tool in tools:
                    return k, len(tools)

    k = next(iter(votes_dict))
    return k, len(votes_dict[k])


# ── Entry collection ───────────────────────────────────────────────────────────

def collect_entries_consensus(tool_best):
    """
    Collect all isoform entries for smart_consensus mode.

    For each relaxed group:
      1. Majority vote on (start, end)  → winning BSJ
      2. Among tools agreeing on BSJ: majority vote on structure → main entry
         Each minority structure at winning BSJ → iso entry
      3. Tools with a different BSJ (within tolerance) → iso entry per distinct BSJ
         (priority structure chosen for each)

    Returns list of entry dicts:
      start, end, struct_rec, bsj_src, struct_src,
      bsj_agree, struct_agree, isoform_label, isoform_tools
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)

    agree_map    = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}
    disagree_bsj = defaultdict(dict)
    for t, r in tool_best.items():
        if bsj_key(r) != winning_bsj:
            disagree_bsj[bsj_key(r)][t] = r

    struct_votes = defaultdict(list)
    for t, r in agree_map.items():
        struct_votes[struct_key(r)].append(t)

    winning_sk, struct_agree = vote_majority(struct_votes, STRUCT_PRIORITY)

    main_tools = struct_votes[winning_sk]
    struct_src = _priority_tool(set(main_tools), STRUCT_PRIORITY)
    entries.append({
        'start':          winning_bsj[0],
        'end':            winning_bsj[1],
        'struct_rec':     agree_map[struct_src],
        'bsj_src':        bsj_src,
        'struct_src':     struct_src,
        'bsj_agree':      bsj_agree,
        'struct_agree':   struct_agree,
        'isoform_label':  'main',
        'isoform_tools':  list(main_tools),
    })

    iso_n = 0

    for sk, tools in struct_votes.items():
        if sk == winning_sk:
            continue
        iso_n += 1
        src = _priority_tool(set(tools), STRUCT_PRIORITY)
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    agree_map[src],
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(tools),
        })

    for bk, bk_tools in sorted(disagree_bsj.items()):
        iso_n += 1
        bsj_src_minor    = _priority_tool(set(bk_tools), BSJ_PRIORITY)
        struct_src_minor  = _priority_tool(set(bk_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    bk_tools[struct_src_minor],
            'bsj_src':       bsj_src_minor,
            'struct_src':    struct_src_minor,
            'bsj_agree':     len(bk_tools),
            'struct_agree':  1,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(bk_tools.keys()),
        })

    return entries


def collect_entries_consensus_xstruct(tool_best, struct_tolerance):
    """
    Like collect_entries_consensus but structure voting uses ALL tools in the
    relaxed-BSJ group, comparing absolute genomic exon coordinates within
    struct_tolerance bp.

    1. BSJ vote → winning BSJ  (identical to smart_consensus)
    2. Group ALL tools by absolute exon coordinate similarity
    3. Plurality vote on those groups → winning structure
    4. Main entry: winning BSJ coords + winning structure
       - If a tool with the exact winning BSJ is in the winning struct group,
         use its record directly (no rebasing needed).
       - Otherwise, rebase the highest-priority tool from the winning group.
    5. Each non-winning struct group → one isoform entry, placed at winning
       BSJ coords (rebased as needed).  No separate minority-BSJ isoforms —
       every tool is accounted for via absolute struct grouping.
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src   = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)
    agree_map = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}

    abs_groups = group_by_abs_struct(tool_best, struct_tolerance)

    winning_idx, struct_agree = vote_struct_groups(abs_groups)
    winning_sg_tools = abs_groups[winning_idx][1]

    bsj_tools_in_winner = [t for t in winning_sg_tools if t in agree_map]
    if bsj_tools_in_winner:
        struct_src = _priority_tool(set(bsj_tools_in_winner), STRUCT_PRIORITY)
        struct_rec = agree_map[struct_src]
    else:
        # Rebase: winning structure has no representative at exact winning BSJ
        struct_src = _priority_tool(set(winning_sg_tools), STRUCT_PRIORITY)
        struct_rec = rebase_struct(tool_best[struct_src],
                                   winning_bsj[0], winning_bsj[1])

    entries.append({
        'start':         winning_bsj[0],
        'end':           winning_bsj[1],
        'struct_rec':    struct_rec,
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(winning_sg_tools),
    })

    iso_n = 0
    for i, (_, sg_tools) in enumerate(abs_groups):
        if i == winning_idx:
            continue
        iso_n += 1
        bsj_tools_in_sg = [t for t in sg_tools if t in agree_map]
        if bsj_tools_in_sg:
            src  = _priority_tool(set(bsj_tools_in_sg), STRUCT_PRIORITY)
            srec = agree_map[src]
        else:
            src  = _priority_tool(set(sg_tools), STRUCT_PRIORITY)
            srec = rebase_struct(tool_best[src],
                                 winning_bsj[0], winning_bsj[1])
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    srec,
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(sg_tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(sg_tools),
        })

    return entries


def collect_entries_consensus_hybrid(tool_best, struct_tolerance):
    """
    BSJ majority vote + absolute-coordinate structure vote restricted to
    exact-winning-BSJ tools only (no cross-BSJ borrowing, no rebasing).

    1. BSJ vote: majority vote across all tools (identical to consensus/xstruct).
    2. Structure vote: group only the exact-BSJ-winning tools by absolute genomic
       exon coordinate similarity (within struct_tolerance bp), then plurality-vote.
    3. Minority-BSJ tools emit separate isoform entries at their own BSJ coords
       (same as smart_consensus — BSJ diversity is preserved in the isoform list).

    Difference from consensus:  merges trivial exon-boundary discrepancies
                                (±struct_tolerance bp) that string equality
                                treats as separate structures.
    Difference from xstruct:    no rebasing — minority-BSJ tools never
                                contribute to the winning-BSJ structure vote.
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)

    agree_map    = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}
    disagree_bsj = defaultdict(dict)
    for t, r in tool_best.items():
        if bsj_key(r) != winning_bsj:
            disagree_bsj[bsj_key(r)][t] = r

    # Absolute-coordinate structure vote — exact-BSJ tools only
    abs_groups = group_by_abs_struct(agree_map, struct_tolerance)
    winning_idx, struct_agree = vote_struct_groups(abs_groups)
    winning_sg_tools = abs_groups[winning_idx][1]

    struct_src = _priority_tool(set(winning_sg_tools), STRUCT_PRIORITY)
    entries.append({
        'start':         winning_bsj[0],
        'end':           winning_bsj[1],
        'struct_rec':    agree_map[struct_src],
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(winning_sg_tools),
    })

    iso_n = 0
    for i, (_, sg_tools) in enumerate(abs_groups):
        if i == winning_idx:
            continue
        iso_n += 1
        src = _priority_tool(set(sg_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    agree_map[src],
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(sg_tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(sg_tools),
        })

    for bk, bk_tools in sorted(disagree_bsj.items()):
        iso_n += 1
        bsj_src_minor    = _priority_tool(set(bk_tools), BSJ_PRIORITY)
        struct_src_minor = _priority_tool(set(bk_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    bk_tools[struct_src_minor],
            'bsj_src':       bsj_src_minor,
            'struct_src':    struct_src_minor,
            'bsj_agree':     len(bk_tools),
            'struct_agree':  1,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(bk_tools.keys()),
        })

    return entries


def collect_entries_priority(tool_best):
    """
    Collect all isoform entries for smart_priority mode.

    Main entry:
      BSJ from highest BSJ-priority tool present.
      Structure from highest struct-priority tool present.

    All unique (BSJ, structure) combinations not matching the main →
    one isoform entry each.

    Returns same entry dict format as collect_entries_consensus.
    """
    entries = []

    bsj_src    = _priority_tool(set(tool_best.keys()), BSJ_PRIORITY)
    struct_src = _priority_tool(set(tool_best.keys()), STRUCT_PRIORITY)

    main_bsj = bsj_key(tool_best[bsj_src])
    main_sk  = struct_key(tool_best[struct_src])

    bsj_agree    = sum(1 for r in tool_best.values() if bsj_key(r)    == main_bsj)
    struct_agree = sum(1 for r in tool_best.values() if struct_key(r) == main_sk)

    # isoform_tools = tools that share BOTH winning BSJ AND winning struct
    main_tools = [t for t, r in tool_best.items()
                  if bsj_key(r) == main_bsj and struct_key(r) == main_sk]
    # bsj_src and struct_src must be represented even if they differ
    main_tools_set = set(main_tools) | {bsj_src, struct_src}

    entries.append({
        'start':         main_bsj[0],
        'end':           main_bsj[1],
        'struct_rec':    tool_best[struct_src],
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(main_tools_set),
    })

    # Collect remaining unique (bsj_key, struct_key) combos
    combos = defaultdict(list)
    for t, r in tool_best.items():
        combo = (bsj_key(r), struct_key(r))
        if combo != (main_bsj, main_sk):
            combos[combo].append(t)

    iso_n = 0
    for (bk, sk), tools in sorted(combos.items()):
        iso_n += 1
        src          = _priority_tool(set(tools), STRUCT_PRIORITY)
        bsj_src_iso  = _priority_tool(set(tools), BSJ_PRIORITY)
        bsj_agree_iso    = sum(1 for r in tool_best.values() if bsj_key(r) == bk)
        struct_agree_iso = sum(1 for r in tool_best.values() if struct_key(r) == sk)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    tool_best[src],
            'bsj_src':       bsj_src_iso,
            'struct_src':    src,
            'bsj_agree':     bsj_agree_iso,
            'struct_agree':  struct_agree_iso,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(tools),
        })

    return entries


# ── Output ─────────────────────────────────────────────────────────────────────

def bed12_line(chrom, start, end, name, score, strand,
               block_count, block_sizes, block_starts):
    return '\t'.join([
        chrom, str(start), str(end), name, str(score), strand,
        str(start), str(start), '0',
        block_count, block_sizes, block_starts,
    ])


def write_outputs(groups, active_tools, sample, outdir, struct_tolerance):
    os.makedirs(outdir, exist_ok=True)

    tool_flags  = active_tools
    tool_blocks = []
    for t in active_tools:
        tool_blocks += [f'{t}_block_sizes', f'{t}_block_starts']

    header = (
        ['#chrom', 'start', 'end', 'strand', 'bsj_id', 'bsj_confidence']
        + tool_flags
        + tool_blocks
        + ['isoform_confidence']
        + ['sel_block_count', 'sel_block_sizes', 'sel_block_starts']
        + ['bsj_source', 'struct_source', 'struct_agree_count', 'max_score']
        + ['isoform_label', 'isoform_tools']
    )

    mode_lines = {
        'consensus':         {'bed': [], 'tsv': ['\t'.join(header)]},
        'consensus_xstruct': {'bed': [], 'tsv': ['\t'.join(header)]},
        'consensus_hybrid':  {'bed': [], 'tsv': ['\t'.join(header)]},
        'priority':          {'bed': [], 'tsv': ['\t'.join(header)]},
    }

    for group_key, tool_map in sorted(groups.items()):
        chrom, _, _, strand = group_key

        # Best record per tool (within this group)
        tool_best = {tool: best_record(recs) for tool, recs in tool_map.items()}
        group_score = max_score_of(tool_best)

        # Per-tool flags and block columns (group-level, same for all isoforms)
        flags  = ['1' if t in tool_best else '0' for t in active_tools]
        blocks = []
        for t in active_tools:
            if t in tool_best:
                blocks += [tool_best[t]['block_sizes'], tool_best[t]['block_starts']]
            else:
                blocks += ['.', '.']

        xstruct_fn = lambda tb: collect_entries_consensus_xstruct(tb, struct_tolerance)
        hybrid_fn  = lambda tb: collect_entries_consensus_hybrid(tb, struct_tolerance)
        for mode, collect_fn in [
            ('consensus',         collect_entries_consensus),
            ('consensus_xstruct', xstruct_fn),
            ('consensus_hybrid',  hybrid_fn),
            ('priority',          collect_entries_priority),
        ]:
            entries = collect_fn(tool_best)

            for entry in entries:
                start   = entry['start']
                end     = entry['end']
                srec    = entry['struct_rec']
                label   = entry['isoform_label']
                iso_tools_str = ','.join(sorted(entry['isoform_tools']))

                bsj_id = make_bsj_id(chrom, start, end, strand,
                                      suffix=(label if label != 'main' else None))

                # Score: max among supporting tools; fall back to group max
                iso_tool_recs = {t: tool_best[t] for t in entry['isoform_tools']
                                 if t in tool_best}
                score = max_score_of(iso_tool_recs) if iso_tool_recs else group_score

                bed = bed12_line(
                    chrom, start, end, bsj_id, score, strand,
                    srec['block_count'], srec['block_sizes'], srec['block_starts']
                )

                tsv_row = '\t'.join(
                    [chrom, str(start), str(end), strand, bsj_id,
                     str(entry['bsj_agree'])]    # bsj_confidence = tool-agreement count
                    + flags
                    + blocks
                    + ['']                        # isoform_confidence placeholder
                    + [srec['block_count'], srec['block_sizes'], srec['block_starts']]
                    + [entry['bsj_src'], entry['struct_src'],
                       str(entry['struct_agree']), str(score)]  # max_score last
                    + [label, iso_tools_str]
                )

                mode_lines[mode]['bed'].append(bed)
                mode_lines[mode]['tsv'].append(tsv_row)

    outputs = {
        f'{sample}_smart_consensus.bed12':                mode_lines['consensus']['bed'],
        f'{sample}_smart_consensus_confidence.tsv':       mode_lines['consensus']['tsv'],
        f'{sample}_smart_consensus_xstruct.bed12':        mode_lines['consensus_xstruct']['bed'],
        f'{sample}_smart_consensus_xstruct_confidence.tsv': mode_lines['consensus_xstruct']['tsv'],
        f'{sample}_smart_consensus_hybrid.bed12':         mode_lines['consensus_hybrid']['bed'],
        f'{sample}_smart_consensus_hybrid_confidence.tsv': mode_lines['consensus_hybrid']['tsv'],
        f'{sample}_smart_priority.bed12':                 mode_lines['priority']['bed'],
        f'{sample}_smart_priority_confidence.tsv':        mode_lines['priority']['tsv'],
    }
    for fname, lines in outputs.items():
        with open(os.path.join(outdir, fname), 'w') as fh:
            if lines:
                fh.write('\n'.join(lines) + '\n')

    n_groups = len(groups)
    for mode in ('consensus', 'consensus_xstruct', 'consensus_hybrid', 'priority'):
        n_all = len(mode_lines[mode]['bed'])
        print(
            f'[smart_merge] {sample}: {n_groups} groups → smart_{mode}: {n_all} isoform entries',
            file=sys.stderr
        )


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if len(args.tool_names) != len(args.bed_files):
        print_error('--tool_names and --bed_files must have the same number of entries')

    all_records  = []
    active_tools = []
    for tool, path in zip(args.tool_names, args.bed_files):
        recs = read_bed12(path, tool)
        if recs:
            all_records.extend(recs)
            active_tools.append(tool)

    if len(active_tools) < 2:
        print(
            f'[smart_merge] Only {len(active_tools)} active tool(s) — need ≥2. Writing empty outputs.',
            file=sys.stderr
        )
        os.makedirs(args.outdir, exist_ok=True)
        for suffix in [
            '_smart_consensus.bed12',
            '_smart_consensus_confidence.tsv',
            '_smart_consensus_xstruct.bed12',
            '_smart_consensus_xstruct_confidence.tsv',
            '_smart_priority.bed12',
            '_smart_priority_confidence.tsv',
        ]:
            open(os.path.join(args.outdir, args.sample + suffix), 'w').close()
        return

    # Keep active_tools in priority order so TSV columns are consistent
    ordered = [t for t in (BSJ_PRIORITY + [t for t in active_tools if t not in BSJ_PRIORITY])
               if t in active_tools]

    struct_tol = args.struct_tolerance if args.struct_tolerance is not None else args.tolerance
    groups = group_relaxed(all_records, args.tolerance)
    write_outputs(groups, ordered, args.sample, args.outdir, struct_tol)


if __name__ == '__main__':
    main()
