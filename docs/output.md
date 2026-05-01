# nf-core/nanocirc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory specified with `--outdir`.

## Pipeline overview

The pipeline processes long-read nanopore FASTQ files through the following steps:

1. **Quality control** — FastQC and NanoPlot assess read quality
2. **circRNA detection** — up to four tools run in parallel (isoCirc, CircFL-seq, CIRI-long, circnick-lrs)
3. **BED12 conversion** — each tool's output is converted to a unified BED12 format
4. **Merging & confidence scoring** — when two or more tools are active, results are merged using the hybrid smart-merge algorithm and scored on two independent confidence axes
5. **MultiQC** — aggregated QC report

---

## Quality control

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `qc/<sample>/fastqc/`
  - `*.html` — FastQC report with quality metrics per sample
  - `*.zip` — Archive containing the report and raw data

</details>

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports per-read quality scores, GC content, sequence length distribution, and adapter content. Useful for flagging low-quality or contaminated samples before analysis.

### NanoPlot

<details markdown="1">
<summary>Output files</summary>

- `qc/<sample>/nanoplot/`
  - `NanoPlot-report.html` — Interactive HTML report
  - `NanoStats.txt` — Summary statistics (N50, mean length, mean quality, etc.)
  - Various PNG plots: read length histogram, quality vs length scatter, etc.

</details>

[NanoPlot](https://github.com/wdecoster/NanoPlot) is designed specifically for nanopore data and provides read length distributions, quality score distributions, and yield-over-time plots.

---

## circRNA detection

Each active tool writes its raw output to a dedicated subdirectory and its BED12 file to a shared `bed12/` directory.

### isoCirc

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/isocirc/`
  - `isocirc_output/isocirc.bed` — Main isoCirc output (BED format with isoform detail)
  - `isocirc_output/` — Full isoCirc output directory

</details>

[isoCirc](https://github.com/Xinglab/isoCirc) detects circRNA isoforms from nanopore long reads using full-length read alignment.

### CircFL-seq

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/circfl_seq/`
  - `circFL_final.bed` — Final CircFL-seq output
  - Full pipeline output directory including intermediate RG, DNSC, cRG, mRG steps

</details>

[CircFL-seq](https://github.com/yangence/circfull) reconstructs full-length circRNA isoforms from nanopore sequencing using a rolling-circle amplification model.

### CIRI-long

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/ciri_long/`
  - `<sample>.info` — Main CIRI-long output with circRNA calls and isoform structure
  - `<sample>.isoforms` — Isoform-level output

</details>

[CIRI-long](https://github.com/bioinfo-biols/CIRI-long) uses a seed-and-extend strategy for circRNA detection and isoform characterisation from long reads.

### circnick-lrs

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/circnick_lrs/`
  - `<sample>/<sample>.circRNA_candidates.annotated.txt` — Annotated circRNA candidates
  - `<sample>/<sample>.circ_circRNA_exon_usage_length_of_exons.txt` — Exon usage per circRNA
  - `<sample>/<sample>.introns...intronCov.bed` — Intron coverage file
- `circrna/<sample>/circnick_lrs/lifted/` _(only if `--circnick_liftover_chain` was provided)_
  - `*_lifted_annotated.txt` — Coordinates lifted to target genome build
  - `*_lifted_exon_usage.txt`
  - `*_lifted_intron_cov.bed`
  - `*_liftover_failed.tsv` — circRNAs excluded due to failed liftover

</details>

[circnick-lrs](https://github.com/dzhang32/circnick) uses built-in mm10 or hg19 references. Provide `--circnick_liftover_chain` if your analysis uses a different genome build.

### BED12 files

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/bed12/`
  - `<sample>_cirilong.bed12` — CIRI-long output in BED12
  - `<sample>_circnick.bed12` — circnick-lrs output in BED12
  - `isocirc.bed` — isoCirc output (already in BED format)
  - `circFL_final.bed` — CircFL-seq output (already in BED format)

</details>

All tool outputs are converted to a 12-column BED format (BED12) for downstream merging. BED12 columns: chrom, start, end, name, score (read count), strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts.

---

## Merged outputs

Merging is performed when **two or more** detection tools are active. All tools are first grouped by relaxed BSJ coordinates (within `--circrna_bsj_tolerance` bp). Within each group, a merge algorithm selects the representative BSJ and exon structure. The result is confidence-scored on two independent axes and optionally filtered before publication.

### Pairwise comparisons

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/pairs/`
  - `<tool_a>_vs_<tool_b>.txt` — bedtools intersect output for each tool pair (used for isoform confidence scoring)

</details>

All pairwise combinations of active tools are compared using `bedtools intersect -split -wo` to identify shared circRNA isoforms.

### Merge algorithms

All tools within a relaxed-BSJ group (coordinates within `--circrna_bsj_tolerance` bp) are treated as candidates for the same circRNA. A merge algorithm selects the representative BSJ and exon structure from those candidates.

**BSJ selection** (all modes except `priority`): majority vote across all tools; ties broken by tool priority CIRI-long > CircFL-seq > IsoCirc > CircNick-LRS.

**Structure selection** differs between modes and is described in the tables below.

#### Default merge mode

The pipeline uses `consensus_hybrid` as its merge algorithm.

| Property | `consensus_hybrid` |
| -------- | ------------------ |
| BSJ | Majority vote across all tools |
| Structure vote participants | Tools sharing the **exact winning BSJ** only |
| Coordinate comparison | Absolute genomic coords (boundaries within `--circrna_bsj_tolerance` bp) |
| Rebasing of minority-BSJ tools | **No** — tools with a different BSJ do not contribute to the structure vote |
| Minority-BSJ tool handling | Emitted as separate isoform entries at their own BSJ coordinates |
| Structure tie-break priority | IsoCirc > CIRI-long > CircNick-LRS > CircFL-seq |

This design selects the most-supported exon structure among tools that agree on the BSJ, without shifting coordinates from tools that landed at a slightly different junction. Minority-BSJ isoforms are preserved in the output rather than silently discarded.

#### Default output files

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/smart/`
  - `<sample>_discovery.bed12` — hybrid, unfiltered (maximum recall)
  - `<sample>_discovery_confidence.tsv`
  - `<sample>_balanced.bed12` — hybrid + no_low filter (best F1)
  - `<sample>_balanced_confidence.tsv`
  - `<sample>_high_confidence.bed12` — hybrid + high_only filter (best precision)
  - `<sample>_high_confidence_confidence.tsv`

</details>

Three confidence-filtered outputs are published by default. After merging, each entry is scored on two independent confidence axes (`bsj_consensus` and `isoform_consensus`, each Low / Medium / High) and one of three filters is applied:

| Output | Filter | Rule | Axes retained |
| ------ | ------ | ---- | ------------- |
| **`discovery`** | none | Keep all entries | any |
| **`balanced`** | `no_low` | Drop entries where either axis is Low | ≥ Medium on both |
| **`high_confidence`** | `high_only` | Keep only entries where both axes are High | High on both |

#### Additional merge modes (`--run_benchmark_modes`)

Three further algorithms are available for research and benchmarking. All were outperformed by `consensus_hybrid` in benchmark evaluation and are not published by default.

| Mode | BSJ selection | Structure vote participants | Minority-BSJ rebasing | Minority-BSJ handling |
| ---- | ------------- | -------------------------- | --------------------- | --------------------- |
| `consensus` | Majority vote | Exact-BSJ tools only | No | Separate isoforms at own coords |
| `consensus_xstruct` | Majority vote | All tools in group | **Yes** — shifted to winning BSJ | Folded into structure vote |
| `priority` | Highest-priority tool | Highest-priority tool | No | Separate isoforms at own coords |

Key differences from `consensus_hybrid`:

- **`consensus`** uses string equality for structure comparison rather than coordinate similarity — two tools reporting the same exon structure with a 1 bp boundary difference are treated as distinct isoforms.
- **`consensus_xstruct`** includes minority-BSJ tools in the structure vote by rebasing their exon coordinates to the winning BSJ. This can incorporate more structural information but may introduce coordinate imprecision when BSJ offset is non-trivial.
- **`priority`** skips voting entirely — BSJ and structure come unconditionally from the single highest-priority tool present.

With `--run_benchmark_modes`, all three additional algorithms are published with four filter variants each (unfiltered, `no_low`, `trusted_only`, `high_only`).

The `trusted_only` filter is only available in benchmark mode:

| Filter | Rule |
| ------ | ---- |
| `trusted_only` | Drop Low-confidence entries unless the source tool is CIRI-long or IsoCirc |

### Confidence TSV format

All `*_confidence.tsv` files share a common format. Confidence is assessed on two **independent axes**:

- **`bsj_consensus`** — quality of BSJ detection: what fraction of active tools agreed on this back-splice junction?
- **`isoform_consensus`** — quality of exon-structure agreement: how well do tools agree on the exon boundaries of this isoform?

Each axis is scored independently (1=Low, 2–3=Medium, 4=High based on binned percentage of tools) allowing a circRNA to have a well-supported BSJ but uncertain isoform structure, or vice versa.

| Column               | Description                                                                 |
| -------------------- | --------------------------------------------------------------------------- |
| `#chrom`             | Chromosome                                                                  |
| `start`              | BSJ start (0-based)                                                         |
| `end`                | BSJ end                                                                     |
| `strand`             | Strand (`+` or `-`)                                                         |
| `bsj_id`             | Unique identifier: `chrom:start-end:strand` (isoforms suffixed `\|iso*`)   |
| `bsj_confidence`     | Number of tools detecting this BSJ (1–4)                                    |
| `<tool>`             | Per-tool presence flag: `1` if detected, `0` if not (one column per tool)  |
| `<tool>_block_sizes` | BED12 block sizes from this tool's call                                     |
| `<tool>_block_starts`| BED12 block starts from this tool's call                                    |
| `isoform_confidence` | Number of tools with confirmed isoform overlap                              |
| `bsj_score`          | Percentage of active tools detecting this BSJ, binned 1–4                  |
| `isoform_score`      | Percentage of active tools with isoform support, binned 1–4 (min 1)        |
| `overlap_score`      | Average pairwise spliced-length overlap fraction, binned 1–4               |
| `bsj_consensus`      | BSJ confidence label: `Low` (score 1), `Medium` (2–3), `High` (4)         |
| `isoform_consensus`  | Isoform confidence label: `Low` (score 1), `Medium` (2–3), `High` (4)     |

**Scoring bins (percentage of active tools):**

| % of active tools | Score | Consensus |
| ----------------- | ----- | --------- |
| ≤ 25%             | 1     | Low       |
| ≤ 50%             | 2     | Medium    |
| ≤ 75%             | 3     | Medium    |
| > 75%             | 4     | High      |

> [!NOTE]
> Consensus labels always reflect agreement among the tools that were actually run.
> A `High` from 2 tools means both tools agreed — it is not mathematically equivalent
> to `High` from 4 tools. The pipeline emits a warning when fewer than 4 tools are active.

---

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html` — Standalone HTML report viewable in any browser
  - `multiqc_data/` — Parsed statistics from all tools
  - `multiqc_plots/` — Static plot images

</details>

[MultiQC](https://multiqc.info) aggregates QC results from FastQC and NanoPlot across all samples into a single report. Skip with `--skip_multiqc`.

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html` — Nextflow execution report (resource usage per process)
  - `execution_timeline_*.html` — Timeline of all processes
  - `execution_trace_*.txt` — Raw trace file with per-task metrics
  - `pipeline_dag_*.html` — Directed acyclic graph of the pipeline
  - `nf_core_nanocirc_software_mqc_versions.yml` — Software versions for all tools

</details>

Nextflow automatically generates execution reports for every run. These are useful for troubleshooting, optimising resource requests, and recording the exact software versions used.
