# nf-core/nanocirc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory specified with `--outdir`.

## Pipeline overview

The pipeline processes long-read nanopore FASTQ files through the following steps:

1. **Quality control** — FastQC and NanoPlot assess read quality
2. **circRNA detection** — up to four tools run in parallel (isoCirc, CircFL-seq, CIRI-long, circnick-lrs)
3. **BED12 conversion** — each tool's output is converted to a unified BED12 format
4. **Merging & confidence scoring** — when two or more tools are active, results are merged and scored
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

Merging is performed when **two or more** detection tools are active. Two merge strategies are applied: BSJ-based (strict and relaxed) and exon-based.

### Pairwise comparisons

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/pairs/`
  - `<tool_a>_vs_<tool_b>.txt` — bedtools intersect output for each tool pair (used for isoform confidence scoring)

</details>

All pairwise combinations of active tools are compared using `bedtools intersect -split -wo` to identify shared circRNA isoforms.

### BSJ-based merge (strict and relaxed)

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/strict/`
  - `<sample>_strict_union.bed12` — All circRNAs detected by any tool (exact BSJ match)
  - `<sample>_strict_union_confidence.tsv` — Confidence table for union
  - `<sample>_strict_intersection.bed12` — circRNAs detected by all active tools
  - `<sample>_strict_intersection_confidence.tsv` — Confidence table for intersection

- `circrna/<sample>/merged/relaxed/`
  - `<sample>_relaxed_union.bed12` — Union with ±`circrna_bsj_tolerance` bp BSJ tolerance
  - `<sample>_relaxed_union_confidence.tsv`
  - `<sample>_relaxed_intersection.bed12`
  - `<sample>_relaxed_intersection_confidence.tsv`

</details>

**Strict mode** groups circRNAs by exact back-splice junction (chrom, start, end, strand).
**Relaxed mode** groups circRNAs whose BSJ coordinates differ by at most `--circrna_bsj_tolerance` bp (default: 5 bp) on each end.

### Exon-based merge

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/exon_based/`
  - `<sample>_exon_union.bed12`
  - `<sample>_exon_union_confidence.tsv`
  - `<sample>_exon_intersection.bed12`
  - `<sample>_exon_intersection_confidence.tsv`

</details>

Groups circRNAs by exon structure similarity (spliced-length overlap ≥ `--circrna_isoform_overlap`, default: 0.95) using a union-find algorithm on pairwise bedtools overlap data. circRNAs with the same BSJ but different isoforms may be grouped differently than in the BSJ-based merge.

### Confidence TSV format

All `*_confidence.tsv` files share a common format:

| Column              | Description                                                                 |
| ------------------- | --------------------------------------------------------------------------- |
| `#chrom`            | Chromosome                                                                  |
| `start`             | BSJ start (0-based)                                                         |
| `end`               | BSJ end                                                                     |
| `strand`            | Strand (`+` or `-`)                                                         |
| `bsj_id`            | Unique BSJ identifier: `chrom:start-end:strand`                             |
| `bsj_confidence`    | Number of tools detecting this BSJ (1–4)                                    |
| `<tool>`            | Per-tool presence flag: `1` if detected, `0` if not (one column per tool)  |
| `<tool>_block_sizes`| BED12 block sizes from this tool's call                                     |
| `<tool>_block_starts`| BED12 block starts from this tool's call                                   |
| `isoform_confidence`| Number of tools with confirmed isoform overlap                              |
| `bsj_score`         | Percentage of active tools detecting this BSJ, binned 1–4                  |
| `isoform_score`     | Percentage of active tools with isoform support, binned 1–4, minimum 1     |
| `overlap_score`     | Average pairwise spliced-length overlap fraction, binned 1–4               |
| `final_score`       | Sum of three scores (3–12)                                                  |
| `tool_consensus`    | Consensus category: `Low` (3–4), `Medium` (5–8), `High` (9–12)            |

**Percentage-based scoring bins (for bsj_score and isoform_score):**

| % of active tools | Score |
| ----------------- | ----- |
| ≤ 25%             | 1     |
| ≤ 50%             | 2     |
| ≤ 75%             | 3     |
| > 75%             | 4     |

**Final score categories:**

| Score range | Category |
| ----------- | -------- |
| 3–4         | Low      |
| 5–8         | Medium   |
| 9–12        | High     |

> [!NOTE]
> `tool_consensus` always reflects agreement among the tools that were actually run.
> A `High` from 2 tools means both tools agreed — it is not mathematically equivalent
> to `High` from 4 tools. The pipeline emits a warning when fewer than 4 tools are active.

A circRNA detected by all active tools with confirmed isoform overlap will receive a `High` consensus score. A circRNA detected by only one tool, or a minority of tools, will receive `Low`.

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
