# nf-core/nanocirc: Usage

## Introduction

**nf-core/nanocirc** is a pipeline for detection and characterisation of circular RNAs (circRNAs) from long-read nanopore sequencing data. It runs up to four detection tools in parallel, converts their outputs to a unified BED12 format, and merges results with confidence scoring when two or more tools are active.

## Samplesheet input

You must provide a samplesheet CSV file with two columns: `sample` and `fastq`.

```csv title="samplesheet.csv"
sample,fastq
SAMPLE1,/path/to/sample1.fastq.gz
SAMPLE2,/path/to/sample2.fastq.gz
```

| Column   | Description                                                                          |
| -------- | ------------------------------------------------------------------------------------ |
| `sample` | Unique sample name. Cannot contain spaces.                                           |
| `fastq`  | Full path to a gzipped FASTQ file (`.fastq.gz` or `.fq.gz`). Must be nanopore data. |

Pass the samplesheet to the pipeline with:

```bash
--input '[path to samplesheet.csv]'
```

## Running the pipeline

### Minimal example

At minimum you need the samplesheet, a reference genome FASTA, a GTF annotation, and at least one detection tool enabled:

```bash
nextflow run nf-core/nanocirc \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results/ \
    --fasta /path/to/genome.fa \
    --gtf   /path/to/annotation.gtf
```

### Full example (all four tools)

```bash
nextflow run nf-core/nanocirc \
    -profile singularity \
    --input        samplesheet.csv \
    --outdir       results/ \
    --fasta        /path/to/genome.fa \
    --gtf          /path/to/annotation.gtf \
    --circrna_db   /path/to/circ_db.bed \
    --run_isocirc  true \
    --run_circfl   true \
    --run_cirilong true \
    --run_circnick true \
    --circnick_species mouse \
    -resume
```

> [!NOTE]
> `--circrna_db` is required when `--run_isocirc` or `--run_cirilong` is set.

### Resuming a run

Add `-resume` to any command to reuse cached results from a previous run:

```bash
nextflow run nf-core/nanocirc ... -resume
```

### Params file

For repeated runs with the same settings, use a params YAML file:

```bash
nextflow run nf-core/nanocirc -profile singularity -params-file params.yaml
```

```yaml title="params.yaml"
input: samplesheet.csv
outdir: results/
fasta: /path/to/genome.fa
gtf: /path/to/annotation.gtf
circrna_db: /path/to/circ_db.bed
circnick_species: mouse
```

---

## Pipeline parameters

### Reference files

| Parameter         | Description                                                           | Required |
| ----------------- | --------------------------------------------------------------------- | -------- |
| `--fasta`         | Reference genome FASTA file                                           | Yes      |
| `--gtf`           | Gene annotation GTF file                                              | Yes      |
| `--circrna_db`    | circRNA database BED file (required for isoCirc and CIRI-long)        | Conditional |
| `--genome_index_dir` | Directory to cache BWA genome index. Defaults to `<outdir>/genome_index` | No |

### Tool selection

By default all four tools are enabled. Disable individual tools with `false`:

| Parameter       | Description               | Default |
| --------------- | ------------------------- | ------- |
| `--run_isocirc`  | Run isoCirc               | `true`  |
| `--run_circfl`   | Run CircFL-seq            | `true`  |
| `--run_cirilong` | Run CIRI-long             | `true`  |
| `--run_circnick` | Run circnick-lrs          | `true`  |

### circnick-lrs options

| Parameter                  | Description                                                                                     | Required when           |
| -------------------------- | ----------------------------------------------------------------------------------------------- | ----------------------- |
| `--circnick_species`       | Species for circnick-lrs built-in reference: `mouse` or `human`                                | `--run_circnick true`   |
| `--circnick_liftover_chain`| UCSC `.chain` file to lift circnick coordinates to the current genome build (optional)          | Never (optional)        |

> [!NOTE]
> circnick-lrs uses built-in mm10 (mouse) or hg19 (human) references internally. If your analysis uses a different genome build, provide `--circnick_liftover_chain` to convert coordinates.

### Merge options

These options control how results from multiple tools are merged and scored:

| Parameter                  | Description                                                              | Default |
| -------------------------- | ------------------------------------------------------------------------ | ------- |
| `--circrna_bsj_tolerance`  | Back-splice junction coordinate tolerance in bp for relaxed merge mode   | `5`     |
| `--circrna_isoform_overlap`| Minimum reciprocal spliced-length overlap for isoform confidence scoring | `0.95`  |

Merging is only performed when **two or more** detection tools are active.

### Consensus scoring

Each circRNA in the merged output receives a `tool_consensus` label (`Low`, `Medium`, or `High`) based on three percentage-based component scores:

- **bsj_score** — what fraction of active tools detected this back-splice junction?
- **isoform_score** — what fraction of active tools confirmed a matching isoform structure?
- **overlap_score** — what is the average pairwise spliced-length overlap across tool pairs?

Each component is binned: ≤25% of tools → 1, ≤50% → 2, ≤75% → 3, >75% → 4. The final score (sum, 3–12) maps to:

| Score | Category |
| ----- | -------- |
| 3–4   | Low      |
| 5–8   | Medium   |
| 9–12  | High     |

> [!WARNING]
> **Fewer than 4 tools reduces scoring resolution.** The pipeline emits a warning when fewer than 4 tools are active. `tool_consensus` always reflects agreement among the tools that _ran_ — a `High` from 2 tools (both agree) is not the same statistical confidence as `High` from all 4 tools. See [scoring examples](#scoring-examples) below.

### QC options

| Parameter        | Description                      | Default |
| ---------------- | -------------------------------- | ------- |
| `--skip_qc`      | Skip all QC steps                | `false` |
| `--skip_fastqc`  | Skip FastQC                      | `false` |
| `--skip_nanoplot`| Skip NanoPlot                    | `false` |
| `--skip_multiqc` | Skip MultiQC report generation   | `false` |

---

## Profiles

Use `-profile` to configure the execution environment. Multiple profiles can be combined, e.g. `-profile singularity,genouest`.

| Profile       | Description                                      |
| ------------- | ------------------------------------------------ |
| `docker`      | Run with Docker containers                       |
| `singularity` | Run with Singularity containers                  |
| `apptainer`   | Run with Apptainer containers                    |
| `test`        | Minimal test run with bundled test data          |

> [!IMPORTANT]
> This pipeline requires containers (Docker, Singularity, or Apptainer). The four circRNA detection tools are only available as container images and **conda is not supported**.

---

## Core Nextflow arguments

> [!NOTE]
> These use a single hyphen (`-`), unlike pipeline parameters which use double hyphen (`--`).

### `-resume`

Restart a pipeline reusing cached results where inputs are unchanged.

### `-work-dir` / `-w`

Directory for Nextflow working files. Defaults to `./work`. On HPC systems it is recommended to set this to a fast scratch filesystem.

### `-c`

Provide a custom Nextflow config file for tuning resource requirements or infrastructure settings.

---

## Resource requirements

Default resource labels used by the pipeline:

| Label              | CPUs | Memory  | Time   |
| ------------------ | ---- | ------- | ------ |
| `process_single`   | 1    | 6 GB    | 4 h    |
| `process_low`      | 2    | 12 GB   | 4 h    |
| `process_medium`   | 6    | 36 GB   | 8 h    |
| `process_high`     | 12   | 72 GB   | 16 h   |
| `process_long`     | 6    | 36 GB   | 120 h  |
| `process_high_memory` | 6 | 200 GB  | 16 h   |

Detection tools (isoCirc, CircFL-seq, CIRI-long, circnick-lrs) run under `process_high`. To override resources for a specific process, add to your config:

```groovy
process {
    withName: 'ISOCIRC' {
        cpus   = 16
        memory = '100.GB'
        time   = '24.h'
    }
}
```

---

## Running in the background

Use `screen`, `tmux`, or the Nextflow `-bg` flag to detach the run from your terminal:

```bash
nextflow run nf-core/nanocirc ... -bg
```

Alternatively, on HPC systems, submit the Nextflow head job itself to the scheduler:

```bash
sbatch --wrap="nextflow run nf-core/nanocirc ..."
```

---

## Troubleshooting

### Pipeline exits with "At least one tool must be active"

All four tool flags are `true` by default. This error only appears if you explicitly set all of them to `false`. Enable at least one tool.

### Pipeline exits with "CircRNA analysis requires '--circrna_db'"

isoCirc and CIRI-long require a circRNA database BED file. Either provide `--circrna_db` or disable those tools with `--run_isocirc false --run_cirilong false`.

### circnick-lrs exits with coordinate warnings

If many circRNAs report exons outside their BSJ boundaries, consider providing `--circnick_liftover_chain` to lift coordinates to the current genome build.

### Memory errors on genome indexing

If PREPARE_GENOME fails due to memory, override its resources in your config:

```groovy
process {
    withName: 'PREPARE_GENOME' {
        memory = '32.GB'
    }
}
```

---

## Scoring examples

The tables below illustrate how `tool_consensus` is assigned under different run configurations. Scores are computed as **bsj_score + isoform_score + overlap_score**, where each component is binned 1–4 based on the percentage of active tools.

### 4-tool run

| Scenario                                | BSJ tools   | Isoform tools | Avg overlap | bsj | iso | ovlp | total | Category |
| --------------------------------------- | ----------- | ------------- | ----------- | --- | --- | ---- | ----- | -------- |
| All 4 agree, full isoform match         | 4/4 (100%)  | 4/4 (100%)    | 90%         | 4   | 4   | 4    | **12** | High    |
| 3 tools agree, good isoform             | 3/4 (75%)   | 3/4 (75%)     | 80%         | 3   | 3   | 4    | **10** | High    |
| All 4 BSJ, but no isoform confirmation  | 4/4 (100%)  | 0/4 (0%)      | —           | 4   | 1   | 1    | **6**  | Medium  |
| 2 tools agree, some isoform             | 2/4 (50%)   | 2/4 (50%)     | 60%         | 2   | 2   | 3    | **7**  | Medium  |
| Only 1 tool detects this circRNA        | 1/4 (25%)   | 0/4 (0%)      | —           | 1   | 1   | 1    | **3**  | Low     |

### 3-tool run

| Scenario                                | BSJ tools   | Isoform tools | Avg overlap | bsj | iso | ovlp | total | Category |
| --------------------------------------- | ----------- | ------------- | ----------- | --- | --- | ---- | ----- | -------- |
| All 3 agree, full isoform match         | 3/3 (100%)  | 3/3 (100%)    | 90%         | 4   | 4   | 4    | **12** | High    |
| 2 tools agree, good isoform             | 2/3 (67%)   | 2/3 (67%)     | 80%         | 3   | 3   | 4    | **10** | High    |
| All 3 BSJ, but no isoform confirmation  | 3/3 (100%)  | 0/3 (0%)      | —           | 4   | 1   | 1    | **6**  | Medium  |
| 2 tools, no isoform                     | 2/3 (67%)   | 0/3 (0%)      | —           | 3   | 1   | 1    | **5**  | Medium  |
| Only 1 tool                             | 1/3 (33%)   | 0/3 (0%)      | —           | 2   | 1   | 1    | **4**  | Low     |

### 2-tool run

| Scenario                                | BSJ tools   | Isoform tools | Avg overlap | bsj | iso | ovlp | total | Category |
| --------------------------------------- | ----------- | ------------- | ----------- | --- | --- | ---- | ----- | -------- |
| Both agree, full isoform match          | 2/2 (100%)  | 2/2 (100%)    | 90%         | 4   | 4   | 4    | **12** | High    |
| Both agree, moderate isoform            | 2/2 (100%)  | 1/2 (50%)     | 55%         | 4   | 2   | 3    | **9**  | High    |
| Both agree, no isoform                  | 2/2 (100%)  | 0/2 (0%)      | —           | 4   | 1   | 1    | **6**  | Medium  |
| Only 1 tool detects                     | 1/2 (50%)   | 0/2 (0%)      | —           | 2   | 1   | 1    | **4**  | Low     |

> [!WARNING]
> `High` from 2 tools means both tools agreed with matching isoforms. With 4 tools, the same label requires independent confirmation from at least 3 tools — a substantially stronger claim. Always consider the `bsj_confidence` column (raw tool count) alongside `tool_consensus`.
