# bcl2fastq_demux
 
WDL workflow to demultiplex Illumina NextSeq BCL data into paired-end gzipped FASTQ files.
 
## Overview
 
Takes a **tar.gz archive** of an Illumina run directory and an Illumina-format **SampleSheet.csv**, and produces per-sample `_R1.fastq.gz` / `_R2.fastq.gz` file pairs.
 
Uses [`bcl_to_fastq`](https://github.com/brwnj/bcl2fastq) (bcl2fastq-nextseq v1.3.0), which wraps Illumina's `bcl2fastq2` and automatically concatenates reads across lanes.
 
**Container:** `quay.io/biocontainers/bcl2fastq-nextseq:1.3.0--pyh5e36f6f_0`
 
## Inputs
 
| Parameter | Type | Default | Description |
|---|---|---|---|
| `bcl_tar_gz` | File | *required* | tar.gz of the Illumina run directory |
| `sample_sheet` | File | *required* | Illumina SampleSheet.csv |
| `reverse_complement` | Boolean | `true` | Reverse-complement index 2 (required for most NextSeq dual-index libraries) |
| `barcode_mismatches` | Int | `0` | Allowed mismatches per index |
| `loading_threads` | Int | `4` | Threads for BCL loading |
| `demux_threads` | Int | `4` | Threads for demultiplexing |
| `processing_threads` | Int | `8` | Threads for FASTQ processing |
| `writing_threads` | Int | `4` | Threads for FASTQ writing |
| `memory_gb` | Int | `16` | RAM in GiB |
| `disk_gb` | Int | `500` | Disk in GiB (must exceed uncompressed BCL size) |
| `preemptible` | Int | `1` | Preemptible VM retries (Google Cloud) |
 
## Outputs
 
| Output | Type | Description |
|---|---|---|
| `fastq_r1` | Array[File] | R1 FASTQ files, one per sample |
| `fastq_r2` | Array[File] | R2 FASTQ files, one per sample |
| `demux_stats` | File | Per-sample read counts (demultiplexing_stats.csv) |
| `bcl2fastq_log` | File | Raw bcl2fastq conversion log |
 
## Usage
 
### Terra.bio
 
1. Import the workflow from Dockstore into your Terra workspace.
2. Upload your run archive and sample sheet to a Google Cloud Storage bucket.
3. Fill in the inputs JSON (see `bcl2fastq_demux_inputs.json` for a template).
4. Launch the workflow.
 
### Cromwell (local)
 
```bash
java -jar cromwell.jar run bcl2fastq_demux.wdl --inputs bcl2fastq_demux_inputs.json
```
 
## Notes
 
- The tool uses `--no-wait` so it does not poll for RTAComplete; it assumes your tar.gz is a completed run.
- Undetermined FASTQ files are excluded from outputs by design.
- For single-index libraries, set `reverse_complement` to `false`.
- `disk_gb` should be at least 3× the size of the compressed tar.gz to accommodate extraction and FASTQ output.
