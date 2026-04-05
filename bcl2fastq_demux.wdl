version 1.0

## bcl2fastq_demux.wdl
##
## Demultiplexes an Illumina BCL run directory (supplied as a tar.gz archive)
## and a CSV sample sheet into paired-end gzipped FASTQ files using
## Illumina BCL Convert v4.4.6.
##
## Binary: bcl-convert (NOT bcl2fastq — this is the newer Illumina tool)
## Docker image: idolawoye/bcl-convert:4.4.6.2 (Docker Hub, public)
##
## Output structure
## ----------------
## The workflow emits three parallel arrays (sample_ids, fastq_r1, fastq_r2)
## that are index-matched — element [i] of each array belongs to the same
## sample. Terra can load these directly into a per-sample data table using
## the "Import from workflow outputs" feature, giving one row per sample.
##
## Sample sheet format
## --------------------
## BCL Convert 4.x expects a v2 sample sheet with sections:
##   [Header], [Reads], [BCLConvert_Settings], [BCLConvert_Data]
## The legacy bcl2fastq v1 format ([Data] with Sample_ID/index columns) is
## also accepted. Index RC handling is read from RunInfo.xml automatically.
##
## Disk sizing strategy
## --------------------
## Disk budget must cover three things simultaneously:
## 1. The compressed tar.gz itself (localized by Cromwell before the script runs)
## 2. The extracted BCL directory (~5-10x compressed size)
## 3. The FASTQ output (~2-3x compressed size)
## Disk is sized dynamically as:
##   ceil(size(bcl_tar_gz, "GiB")) * disk_multiplier + disk_overhead_gb
## Defaults: multiplier=20, overhead=100 GiB.
## Raise disk_multiplier if localization itself fails with "No space left on device".
##   A good formula: set disk_multiplier = ceil(total_needed_GiB / compressed_GiB)
##
## Inputs
## ------
##   bcl_tar_gz        : tar.gz of the Illumina run directory
##   sample_sheet      : Illumina SampleSheet.csv (v1 or v2 format)
##   no_lane_splitting : merge all lanes into one FASTQ per sample [default true]
##   memory_gb         : RAM in GiB; 32 GiB recommended for large runs [default 32]
##   disk_multiplier   : multiplier on compressed input size [default 20]
##   disk_overhead_gb  : flat GiB overhead [default 100]
##   preemptible       : preemptible VM retries [default 1]
##
## Outputs
## -------
##   sample_ids      : Array[String] — sample IDs parsed from FASTQ filenames
##   fastq_r1        : Array[File]   — R1 FASTQ per sample (index-matched to sample_ids)
##   fastq_r2        : Array[File]   — R2 FASTQ per sample (index-matched to sample_ids)
##   bcl_convert_log : bcl-convert stdout+stderr log

workflow bcl2fastq_demux {

    meta {
        author: "Idowu Olawoye"
        email: "idowuolawoye@gmail.com"
        description: "Demultiplex an Illumina BCL run (tar.gz) and SampleSheet.csv into paired-end gzipped FASTQ files using BCL Convert v4.4.6. Outputs parallel sample_ids/fastq_r1/fastq_r2 arrays for direct Terra data table import."
    }

    parameter_meta {
        bcl_tar_gz: {
            description: "tar.gz archive of the Illumina run directory",
            localization_optional: false
        }
        sample_sheet: {
            description: "Illumina SampleSheet.csv (v1 or v2 format)",
            localization_optional: false
        }
        no_lane_splitting: {
            description: "Merge all lanes into a single FASTQ per sample. When true, output filenames omit the L00N lane component.",
            default: true
        }
        memory_gb: {
            description: "RAM in GiB. bcl-convert requires at least 2 GiB free after OS overhead; 32 GiB recommended for large runs.",
            default: 32
        }
        disk_multiplier: {
            description: "Multiplier on compressed BCL input size for disk estimate. Must cover: the tar.gz itself (1x, localized before script runs) + BCL extraction (~5-10x) + FASTQ output (~2-3x). Default 20 is conservative; raise further if localization fails with No space left on device.",
            default: 20
        }
        disk_overhead_gb: {
            description: "Flat GiB added on top of the scaled disk estimate for OS, container layers, and tmp files.",
            default: 100
        }
        preemptible: {
            description: "Number of preemptible VM retries (Google Cloud).",
            default: 1
        }
    }

    input {
        File    bcl_tar_gz
        File    sample_sheet
        Boolean no_lane_splitting = true
        Int     memory_gb         = 32
        Int     disk_multiplier   = 20
        Int     disk_overhead_gb  = 100
        Int     preemptible       = 1
    }

    # Dynamic disk sizing based on actual compressed input size
    Int disk_gb = ceil(size(bcl_tar_gz, "GiB")) * disk_multiplier + disk_overhead_gb

    call BclConvert {
        input:
            bcl_tar_gz        = bcl_tar_gz,
            sample_sheet      = sample_sheet,
            no_lane_splitting = no_lane_splitting,
            memory_gb         = memory_gb,
            disk_gb           = disk_gb,
            preemptible       = preemptible
    }

    # Parse sample IDs from FASTQ filenames for Terra data table import.
    # Inputs are typed Array[String] — Cromwell does NOT localize (download)
    # String inputs, so the 378 FASTQ files are never re-downloaded here.
    # Only the GCS path strings are passed; the Python script works on basenames.
    call ParseSampleIds {
        input:
            fastq_r1_paths = BclConvert.fastq_r1,
            fastq_r2_paths = BclConvert.fastq_r2,
            preemptible    = preemptible
    }

    output {
        Array[String] sample_ids      = ParseSampleIds.sample_ids
        Array[String] fastq_r1        = ParseSampleIds.fastq_r1_sorted
        Array[String] fastq_r2        = ParseSampleIds.fastq_r2_sorted
        File          bcl_convert_log = BclConvert.bcl_convert_log
    }
}

# ── Task 1: BCL conversion ────────────────────────────────────────────────────

task BclConvert {

    meta {
        description: "Extracts the BCL tar.gz and runs bcl-convert 4.4.6 to produce paired-end FASTQ files."
    }

    input {
        File    bcl_tar_gz
        File    sample_sheet
        Boolean no_lane_splitting
        Int     memory_gb
        Int     disk_gb
        Int     preemptible
    }

    String lane_splitting_arg = if no_lane_splitting then "true" else "false"

    command <<<
        set -euo pipefail

        WORKDIR="$(pwd)"

        echo "[INFO] bcl-convert version:"
        bcl-convert --version

        # Report provisioned disk size so future failures are easier to diagnose
        echo "[INFO] Provisioned disk (~{disk_gb} GiB requested):"
        df -h .

        mkdir -p "${WORKDIR}/run_dir"
        echo "[INFO] Extracting: ~{bcl_tar_gz}"
        tar -xzf "~{bcl_tar_gz}" -C "${WORKDIR}/run_dir" --strip-components=1
        echo "[INFO] Disk after extraction:"
        df -h .

        cp "~{sample_sheet}" "${WORKDIR}/run_dir/SampleSheet.csv"

        mkdir -p "${WORKDIR}/fastqs_out"
        mkdir -p "${WORKDIR}/out_r1"
        mkdir -p "${WORKDIR}/out_r2"

        echo "[INFO] Starting bcl-convert"
        BCL_EXIT=0
        bcl-convert \
            --bcl-input-directory    "${WORKDIR}/run_dir" \
            --output-directory       "${WORKDIR}/fastqs_out" \
            --sample-sheet           "${WORKDIR}/run_dir/SampleSheet.csv" \
            --no-lane-splitting      ~{lane_splitting_arg} \
            --bcl-only-matched-reads true \
            --force \
            2>&1 | tee "${WORKDIR}/bcl-convert.log" || BCL_EXIT=${PIPESTATUS[0]}

        if [ "${BCL_EXIT}" -ne 0 ]; then
            echo "[ERROR] bcl-convert exited with code ${BCL_EXIT}."
            exit "${BCL_EXIT}"
        fi

        echo "[INFO] Disk after bcl-convert:"
        df -h .

        # Move FASTQs — no findutils in this container, use bash globbing
        echo "[INFO] Moving R1 FASTQs"
        for f in "${WORKDIR}"/fastqs_out/*_R1_001.fastq.gz; do
            case "${f##*/}" in Undetermined*) continue ;; esac
            mv "${f}" "${WORKDIR}/out_r1/"
        done

        echo "[INFO] Moving R2 FASTQs"
        for f in "${WORKDIR}"/fastqs_out/*_R2_001.fastq.gz; do
            case "${f##*/}" in Undetermined*) continue ;; esac
            mv "${f}" "${WORKDIR}/out_r2/"
        done

        R1_COUNT=$(ls "${WORKDIR}/out_r1/"*_R1_001.fastq.gz 2>/dev/null | wc -l)
        R2_COUNT=$(ls "${WORKDIR}/out_r2/"*_R2_001.fastq.gz 2>/dev/null | wc -l)
        echo "[INFO] R1 files: ${R1_COUNT}  R2 files: ${R2_COUNT}"

        if [ "${R1_COUNT}" -eq 0 ] || [ "${R2_COUNT}" -eq 0 ]; then
            echo "[ERROR] No FASTQ files found — bcl-convert may have been killed (OOM)."
            ls -lh "${WORKDIR}/fastqs_out/" || true
            exit 1
        fi

        ls -lh "${WORKDIR}/out_r1/"
        ls -lh "${WORKDIR}/out_r2/"
        echo "[INFO] Final disk:"; df -h .
    >>>

    output {
        Array[File] fastq_r1        = glob("out_r1/*_R1_001.fastq.gz")
        Array[File] fastq_r2        = glob("out_r2/*_R2_001.fastq.gz")
        File        bcl_convert_log = "bcl-convert.log"
    }

    runtime {
        docker:      "idolawoye/bcl-convert:4.4.6.2"
        memory:      "~{memory_gb} GiB"
        disks:       "local-disk ~{disk_gb} SSD"
        cpu:         8
        preemptible: preemptible
        maxRetries:  1
    }
}

# ── Task 2: Parse sample IDs from FASTQ filenames ────────────────────────────
#
# bcl-convert names FASTQs as:
#   <Sample_ID>_S<N>_R[12]_001.fastq.gz          (no-lane-splitting)
#   <Sample_ID>_S<N>_L00<N>_R[12]_001.fastq.gz   (per-lane)
#
# IMPORTANT: inputs are Array[String], NOT Array[File].
# Cromwell localizes (re-downloads) File inputs onto the task VM.
# With 189 samples that means downloading ~378 large FASTQ files onto a
# small disk just to read their filenames — which caused the disk-full error.
# Passing paths as String prevents any localization. The Python script only
# needs os.path.basename() on the path string; the files stay in GCS.
#
# This task:
#   1. Sorts R1 and R2 path arrays by filename so they are index-matched.
#   2. Extracts the Sample_ID prefix (everything before _S<digits>_).
#   3. Writes the three parallel arrays as output files that WDL reads back.
#
# Uses python:3.11-slim — tiny image, no downloads needed.

task ParseSampleIds {

    meta {
        description: "Parses sample IDs from FASTQ filenames and emits index-matched sample_ids, fastq_r1, fastq_r2 arrays for Terra data table import."
    }

    input {
        Array[String] fastq_r1_paths
        Array[String] fastq_r2_paths
        Int           preemptible
    }

    command <<<
        set -euo pipefail

        python3 - << 'PYEOF'
import os
import re

# WDL write_json() writes the array to a temp file; read it back as lines.
# In WDL 1.0 command blocks, ~{sep="\n" fastq_r1} interpolates the array
# as a newline-separated string — use that directly instead of JSON.
r1_sorted = sorted(
    [p.strip() for p in """~{sep="\n" fastq_r1_paths}""".splitlines() if p.strip()],
    key=os.path.basename
)
r2_sorted = sorted(
    [p.strip() for p in """~{sep="\n" fastq_r2_paths}""".splitlines() if p.strip()],
    key=os.path.basename
)

if len(r1_sorted) != len(r2_sorted):
    raise ValueError(
        f"R1 count ({len(r1_sorted)}) != R2 count ({len(r2_sorted)}). "
        "Arrays cannot be index-matched."
    )

def extract_sample_id(path):
    name = os.path.basename(path)
    m = re.match(r'^(.+?)_S\d+_', name)
    if not m:
        raise ValueError(f"Cannot parse sample ID from filename: {name}")
    return m.group(1)

sample_ids = [extract_sample_id(p) for p in r1_sorted]
r2_ids     = [extract_sample_id(p) for p in r2_sorted]

mismatches = [
    (i, sample_ids[i], r2_ids[i])
    for i in range(len(sample_ids))
    if sample_ids[i] != r2_ids[i]
]
if mismatches:
    raise ValueError(f"R1/R2 sample ID mismatches: {mismatches}")

with open("sample_ids.txt", "w") as fh:
    fh.write("\n".join(sample_ids) + "\n")

with open("r1_sorted.txt", "w") as fh:
    fh.write("\n".join(r1_sorted) + "\n")

with open("r2_sorted.txt", "w") as fh:
    fh.write("\n".join(r2_sorted) + "\n")

print(f"[INFO] Parsed {len(sample_ids)} samples.")
for sid, r1, r2 in zip(sample_ids, r1_sorted, r2_sorted):
    print(f"  {sid}")
    print(f"    R1: {os.path.basename(r1)}")
    print(f"    R2: {os.path.basename(r2)}")
PYEOF
    >>>

    output {
        Array[String] sample_ids      = read_lines("sample_ids.txt")
        Array[String] fastq_r1_sorted = read_lines("r1_sorted.txt")
        Array[String] fastq_r2_sorted = read_lines("r2_sorted.txt")
    }

    runtime {
        docker:      "python:3.11-slim"
        memory:      "2 GiB"
        disks:       "local-disk 10 HDD"
        cpu:         1
        preemptible: preemptible
        maxRetries:  1
    }
}
