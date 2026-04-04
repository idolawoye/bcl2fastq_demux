version 1.0

## bcl2fastq_demux.wdl
##
## Demultiplexes an Illumina BCL run directory (supplied as a tar.gz archive)
## and a CSV sample sheet into paired-end gzipped FASTQ files using
## Illumina BCL Convert v4.4.6.
##
## Binary: bcl-convert (NOT bcl2fastq — this is the newer Illumina tool)
##
##
## Sample sheet format
## --------------------
## BCL Convert 4.x expects a v2 sample sheet with sections:
##   [Header], [Reads], [BCLConvert_Settings], [BCLConvert_Data]
## The legacy bcl2fastq v1 format ([Data] with Sample_ID/index columns) is
## also accepted. Index RC handling is read from RunInfo.xml automatically;
## you do not need to set it manually for NextSeq/NovaSeq runs.
##
## Disk sizing strategy
## --------------------
## BCL archives expand ~5-10x when extracted; FASTQ output adds ~2-3x more.
## Disk is sized dynamically as:
##   ceil(size(bcl_tar_gz, "GiB")) * disk_multiplier + disk_overhead_gb
## Defaults: multiplier=15, overhead=50 GiB.
##
## Inputs
## ------
##   bcl_tar_gz         : tar.gz of the Illumina run directory
##   sample_sheet       : Illumina SampleSheet.csv (v1 or v2 format)
##   no_lane_splitting  : merge all lanes into one FASTQ per sample [default true]
##   memory_gb          : RAM in GiB [default 16]
##   disk_multiplier    : multiplier on compressed input size [default 15]
##   disk_overhead_gb   : flat GiB overhead [default 50]
##   preemptible        : preemptible VM retries [default 1]
##
## Outputs
## -------
##   fastq_r1          : Array of R1 FASTQ files (one per sample)
##   fastq_r2          : Array of R2 FASTQ files (one per sample)
##   bcl_convert_log   : bcl-convert stdout+stderr log

workflow bcl2fastq_demux {

    meta {
        author: "Idowu Olawoye"
        email:	"idowuolawoye@gmail.com"
        description: "Demultiplex an Illumina BCL run (tar.gz) and SampleSheet.csv into paired-end gzipped FASTQ files using BCL Convert v4.4.6."
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
            description: "RAM to allocate in GiB.",
            default: 16
        }
        disk_multiplier: {
            description: "Multiplier on compressed BCL input size for disk estimate.",
            default: 15
        }
        disk_overhead_gb: {
            description: "Flat GiB added on top of the scaled disk estimate.",
            default: 50
        }
        preemptible: {
            description: "Number of preemptible VM retries (Google Cloud).",
            default: 1
        }
    }

    input {
        File    bcl_tar_gz
        File    sample_sheet
        Boolean no_lane_splitting  = true
        Int     memory_gb          = 16
        Int     disk_multiplier    = 15
        Int     disk_overhead_gb   = 50
        Int     preemptible        = 1
    }

    # Dynamic disk sizing based on actual compressed input size
    Int disk_gb = ceil(size(bcl_tar_gz, "GiB")) * disk_multiplier + disk_overhead_gb

    call BclConvert {
        input:
            bcl_tar_gz         = bcl_tar_gz,
            sample_sheet       = sample_sheet,
            no_lane_splitting  = no_lane_splitting,
            memory_gb          = memory_gb,
            disk_gb            = disk_gb,
            preemptible        = preemptible
    }

    output {
        Array[File] fastq_r1        = BclConvert.fastq_r1
        Array[File] fastq_r2        = BclConvert.fastq_r2
        File        bcl_convert_log = BclConvert.bcl_convert_log
    }
}

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

    # bcl-convert --no-lane-splitting takes an explicit boolean argument: true/false
    String lane_splitting_arg = if no_lane_splitting then "true" else "false"

    command <<<
        set -euo pipefail

        WORKDIR="$(pwd)"

        # ── 0. Confirm binary ─────────────────────────────────────────────────
        echo "[INFO] bcl-convert version:"
        bcl-convert --version

        # ── 1. Disk before extraction ─────────────────────────────────────────
        echo "[INFO] Disk before extraction:"
        df -h .

        # ── 2. Extract BCL run archive ────────────────────────────────────────
        mkdir -p "${WORKDIR}/run_dir"
        echo "[INFO] Extracting: ~{bcl_tar_gz}"
        tar -xzf "~{bcl_tar_gz}" -C "${WORKDIR}/run_dir" --strip-components=1
        echo "[INFO] Disk after extraction:"
        df -h .

        # ── 3. Stage the sample sheet at the run directory root ───────────────
        # bcl-convert searches --bcl-input-directory for SampleSheet.csv by
        # default, so placing it there means we don't need --sample-sheet.
        # We pass --sample-sheet explicitly anyway for clarity.
        cp "~{sample_sheet}" "${WORKDIR}/run_dir/SampleSheet.csv"

        # ── 4. Create output directory ────────────────────────────────────────
        mkdir -p "${WORKDIR}/fastqs_out"

        # ── 5. Run bcl-convert ────────────────────────────────────────────────
        echo "[INFO] Starting bcl-convert"
        bcl-convert \
            --bcl-input-directory   "${WORKDIR}/run_dir" \
            --output-directory      "${WORKDIR}/fastqs_out" \
            --sample-sheet          "${WORKDIR}/run_dir/SampleSheet.csv" \
            --no-lane-splitting     ~{lane_splitting_arg} \
            --bcl-only-matched-reads true \
            --force \
            2>&1 | tee "${WORKDIR}/bcl-convert.log"

        echo "[INFO] Disk after bcl-convert:"
        df -h .

        # ── 6. Collect FASTQ outputs ──────────────────────────────────────────
        # With --no-lane-splitting true:  <Sample_ID>_S#_R[12]_001.fastq.gz
        # With --no-lane-splitting false: <Sample_ID>_S#_L00#_R[12]_001.fastq.gz
        # In both cases the glob *_R1_001.fastq.gz / *_R2_001.fastq.gz matches.
        # FASTQs land directly in fastqs_out/
        mkdir -p "${WORKDIR}/out_r1" "${WORKDIR}/out_r2"

        find "${WORKDIR}/fastqs_out" -maxdepth 1 \
            -name "*_R1_001.fastq.gz" ! -name "Undetermined*" \
            -exec mv {} "${WORKDIR}/out_r1/" \;

        find "${WORKDIR}/fastqs_out" -maxdepth 1 \
            -name "*_R2_001.fastq.gz" ! -name "Undetermined*" \
            -exec mv {} "${WORKDIR}/out_r2/" \;

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
