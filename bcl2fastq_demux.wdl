version 1.0

## bcl2fastq_demux.wdl
##
## Demultiplexes an Illumina NextSeq BCL run directory (supplied as a tar.gz
## archive) and a CSV sample sheet into paired-end gzipped FASTQ files.
##
## Tool: bcl_to_fastq (brwnj/bcl2fastq) v1.3.0
## Container: quay.io/biocontainers/bcl2fastq-nextseq:1.3.0--pyh5e36f6f_0
##
## Disk sizing strategy
## --------------------
## BCL archives expand to ~5-10x their compressed size on disk, and the FASTQ
## output adds roughly another 2-3x of the compressed input. The disk is
## therefore sized dynamically as:
##
##   ceil(bcl_tar_gz_bytes / 1 GiB) * disk_multiplier + disk_overhead_gb
##
## with safe defaults of multiplier=15 and overhead=50 GiB. This usually
## avoids "No space left on device" errors without over-provisioning.
## Override disk_multiplier or disk_overhead_gb if your run is unusually
## large or your FASTQ output is much bigger than the BCL input.
##
## Inputs
## ------
##   bcl_tar_gz        : tar.gz archive of the Illumina run directory
##                       (must contain RunInfo.xml, RTAComplete.txt, etc.)
##   sample_sheet      : Illumina-format SampleSheet.csv
##   reverse_complement: reverse-complement index 2 (required for most
##                       NextSeq dual-index libraries) [default true]
##   barcode_mismatches: allowed mismatches per index [default 0]
##   loading_threads   : BCL loading threads  [default 4]
##   demux_threads     : demultiplexing threads [default 4]
##   processing_threads: FASTQ processing threads [default 8]
##   writing_threads   : FASTQ writing threads  [default 4]
##   memory_gb         : RAM to allocate (GiB) [default 16]
##   disk_multiplier   : multiplier on compressed input size for disk [default 15]
##   disk_overhead_gb  : flat GiB added on top of the scaled disk [default 50]
##   preemptible       : number of preemptible retries [default 1]
##
## Outputs
## -------
##   fastq_r1     : Array of R1 FASTQ files (one per sample)
##   fastq_r2     : Array of R2 FASTQ files (one per sample)
##   demux_stats  : demultiplexing_stats.csv produced by bcl_to_fastq
##   bcl2fastq_log: raw bcl2fastq log

workflow bcl2fastq_demux {

    meta {
        author: "Idowu Olawoye"
        email: "idowuolawoye@gmail.com"
        description: "Demultiplex a NextSeq BCL run (tar.gz) and SampleSheet.csv into paired-end gzipped FASTQ files using bcl_to_fastq (bcl2fastq-nextseq)."
    }

    parameter_meta {
        bcl_tar_gz: {
            description: "tar.gz archive of the Illumina run directory",
            localization_optional: false
        }
        sample_sheet: {
            description: "Illumina SampleSheet.csv",
            localization_optional: false
        }
        reverse_complement: {
            description: "Reverse-complement index 2 of the sample sheet. Required for most NextSeq dual-index runs.",
            default: true
        }
        barcode_mismatches: {
            description: "Number of allowed mismatches per index barcode.",
            default: 0
        }
        loading_threads: {
            description: "Number of threads for loading BCL data.",
            default: 4
        }
        demux_threads: {
            description: "Number of threads for demultiplexing.",
            default: 4
        }
        processing_threads: {
            description: "Number of threads for processing demultiplexed data.",
            default: 8
        }
        writing_threads: {
            description: "Number of threads for writing FASTQ data.",
            default: 4
        }
        memory_gb: {
            description: "RAM to allocate in GiB.",
            default: 16
        }
        disk_multiplier: {
            description: "Multiplier applied to the compressed BCL input size to estimate total disk needed. Default 15 covers extraction (~5-10x) plus FASTQ output (~2-3x).",
            default: 15
        }
        disk_overhead_gb: {
            description: "Flat GiB added on top of the scaled disk estimate for OS, container layers, and log files.",
            default: 50
        }
        preemptible: {
            description: "Number of times to allow preemptible VM retries (Google Cloud).",
            default: 1
        }
    }

    input {
        File    bcl_tar_gz
        File    sample_sheet
        Boolean reverse_complement    = true
        Int     barcode_mismatches    = 0
        Int     loading_threads       = 4
        Int     demux_threads         = 4
        Int     processing_threads    = 8
        Int     writing_threads       = 4
        Int     memory_gb             = 16
        Int     disk_multiplier       = 15
        Int     disk_overhead_gb      = 50
        Int     preemptible           = 1
    }

    # Compute disk size from actual input file size at workflow-parse time.
    # size() returns GiB; ceil() rounds up to the nearest whole GiB.
    Int disk_gb = ceil(size(bcl_tar_gz, "GiB")) * disk_multiplier + disk_overhead_gb

    call Demultiplex {
        input:
            bcl_tar_gz         = bcl_tar_gz,
            sample_sheet       = sample_sheet,
            reverse_complement = reverse_complement,
            barcode_mismatches = barcode_mismatches,
            loading_threads    = loading_threads,
            demux_threads      = demux_threads,
            processing_threads = processing_threads,
            writing_threads    = writing_threads,
            memory_gb          = memory_gb,
            disk_gb            = disk_gb,
            preemptible        = preemptible
    }

    output {
        Array[File] fastq_r1    = Demultiplex.fastq_r1
        Array[File] fastq_r2    = Demultiplex.fastq_r2
        File        demux_stats = Demultiplex.demux_stats
        File        bcl2fastq_log = Demultiplex.bcl2fastq_log
    }
}

task Demultiplex {

    meta {
        description: "Extracts the BCL tar.gz, places the sample sheet, runs bcl_to_fastq, and collects paired-end FASTQ outputs."
    }

    input {
        File    bcl_tar_gz
        File    sample_sheet
        Boolean reverse_complement
        Int     barcode_mismatches
        Int     loading_threads
        Int     demux_threads
        Int     processing_threads
        Int     writing_threads
        Int     memory_gb
        Int     disk_gb
        Int     preemptible
    }

    # Build the optional --reverse-complement flag string once
    String rc_flag = if reverse_complement then "--reverse-complement" else ""

    command <<<
        set -euo pipefail

        # ── 0. Report available disk space for debugging ─────────────────────
        echo "[INFO] Disk available before extraction:"
        df -h .

        # ── 1. Extract the BCL run archive ──────────────────────────────────
        mkdir -p run_dir
        echo "[INFO] Extracting BCL archive: ~{bcl_tar_gz}"
        tar -xzf "~{bcl_tar_gz}" -C run_dir --strip-components=1
        echo "[INFO] Disk after extraction:"
        df -h .

        # ── 2. Place the sample sheet inside the run directory ──────────────
        # bcl_to_fastq expects SampleSheet.csv at the root of the run folder
        cp "~{sample_sheet}" run_dir/SampleSheet.csv

        # ── 3. Run bcl_to_fastq ─────────────────────────────────────────────
        echo "[INFO] Starting bcl_to_fastq"
        bcl_to_fastq \
            --runfolder run_dir \
            --loading        ~{loading_threads} \
            --demultiplexing ~{demux_threads} \
            --processing     ~{processing_threads} \
            --writing        ~{writing_threads} \
            --barcode-mismatches ~{barcode_mismatches} \
            ~{rc_flag} \
            --no-wait

        echo "[INFO] Disk after bcl_to_fastq:"
        df -h .

        # ── 4. Collect outputs ──────────────────────────────────────────────
        BASECALLS="run_dir/Data/Intensities/BaseCalls"

        # Move (not copy) FASTQs to avoid doubling disk usage.
        # bcl_to_fastq names them <sample>_R1.fastq.gz / <sample>_R2.fastq.gz
        mkdir -p fastqs_out

        find "${BASECALLS}" -maxdepth 1 -name "*_R1.fastq.gz" \
            ! -name "Undetermined*" \
            -exec mv {} fastqs_out/ \;

        find "${BASECALLS}" -maxdepth 1 -name "*_R2.fastq.gz" \
            ! -name "Undetermined*" \
            -exec mv {} fastqs_out/ \;

        # Move stats and log to the Cromwell working directory
        mv run_dir/demultiplexing_stats.csv demultiplexing_stats.csv
        mv run_dir/bcl2fastq.log            bcl2fastq.log

        echo "[INFO] Done. FASTQ files:"
        ls -lh fastqs_out/
        echo "[INFO] Final disk usage:"
        df -h .
    >>>

    output {
        Array[File] fastq_r1      = glob("fastqs_out/*_R1.fastq.gz")
        Array[File] fastq_r2      = glob("fastqs_out/*_R2.fastq.gz")
        File        demux_stats   = "demultiplexing_stats.csv"
        File        bcl2fastq_log = "bcl2fastq.log"
    }

    runtime {
        docker:      "quay.io/biocontainers/bcl2fastq-nextseq:1.3.0--pyh5e36f6f_0"
        memory:      "~{memory_gb} GiB"
        disks:       "local-disk ~{disk_gb} SSD"
        cpu:         processing_threads
        preemptible: preemptible
        maxRetries:  1
    }
}
