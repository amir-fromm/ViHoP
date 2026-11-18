#!/bin/bash
# Search candidate proteins against a general DIAMOND database and keep
# viral candidates: a best hit to Viruses, or no hit at all.
# Bacteria, Eukaryotes, Others, and Phages are written as false positives.
#
# Example, with config and summarize_diamonds.py next to this script:
#   diamond_to_general.sh --input-fasta candidates.faa

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ORIGINAL_ARGS=("$@")

usage() {
    cat <<EOF
Usage: diamond_to_general.sh [options]

Required, unless set in the config file:
  --input-fasta FILE         Candidate FASTA, protein or nucleotide
  --database FILE            DIAMOND database (.dmnd, or the name used with makedb)
  --catalog FILE             RefSeq taxid catalog (taxid, accession)
  --nodes FILE               taxid and NCBI division
  --merged FILE              old taxid and replacement taxid

Defaults:
  --working-dir DIR          Base directory. A relative --output-dir is created here
                             (default: the current directory)
  --output-dir DIR           Folder created for every result file, the log, and
                             the DIAMOND temporary files
                             (default: WORKING_DIR/mapped_to_db)
  --summarize-script FILE    summarize_diamonds.py (assigns NCBI divisions)
                             (default: ${SCRIPT_DIR}/summarize_diamonds.py)
  --config FILE              Shell file of KEY=value settings. Flags on the
                             command line override the file.
                             (default: ${SCRIPT_DIR}/config)
                             Copy config.example to config and fill in the paths.

Search:
  --diamond FILE             DIAMOND executable (default: diamond on PATH)
  --python FILE              Python executable (default: python)
  --threads N                DIAMOND threads (default: 16)
  --sensitivity MODE         fast, faster, mid-sensitive, sensitive,
                             more-sensitive, very-sensitive, ultra-sensitive
                             (default: very-sensitive)
  --max-target-seqs N        Best hits to keep per protein (default: 1)
  --evalue FLOAT             Maximum e-value (default: 0.001)
  --block-size FLOAT         DIAMOND sequence block size, billions of letters
  --tmpdir DIR               DIAMOND temporary directory
                             (default: OUTPUT_DIR/diamond_tmp)
  --chunks N                 Passed through to the summary script (default: 100000)

A protein FASTA is searched with blastp. A nucleotide FASTA is searched with
blastx, so the reads are translated. Letters E, F, I, L, P, Q, J, Z, O, or a
stop (*) mark the file as protein. Otherwise a sequence that is at least 90%
A, C, G, T, U, and N is nucleotide.

Cluster, optional:
  --cluster none|lsf         none runs DIAMOND in this shell (default: none).
                             lsf submits one job and returns. The compute
                             node runs the search and the summary.
  --queue NAME               LSF queue (default: short)
  --memory-mb N              LSF memory reserved per task, in MB (default: 8000)
  --tasks N                  LSF tasks, passed as bsub -n (default: --threads).
                             span[ptile] is set to this same N, so every slot
                             is on one host. DIAMOND does not run across hosts.

Underscore names are accepted too: --working_dir, --input_fasta, and so on.
An existing non-empty alignment table is reused and DIAMOND is not run again.
EOF
}

WORKING_DIR=""
OUTPUT_DIR=""
INPUT_FASTA=""
DATABASE=""
SUMMARIZE_SCRIPT=""
DIAMOND_BIN="diamond"
PYTHON_BIN="python"
THREADS=16
SENSITIVITY="very-sensitive"
MAX_TARGET_SEQS=1
EVALUE="0.001"
BLOCK_SIZE=""
TMPDIR_PATH=""
CHUNKS=100000
CATALOG=""
NODES=""
MERGED=""
CLUSTER="none"
QUEUE="short"
MEMORY_MB=8000
TASKS=""

CONFIG_FILE=""
CONFIG_SET=0
previous=""
for argument in "$@"; do
    if [ "$previous" = "--config" ]; then
        CONFIG_FILE=$argument
        CONFIG_SET=1
    fi
    previous=$argument
done
if [ "$CONFIG_SET" -eq 0 ]; then
    CONFIG_FILE="$SCRIPT_DIR/config"
fi
if [ -f "$CONFIG_FILE" ]; then
    # shellcheck source=config.example
    source "$CONFIG_FILE"
elif [ "$CONFIG_SET" -eq 1 ]; then
    echo "Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --working-dir|--working_dir) WORKING_DIR=$2; shift 2 ;;
        --output-dir|--output_dir) OUTPUT_DIR=$2; shift 2 ;;
        --input-fasta|--input_fasta) INPUT_FASTA=$2; shift 2 ;;
        --database) DATABASE=$2; shift 2 ;;
        --summarize-script|--summarize_script) SUMMARIZE_SCRIPT=$2; shift 2 ;;
        --diamond) DIAMOND_BIN=$2; shift 2 ;;
        --python) PYTHON_BIN=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --sensitivity) SENSITIVITY=$2; shift 2 ;;
        --max-target-seqs|--max_target_seqs) MAX_TARGET_SEQS=$2; shift 2 ;;
        --evalue) EVALUE=$2; shift 2 ;;
        --block-size|--block_size) BLOCK_SIZE=$2; shift 2 ;;
        --tmpdir) TMPDIR_PATH=$2; shift 2 ;;
        --chunks) CHUNKS=$2; shift 2 ;;
        --catalog) CATALOG=$2; shift 2 ;;
        --nodes) NODES=$2; shift 2 ;;
        --merged) MERGED=$2; shift 2 ;;
        --config) shift 2 ;;
        --cluster) CLUSTER=$2; shift 2 ;;
        --queue) QUEUE=$2; shift 2 ;;
        --memory-mb|--memory_mb) MEMORY_MB=$2; shift 2 ;;
        --tasks) TASKS=$2; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [ -z "$WORKING_DIR" ]; then
    WORKING_DIR=$(pwd)
fi
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$WORKING_DIR/mapped_to_db"
elif [ "${OUTPUT_DIR#/}" = "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$WORKING_DIR/$OUTPUT_DIR"
fi
OUTPUT_DIR=${OUTPUT_DIR%/}
if [ -z "$SUMMARIZE_SCRIPT" ]; then
    SUMMARIZE_SCRIPT="$SCRIPT_DIR/summarize_diamonds.py"
fi

if [ -z "$INPUT_FASTA" ] || [ -z "$DATABASE" ] || [ -z "$CATALOG" ] || [ -z "$NODES" ] || [ -z "$MERGED" ]; then
    echo "Missing a required option." >&2
    usage >&2
    exit 1
fi

case "$SENSITIVITY" in
    fast|faster|mid-sensitive|sensitive|more-sensitive|very-sensitive|ultra-sensitive) ;;
    *)
        echo "Unknown sensitivity: $SENSITIVITY" >&2
        exit 1
        ;;
esac

case "$CLUSTER" in
    none|lsf) ;;
    *)
        echo "Unknown cluster: $CLUSTER (use none or lsf)" >&2
        exit 1
        ;;
esac

if [ -z "$TASKS" ]; then
    TASKS=$THREADS
fi
if [ -z "$TMPDIR_PATH" ]; then
    TMPDIR_PATH="$OUTPUT_DIR/diamond_tmp"
fi

LOG_DIR="$OUTPUT_DIR/log"
ALIGNMENT="$OUTPUT_DIR/diamond_candidates_all.tsv"
EDITED="$OUTPUT_DIR/edited_diamond_all.tsv"
FALSE_POSITIVES="$OUTPUT_DIR/false_positives.list"
SUMMARY="$OUTPUT_DIR/summary_diamond.tsv"
NO_MATCHES="$OUTPUT_DIR/no_matches.list"
echo "Output directory: $OUTPUT_DIR"

require_file() {
    if [ ! -s "$1" ]; then
        echo "Missing or empty file: $1" >&2
        exit 1
    fi
}

require_file "$INPUT_FASTA"
require_file "$SUMMARIZE_SCRIPT"
require_file "$CATALOG"
require_file "$NODES"
require_file "$MERGED"
if [ ! -x "$DIAMOND_BIN" ] && ! command -v "$DIAMOND_BIN" >/dev/null 2>&1; then
    echo "DIAMOND executable not found: $DIAMOND_BIN" >&2
    exit 1
fi
if [ ! -s "$DATABASE" ] && [ ! -s "${DATABASE}.dmnd" ]; then
    echo "DIAMOND database not found: $DATABASE" >&2
    exit 1
fi

mkdir -p "$LOG_DIR" "$TMPDIR_PATH"

fasta_molecule() {
    awk '
        BEGIN { protein_only = 0; nucleotide = 0; seen = 0; limit = 200000 }
        /^>/ { next }
        {
            sequence = toupper($0)
            length_sequence = length(sequence)
            for (i = 1; i <= length_sequence; i++) {
                letter = substr(sequence, i, 1)
                if (letter !~ /[A-Z*]/) continue
                seen++
                if (index("EFILPQJZO*", letter)) protein_only++
                else if (index("ACGTUN", letter)) nucleotide++
                if (seen >= limit) break
            }
            if (seen >= limit) exit
        }
        END {
            if (seen == 0) {
                print "empty"
            } else if (protein_only > 0) {
                print "protein"
            } else if (nucleotide / seen >= 0.9) {
                print "nucleotide"
            } else {
                print "protein"
            }
        }
    ' "$1"
}

MOLECULE=$(fasta_molecule "$INPUT_FASTA")
case "$MOLECULE" in
    protein)
        DIAMOND_COMMAND=blastp
        echo "Input is protein. Running a direct DIAMOND blastp search."
        ;;
    nucleotide)
        DIAMOND_COMMAND=blastx
        echo "Input is nucleotide. Running a translated DIAMOND blastx search."
        ;;
    *)
        echo "No sequence letters found in $INPUT_FASTA" >&2
        exit 1
        ;;
esac

run_diamond() {
    local -a command
    command=(
        "$DIAMOND_BIN" "$DIAMOND_COMMAND"
        --db "$DATABASE"
        --query "$INPUT_FASTA"
        --out "$ALIGNMENT"
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        --threads "$THREADS"
        --max-target-seqs "$MAX_TARGET_SEQS"
        --max-hsps 1
        --evalue "$EVALUE"
        --tmpdir "$TMPDIR_PATH"
        "--${SENSITIVITY}"
        --verbose
    )
    if [ -n "$BLOCK_SIZE" ]; then
        command+=(--block-size "$BLOCK_SIZE")
    fi

    if [ "$CLUSTER" = "none" ]; then
        echo "Running DIAMOND $DIAMOND_COMMAND"
        "${command[@]}"
        return
    fi

    local runner="$LOG_DIR/run_pipeline.sh"
    local job_file="$LOG_DIR/diamond.jobid"
    if [ -s "$job_file" ]; then
        local previous
        previous=$(cat "$job_file")
        if bjobs -noheader "$previous" 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; then
            echo "LSF job $previous is already queued for this output directory."
            echo "This session does not wait. Results will be written to $OUTPUT_DIR"
            exit 0
        fi
    fi

    local -a forwarded=()
    local skip_next=0
    local arg
    for arg in "${ORIGINAL_ARGS[@]}"; do
        if [ "$skip_next" -eq 1 ]; then
            skip_next=0
            continue
        fi
        case "$arg" in
            --cluster) skip_next=1 ;;
            --cluster=*) ;;
            *) forwarded+=("$arg") ;;
        esac
    done

    {
        printf '#!/bin/bash\nset -euo pipefail\n'
        printf '%q ' "$SCRIPT_DIR/diamond_to_general.sh" --cluster none "${forwarded[@]}"
        printf '\n'
    } > "$runner"
    chmod +x "$runner"

    echo "Submitting DIAMOND to LSF queue $QUEUE"
    local job_id
    job_id=$(bsub -q "$QUEUE" \
        -J "$LOG_DIR/diamond_refseq." \
        -o "$LOG_DIR/diamond_refseq.%J.o" \
        -e "$LOG_DIR/diamond_refseq.%J.e" \
        -n "$TASKS" \
        -R "rusage[mem=${MEMORY_MB}] span[ptile=${TASKS}]" \
        "$runner" | awk '/^Job </ { gsub(/[<>]/, "", $2); print $2 }')
    if [ -z "$job_id" ]; then
        echo "LSF did not return a job id." >&2
        exit 1
    fi
    echo "$job_id" > "$job_file"
    echo "LSF job $job_id submitted. This session does not wait."
    echo "The search and the summary run on the compute node. Results: $OUTPUT_DIR"
    exit 0
}

if [ -s "$ALIGNMENT" ]; then
    echo "Found existing DIAMOND results ($ALIGNMENT); skipping the search"
else
    run_diamond
fi

if [ ! -f "$ALIGNMENT" ]; then
    echo "DIAMOND did not write $ALIGNMENT" >&2
    exit 1
fi

if [ ! -s "$ALIGNMENT" ]; then
    echo "DIAMOND produced no alignments. Every input protein is recorded as no hit."
    : > "$FALSE_POSITIVES"
    : > "$SUMMARY"
else
    line_count=$(wc -l < "$ALIGNMENT")
    echo "Assigning a division to the best hit of each protein"
    summarize_cmd=(
        "$PYTHON_BIN" "$SUMMARIZE_SCRIPT"
        --chunks "$CHUNKS"
        --size_table "$line_count"
        --output_path "$EDITED"
        --alignment "$ALIGNMENT"
        --catalog "$CATALOG"
        --nodes "$NODES"
        --merged "$MERGED"
    )
    "${summarize_cmd[@]}"

    awk -F '\t' '
        !/^#/ && ($14 == "Bacteria" || $14 == "Eukaryotes" || $14 == "Others" || $14 == "Phages") {
            if (!seen[$1]++) print $1
        }
    ' "$EDITED" > "$FALSE_POSITIVES"

    awk -F '\t' '
        !/^#/ && $14 == "Viruses" {
            OFS = "\t"
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
        }
    ' "$EDITED" > "$SUMMARY"
fi

input_ids=$(mktemp)
matched_ids=$(mktemp)
awk '/^>/ { id = substr($1, 2); if (id != "") print id }' "$INPUT_FASTA" | sort -u > "$input_ids"
awk -F '\t' '!/^#/ && $1 != "" { print $1 }' "$ALIGNMENT" | sort -u > "$matched_ids"
comm -23 "$input_ids" "$matched_ids" > "$NO_MATCHES"
awk -F '\t' 'BEGIN { OFS = "\t" } { print $0, "no hit", "", "", "", "", "", "", "", "", "", "" }' \
    "$NO_MATCHES" >> "$SUMMARY"
echo "Proteins with no match: $(wc -l < "$NO_MATCHES")"
echo "False positives: $(wc -l < "$FALSE_POSITIVES")"
echo "Viral candidates (virus hits plus no hits): $SUMMARY"
rm -f "$input_ids" "$matched_ids"
