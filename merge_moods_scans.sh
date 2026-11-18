#!/usr/bin/env bash
# merge_moods.sh: merge sorted per-motif MOODS bed.gz files into ONE bgzipped,
# tabix-indexed BED, using bounded two-level k-way merges (no re-sort).
#
# Run it on the login node; it submits itself to Slurm with the resources you ask for:
#   ./merge_moods.sh \
#       --input-dir /net/seq/data2/projects/sabramov/ENCODE4/motif_scans.p1e_minus_5/mouse/output/moods_scans_ref \
#       --output mouse.moods_all.bed.gz --cpus 16 --mem 16G --time 12:00:00
# Inside an existing allocation (srun --pty bash, or sbatch) it runs in place instead.
# See ./merge_moods.sh --help for every option.

set -euo pipefail
export LC_ALL=C

usage() {
    cat <<'EOF'
Usage: merge_moods.sh --input-dir DIR --output FILE [options]

Data
  --input-dir DIR       Folder containing the per-motif .bed.gz files (required)
  --output FILE         Final .bed.gz path; its .tbi is written next to it (required)
  --pattern GLOB        Input filename pattern (default: *.moods.log.bed.gz)

Slurm
  --partition NAME      (default: hpcz-test)
  --cpus N              CPUs for the single task (default: 16)
  --mem SIZE            (default: 16G)
  --time D-HH:MM:SS     (default: 12:00:00)
  --job-name NAME       (default: merge_<output name>)
  --log FILE            Slurm log; %j = job id (default: ./merge_<output name>_%j.log)
  --account NAME, --qos NAME
  --sbatch-opt OPT      Any extra sbatch option, repeatable (e.g. --sbatch-opt=--mail-type=END)
  --dry-run             Print the sbatch command instead of submitting
  --local               Run in the current shell, never submit
  --submit              Submit even when already inside a Slurm allocation

Tuning
  --workers N           Concurrent stage-1 batches (default: 3/8 of CPUs, min 1)
  --batch-size N        Files per stage-1 batch (default: 48)
  --bgzip-threads N     (default: half of CPUs, min 1)
  --workdir DIR         Scratch for intermediates (default: $SLURM_TMPDIR, $TMPDIR or /tmp)
  --limit N             Test run on the first N files; output gets a .testN suffix
  --no-check-input      Skip the per-file sort/integrity check
  --no-check-counts     Skip the input/output record-count comparison
EOF
}

# ---------------------------------------------------------------- arguments
IN=""; OUTPUT=""
PATTERN="*.moods.log.bed.gz"
PARTITION="hpcz-test"; CPUS=16; MEM="16G"; TIME="12:00:00"
JOB_NAME=""; LOG=""; ACCOUNT=""; QOS=""; SBATCH_OPTS=()
MODE="auto"; DRY_RUN=0
WORKERS=""; BATCH_SIZE=48; BGZIP_THREADS=""; WORKDIR=""; LIMIT=""
CHECK_INPUT=1; CHECK_COUNTS=1

need() { [[ $# -ge 2 && -n "${2:-}" ]] || { echo "Option $1 needs a value" >&2; exit 2; }; }
while (( $# )); do
    case "$1" in
        --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}"; continue ;;
        --input-dir)     need "$@"; IN=$2; shift 2 ;;
        --output)        need "$@"; OUTPUT=$2; shift 2 ;;
        --pattern)       need "$@"; PATTERN=$2; shift 2 ;;
        --partition)     need "$@"; PARTITION=$2; shift 2 ;;
        --cpus)          need "$@"; CPUS=$2; shift 2 ;;
        --mem)           need "$@"; MEM=$2; shift 2 ;;
        --time)          need "$@"; TIME=$2; shift 2 ;;
        --job-name)      need "$@"; JOB_NAME=$2; shift 2 ;;
        --log)           need "$@"; LOG=$2; shift 2 ;;
        --account)       need "$@"; ACCOUNT=$2; shift 2 ;;
        --qos)           need "$@"; QOS=$2; shift 2 ;;
        --sbatch-opt)    need "$@"; SBATCH_OPTS+=("$2"); shift 2 ;;
        --workers)       need "$@"; WORKERS=$2; shift 2 ;;
        --batch-size)    need "$@"; BATCH_SIZE=$2; shift 2 ;;
        --bgzip-threads) need "$@"; BGZIP_THREADS=$2; shift 2 ;;
        --workdir)       need "$@"; WORKDIR=$2; shift 2 ;;
        --limit)         need "$@"; LIMIT=$2; shift 2 ;;
        --no-check-input)  CHECK_INPUT=0; shift ;;
        --no-check-counts) CHECK_COUNTS=0; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --local)   MODE="local"; shift ;;
        --submit)  MODE="submit"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done

[[ -n "$IN" ]] || { echo "--input-dir is required (see --help)" >&2; exit 2; }
[[ -n "$OUTPUT" ]] || { echo "--output is required (see --help)" >&2; exit 2; }
is_pos_int() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
for pair in "cpus:$CPUS" "batch-size:$BATCH_SIZE" "workers:${WORKERS:-1}" \
            "bgzip-threads:${BGZIP_THREADS:-1}" "limit:${LIMIT:-1}"; do
    is_pos_int "${pair#*:}" || { echo "--${pair%%:*} must be a positive integer" >&2; exit 2; }
done

[[ -d "$IN" ]] || { echo "Missing input directory: $IN" >&2; exit 1; }
IN="$(realpath -e -- "$IN")"
OUTPUT="$(realpath -m -- "$OUTPUT")"
[[ "$OUTPUT" == *.gz ]] || { echo "--output must end in .gz (tabix needs bgzip)" >&2; exit 2; }
[[ "$OUTPUT" != "$IN"/* ]] || { echo "--output must not be inside --input-dir" >&2; exit 2; }
# Short, filesystem-safe label from the output name (used for job name, log and scratch dir).
NAME="$(basename "$OUTPUT")"; NAME="${NAME%.gz}"; NAME="${NAME%.bed}"
NAME="$(printf '%s' "$NAME" | tr -c 'A-Za-z0-9._-' '_')"
[[ -z "$WORKDIR" ]] || WORKDIR="$(realpath -m -- "$WORKDIR")"

# ---------------------------------------------------------------- submit mode
if [[ "$MODE" == submit || ( "$MODE" == auto && -z "${SLURM_JOB_ID:-}" ) ]]; then
    command -v sbatch >/dev/null || { echo "sbatch not found; use --local" >&2; exit 1; }
    JOB_NAME="${JOB_NAME:-merge_${NAME}}"
    LOG="$(realpath -m -- "${LOG:-$PWD/merge_${NAME}_%j.log}")"
    mkdir -p "$(dirname "$LOG")"

    # Settings handed to the job, with all paths already absolute.
    JOB_ARGS=(--local --input-dir "$IN" --output "$OUTPUT"
              --pattern "$PATTERN" --cpus "$CPUS" --batch-size "$BATCH_SIZE")
    if [[ -n "$WORKERS" ]]; then JOB_ARGS+=(--workers "$WORKERS"); fi
    if [[ -n "$BGZIP_THREADS" ]]; then JOB_ARGS+=(--bgzip-threads "$BGZIP_THREADS"); fi
    if [[ -n "$WORKDIR" ]]; then JOB_ARGS+=(--workdir "$WORKDIR"); fi
    if [[ -n "$LIMIT" ]]; then JOB_ARGS+=(--limit "$LIMIT"); fi
    if [[ "$CHECK_INPUT" == 0 ]]; then JOB_ARGS+=(--no-check-input); fi
    if [[ "$CHECK_COUNTS" == 0 ]]; then JOB_ARGS+=(--no-check-counts); fi

    # One node, one task: every step runs on a single machine, so extra nodes would sit idle.
    SB=(sbatch --job-name="$JOB_NAME" --partition="$PARTITION" --nodes=1 --ntasks=1
        --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME" --output="$LOG" --export=ALL)
    if [[ -n "$ACCOUNT" ]]; then SB+=(--account="$ACCOUNT"); fi
    if [[ -n "$QOS" ]]; then SB+=(--qos="$QOS"); fi
    SB+=(${SBATCH_OPTS[@]+"${SBATCH_OPTS[@]}"})
    SB+=("$(realpath -e -- "${BASH_SOURCE[0]}")" "${JOB_ARGS[@]}")

    if (( DRY_RUN )); then
        printf '%q ' "${SB[@]}"; echo
        exit 0
    fi
    "${SB[@]}"
    echo "Log: $LOG   (watch with: squeue -u ${USER:-$(id -un)})"
    exit 0
fi

# ---------------------------------------------------------------- run mode
NCPU="${SLURM_CPUS_PER_TASK:-$CPUS}"
WORKERS="${WORKERS:-$(( NCPU * 3 / 8 > 0 ? NCPU * 3 / 8 : 1 ))}"
BGZIP_THREADS="${BGZIP_THREADS:-$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))}"
(( WORKERS <= NCPU )) || { echo "--workers exceeds allocated CPUs ($NCPU)" >&2; exit 2; }

for tool in gzip zstd bgzip tabix sort split xargs mkfifo wc df awk; do
    command -v "$tool" >/dev/null || {
        echo "Missing $tool. Load your cluster's relevant modules/conda environment." >&2
        exit 1
    }
done

if [[ -n "$LIMIT" ]]; then
    if [[ "$OUTPUT" == *.bed.gz ]]; then FINAL="${OUTPUT%.bed.gz}.test${LIMIT}.bed.gz"
    else FINAL="${OUTPUT%.gz}.test${LIMIT}.gz"; fi
else
    FINAL="$OUTPUT"
fi
OUT_DIR="$(dirname "$FINAL")"
mkdir -p "$OUT_DIR"
if [[ -e "$FINAL" || -e "$FINAL.tbi" ]]; then
    echo "Output already exists: $FINAL (move/delete it before rerunning)" >&2
    exit 1
fi

# Cluster-specific SLURM_TMPDIR takes priority; /tmp is only a fallback.
SCRATCH_BASE="${WORKDIR:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}"
mkdir -p "$SCRATCH_BASE"
WORK="$(mktemp -d "${SCRATCH_BASE%/}/merge_${NAME}_${SLURM_JOB_ID:-$$}_XXXXXX")"
mkdir -p "$WORK/lists" "$WORK/chunks" "$WORK/fifos" "$WORK/counts"
TEMP_FINAL="$OUT_DIR/.$(basename "$FINAL" .gz).${SLURM_JOB_ID:-$$}.tmp.gz"
cleanup() {
    rm -f "$TEMP_FINAL" "$TEMP_FINAL.tbi"
    rm -rf "$WORK"
}
trap cleanup EXIT
# Slurm sends SIGTERM at the time limit; bash skips the EXIT trap unless the signal is trapped.
trap 'exit 143' TERM
trap 'exit 130' INT

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
log "Host=$(hostname); job=${SLURM_JOB_ID:-none}"
log "Name=$NAME; CPUs=$NCPU; concurrent batches=$WORKERS; batch size=$BATCH_SIZE; bgzip threads=$BGZIP_THREADS"
log "Input=$IN ($PATTERN)"
log "Output=$FINAL; scratch=$WORK"
df -h "$WORK"
if [[ "$SCRATCH_BASE" == /tmp ]]; then
    log "NOTE: using /tmp; it must have room for roughly the compressed input size."
fi

# Input filenames in this project have no newlines.
shopt -s nullglob
# shellcheck disable=SC2206  # PATTERN is intentionally glob-expanded
FILES=("$IN"/$PATTERN)
(( ${#FILES[@]} > 0 )) || { echo "No input files match $IN/$PATTERN" >&2; exit 1; }
if [[ -n "$LIMIT" ]]; then
    FILES=("${FILES[@]:0:LIMIT}")
fi
printf '%s\n' "${FILES[@]}" > "$WORK/inputs.txt"
log "Found ${#FILES[@]} input motif files"

# Each sort process opens all of its inputs at once: BATCH_SIZE in stage 1, NCHUNKS in stage 2.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
NCHUNKS=$(( (${#FILES[@]} + BATCH_SIZE - 1) / BATCH_SIZE ))
NEED_FD=$(( (NCHUNKS > BATCH_SIZE ? NCHUNKS : BATCH_SIZE) + 32 ))
FD_LIMIT="$(ulimit -n)"
if [[ "$FD_LIMIT" != unlimited ]] && (( FD_LIMIT < NEED_FD )); then
    echo "Open-file limit $FD_LIMIT < needed $NEED_FD; raise --batch-size or the limit" >&2
    exit 1
fi
log "Open-file limit=$FD_LIMIT (need ~$NEED_FD); expected stage-1 chunks=$NCHUNKS"

# IMPORTANT: this algorithm merges; it does not repair unsorted input.
# One parallel pass checks gzip integrity + sort order (same comparator as the merges:
# chrom as C-locale string, then start, then end, numerically) and records each file's
# line count, so the final completeness check needs no second read of the inputs.
CHECK_AWK='
BEGIN { FS = "\t" }
NR > 1 {
    if ($1 < c) bad = 1
    else if ($1 == c) {
        if ($2 + 0 < s) bad = 1
        else if ($2 + 0 == s && $3 + 0 < e) bad = 1
    }
    if (bad) { printf "unsorted at line %d\n", NR > "/dev/stderr"; exit 1 }
}
{ c = $1; s = $2 + 0; e = $3 + 0 }
END { if (!bad) print NR }'
export CHECK_AWK WORK

if [[ "$CHECK_INPUT" == 1 ]]; then
    log "Checking sorted order + gzip integrity of ALL selected files (and counting records)"
    check_one() {
        set -euo pipefail
        local n
        if ! n=$(gzip -dc -- "$1" | awk "$CHECK_AWK"); then
            echo "ERROR: corrupt or unsorted input: $1" >&2
            return 1
        fi
        printf '%s\n' "$n" > "$WORK/counts/${1##*/}.n"
    }
    export -f check_one
    printf '%s\0' "${FILES[@]}" |
        xargs -0 -r -n 1 -P "$NCPU" bash -c 'check_one "$1"' _
    log "All selected inputs passed the sort-order check"
else
    log "WARNING: --no-check-input; merge assumes the inputs are sorted"
fi

# Launch <= BATCH_SIZE gzip producers and merge them to one .zst chunk.
# Background exit codes are explicitly checked (not a bare `wait`).
merge_batch() {
    set -euo pipefail
    local batch_list="$1" name fifo_dir f fifo pid rc=0 i=0
    local -a paths=() pids=() pipes=()
    name="${batch_list##*/}"
    fifo_dir="$WORK/fifos/$name"
    mkdir -p "$fifo_dir"
    mapfile -t paths < "$batch_list"
    (( ${#paths[@]} > 0 )) || { echo "Empty batch $name" >&2; return 1; }

    # This is a separate xargs-launched bash, so EXIT cleanup is local to this worker.
    cleanup_worker() {
        local child
        for child in "${pids[@]}"; do kill "$child" 2>/dev/null || true; done
        wait 2>/dev/null || true
        rm -rf "$fifo_dir"
    }
    trap cleanup_worker EXIT
    for f in "${paths[@]}"; do
        fifo="$fifo_dir/$i"
        mkfifo "$fifo"
        gzip -dc -- "$f" > "$fifo" &
        pids+=("$!")
        pipes+=("$fifo")
        i=$((i+1))
    done

    sort --merge --stable --batch-size="$(( ${#pipes[@]} > 1 ? ${#pipes[@]} : 2 ))" \
        -k1,1 -k2,2n -k3,3n "${pipes[@]}" |
        zstd -q -1 --single-thread -c > "$WORK/chunks/$name.zst" || rc=1

    if (( rc != 0 )); then
        for pid in "${pids[@]}"; do kill "$pid" 2>/dev/null || true; done
    fi
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then rc=1; fi
    done
    trap - EXIT
    rm -rf "$fifo_dir"
    (( rc == 0 )) || { echo "Batch failed: $name" >&2; return 1; }
}
export -f merge_batch

split -l "$BATCH_SIZE" -d -a 4 "$WORK/inputs.txt" "$WORK/lists/batch_"
log "Stage 1: $(find "$WORK/lists" -type f | wc -l) batches"
find "$WORK/lists" -name 'batch_*' -type f -print0 |
    sort -z |
    xargs -0 -r -n 1 -P "$WORKERS" bash -c 'merge_batch "$1"' _
log "Stage 1 completed"
du -sh "$WORK/chunks"

# Stage 2: merge sorted .zst chunks directly into ONE BGZF stream.
CHUNKS=("$WORK"/chunks/batch_*.zst)
(( ${#CHUNKS[@]} > 0 )) || { echo "No intermediate chunks" >&2; exit 1; }
log "Stage 2: merging ${#CHUNKS[@]} compressed chunks"
mkdir -p "$WORK/fifos/final"
PIDS=()
PIPES=()
for i in "${!CHUNKS[@]}"; do
    fifo="$WORK/fifos/final/$i"
    mkfifo "$fifo"
    zstd -q -dc -- "${CHUNKS[$i]}" > "$fifo" &
    PIDS+=("$!")
    PIPES+=("$fifo")
done

rc=0
sort --merge --stable --batch-size="$(( ${#PIPES[@]} > 1 ? ${#PIPES[@]} : 2 ))" \
    -k1,1 -k2,2n -k3,3n "${PIPES[@]}" |
    bgzip -@ "$BGZIP_THREADS" -c > "$TEMP_FINAL" || rc=1
if (( rc != 0 )); then
    for pid in "${PIDS[@]}"; do kill "$pid" 2>/dev/null || true; done
fi
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then rc=1; fi
done
(( rc == 0 )) || { echo "Final merge/decompression failed" >&2; exit 1; }
log "Stage 2 completed"

# One multithreaded read of the output verifies BGZF integrity and counts records.
# Sort order is enforced by tabix below (it refuses unsorted or split chromosome blocks).
log "Validating final BGZF (integrity + record count)"
count_out=$(bgzip -@ "$BGZIP_THREADS" -dc "$TEMP_FINAL" | wc -l)
log "Output records=$count_out"
if [[ "$CHECK_COUNTS" == 1 ]]; then
    if [[ "$CHECK_INPUT" == 1 ]]; then
        count_in=$(cat "$WORK"/counts/*.n | awk '{ s += $1 } END { print s + 0 }')
    else
        log "--no-check-input was set, so counting input records with a serial read"
        count_in=$(gzip -dc -- "${FILES[@]}" | wc -l)
    fi
    log "Input records=$count_in; output records=$count_out"
    [[ "$count_in" -eq "$count_out" ]] || { echo "ERROR: missing/extra records!" >&2; exit 1; }
fi

log "Building tabix index"
tabix -p bed "$TEMP_FINAL"
[[ -s "$TEMP_FINAL.tbi" ]] || { echo "Tabix index missing/empty" >&2; exit 1; }

mv "$TEMP_FINAL" "$FINAL"
mv "$TEMP_FINAL.tbi" "$FINAL.tbi"
log "SUCCESS: $FINAL"
ls -lh "$FINAL" "$FINAL.tbi"
