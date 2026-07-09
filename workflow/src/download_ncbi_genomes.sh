#!/bin/bash

################################################################################
# Script modified from Claude.ai generated code
# NCBI Genome Download Script
# Downloads genome sequences, annotations, and protein sequences for a list of
# accessions into a shared genome database, skipping already-downloaded genomes
################################################################################

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
print_warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }
print_debug() { echo -e "${BLUE}[DEBUG]${NC} $1"; }

# ── Arguments ─────────────────────────────────────────────────────────────────

ACCESSIONS_FILE="${1:-}"
GENOME_DB="${2:-}"
GENOME_LOG="${3:-}"
AVAILABILITY="${4:-}"
FORMATS="${5:-gff3,genome,protein}"
BATCH_SIZE="${6:-50}"
TMP_BASE="${7:-.}"

# ── Validate arguments ────────────────────────────────────────────────────────

validate_args() {
    if [[ -z "$ACCESSIONS_FILE" || -z "$GENOME_DB" || -z "$GENOME_LOG" || -z "$AVAILABILITY" ]]; then
        print_error "Missing required arguments"
        echo ""
        echo "Usage: $0 <accessions_file> <genome_db> <genome_log> <availability> [formats] [batch_size] [tmp_dir]"
        echo ""
        echo "  accessions_file   Path to accessions file (plain text or JSON array)"
        echo "  genome_db         Path to shared genome database directory"
        echo "  genome_log        Path to genome download log JSON file"
        echo "  availability      Path to output availability JSON file"
        echo "  formats           Comma-separated download formats (default: gff3,genome,protein)"
        echo "  batch_size        Accessions per download batch (default: 50)"
        echo "  tmp_dir           Directory for temporary files (default: .)"
        exit 1
    fi

    if [[ ! -f "$ACCESSIONS_FILE" ]]; then
        print_error "Accessions file not found: $ACCESSIONS_FILE"
        exit 1
    fi
}

# ── Temp directory setup ──────────────────────────────────────────────────────

setup_tmpdir() {
    TMPDIR=$(mktemp -d "${TMP_BASE}/genome_download_XXXXXX")
    trap 'rm -rf "$TMPDIR"' EXIT
    print_info "Temp directory: $TMPDIR"
}

# ── Parse accessions ──────────────────────────────────────────────────────────

parse_accessions() {
    if jq empty "$ACCESSIONS_FILE" 2>/dev/null; then
        jq -r '.[]' "$ACCESSIONS_FILE" > "$TMPDIR/all_accessions.txt"
    else
        grep -v '^#' "$ACCESSIONS_FILE" | grep -v '^[[:space:]]*$' | sed 's/[[:space:]]*$//' > "$TMPDIR/all_accessions.txt"
    fi

    local count
    count=$(wc -l < "$TMPDIR/all_accessions.txt")
    if [[ "$count" -eq 0 ]]; then
        print_error "No accessions found in $ACCESSIONS_FILE"
        exit 1
    fi
    print_info "Parsed $count accessions"
}

# ── Partition into already-available vs to-download ───────────────────────────

partition_accessions() {
    > "$TMPDIR/available.txt"
    > "$TMPDIR/to_download.txt"

    while IFS= read -r accession; do
        local acc_dir="$GENOME_DB/$accession"
        if [[ -f "$acc_dir/genomic.fna" && -f "$acc_dir/geneannotation.gff" && -f "$acc_dir/protein.faa" ]]; then
            echo "$accession" >> "$TMPDIR/available.txt"
        else
            echo "$accession" >> "$TMPDIR/to_download.txt"
        fi
    done < "$TMPDIR/all_accessions.txt"

    local n_available n_to_download
    n_available=$(wc -l < "$TMPDIR/available.txt")
    n_to_download=$(wc -l < "$TMPDIR/to_download.txt")
    print_info "Already in database: $n_available | To download: $n_to_download"
}

# ── Download a single batch ───────────────────────────────────────────────────

download_batch() {
    local batch_file="$1"
    local batch_num="$2"
    local max_retries=3
    local retry_delay=10
    local batch_zip="$TMPDIR/batch_${batch_num}.zip"
    local batch_dir="$TMPDIR/batch_${batch_num}"

    local attempt=1
    while [[ $attempt -le $max_retries ]]; do
        print_info "Downloading batch $batch_num (attempt $attempt/$max_retries)..."
        if datasets download genome accession \
            --inputfile "$batch_file" \
            --include "$FORMATS" \
            --filename "$batch_zip" \
            --no-progressbar 2>&1; then
            break
        fi

        if [[ $attempt -eq $max_retries ]]; then
            print_error "Batch $batch_num failed after $max_retries attempts"
            return 1
        fi

        print_warn "Batch $batch_num failed, retrying in ${retry_delay}s..."
        sleep "$retry_delay"
        attempt=$((attempt + 1))
    done

    unzip -q -o "$batch_zip" -d "$batch_dir"
    rm "$batch_zip"

    for acc_dir in "$batch_dir/ncbi_dataset/data"/*/; do
        [[ -d "$acc_dir" ]] || continue
        local accession
        accession=$(basename "$acc_dir")
        local target_dir="$GENOME_DB/$accession"
        mkdir -p "$target_dir"

        find "$acc_dir" -name "*_genomic.fna"  -exec cp {} "$target_dir/genomic.fna" \;
        find "$acc_dir" -name "*.gff"          -exec cp {} "$target_dir/geneannotation.gff" \;
        find "$acc_dir" -name "*protein.faa"   -exec cp {} "$target_dir/protein.faa" \;

        if [[ -f "$target_dir/genomic.fna" && -f "$target_dir/geneannotation.gff" && -f "$target_dir/protein.faa" ]]; then
            echo "$accession" >> "$TMPDIR/newly_downloaded.txt"
            print_debug "Stored: $accession"
        else
            print_warn "Incomplete download for $accession, removing"
            rm -rf "$target_dir"
        fi
    done

    rm -rf "$batch_dir"
}

# ── Download all pending accessions ──────────────────────────────────────────

download_all() {
    local n_to_download
    n_to_download=$(wc -l < "$TMPDIR/to_download.txt")

    if [[ "$n_to_download" -eq 0 ]]; then
        print_info "Nothing to download"
        return
    fi

    local n_batches=$(( (n_to_download + BATCH_SIZE - 1) / BATCH_SIZE ))
    print_info "Splitting $n_to_download accessions into $n_batches batches"
    split -l "$BATCH_SIZE" -d -a 3 "$TMPDIR/to_download.txt" "$TMPDIR/batch_"

    local failed=0
    local batch_num=0

    for batch_file in "$TMPDIR"/batch_[0-9]*; do
        [[ -f "$batch_file" ]] || continue
        batch_num=$((batch_num + 1))

        if ! download_batch "$batch_file" "$batch_num"; then
            failed=$((failed + 1))
        fi

        sleep 2
    done

    print_info "Batches complete: $((batch_num - failed)) succeeded, $failed failed"
}

# ── Write genome download log ─────────────────────────────────────────────────

write_genome_log() {
    if [[ ! -f "$TMPDIR/newly_downloaded.txt" ]]; then
        print_info "No new genomes downloaded, genome log unchanged"
        return
    fi

    local log_content="{}"
    if [[ -f "$GENOME_LOG" ]]; then
        log_content=$(cat "$GENOME_LOG")
    fi

    local download_date
    download_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    while IFS= read -r accession; do
        log_content=$(echo "$log_content" | jq \
            --arg acc "$accession" \
            --arg ts  "$download_date" \
            --arg src "NCBI" \
            --arg fmt "$FORMATS" \
            '.[$acc] = {"downloaded": $ts, "source": $src, "formats": $fmt, "status": "success"}')
    done < "$TMPDIR/newly_downloaded.txt"

    echo "$log_content" | jq '.' > "$GENOME_LOG"
    print_info "Updated genome log: $GENOME_LOG"
}

# ── Write availability file ───────────────────────────────────────────────────

write_availability() {
    local all_available
    all_available=$(cat "$TMPDIR/available.txt" \
            "${TMPDIR}/newly_downloaded.txt" 2>/dev/null | sort -u | jq -R . | jq -s .)

    local analysis_date
    analysis_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    jq -n \
        --arg  date "$analysis_date" \
        --arg  db   "$GENOME_DB" \
        --argjson accs "$all_available" \
        '{"analysis_date": $date, "genome_database": $db, "accessions": $accs}' \
        > "$AVAILABILITY"

    print_info "Wrote availability: $AVAILABILITY"
}

# ── Main ──────────────────────────────────────────────────────────────────────

main() {
    validate_args
    setup_tmpdir

    print_info "Accessions file:  $ACCESSIONS_FILE"
    print_info "Genome database:  $GENOME_DB"
    print_info "Genome log:       $GENOME_LOG"
    print_info "Availability:     $AVAILABILITY"
    print_info "Formats:          $FORMATS"
    print_info "Batch size:       $BATCH_SIZE"

    parse_accessions
    partition_accessions
    download_all
    write_genome_log
    write_availability

    echo ""
    print_info "Done"
}

main
