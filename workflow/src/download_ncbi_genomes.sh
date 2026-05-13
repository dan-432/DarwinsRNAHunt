#!/bin/bash

################################################################################
# Script modified from Claude.ai generated code
# NCBI Genome Download Script
# Downloads genome sequences and GFF3 annotations for a list of accessions
################################################################################

set -euo pipefail

# Configuration
GENOME_ACCESSIONS_FILE="${1}"
OUTPUT_DIR="${2}"
INCLUDE_FORMATS="${3:-gff3,genome,protein}"
BATCH_SIZE="${4:-50}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
print_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

check_dependencies() {
    if ! command -v datasets &> /dev/null; then
        print_error "NCBI datasets CLI not found! Install: conda install -c conda-forge ncbi-datasets-cli"
        exit 1
    fi
    if ! command -v jq &> /dev/null; then
        print_error "jq not found! Install: conda install -c conda-forge jq"
        exit 1
    fi
}

extract_accessions() {
    local input_file=$1
    local output_file=$2
    
    if jq empty "$input_file" 2>/dev/null; then
        jq -r '.[]' "$input_file" > "$output_file"
    else
        grep -v '^#' "$input_file" | grep -v '^$' | sed 's/[[:space:]]*$//' > "$output_file"
    fi
    
    local num=$(wc -l < "$output_file")
    if [[ $num -eq 0 ]]; then
        print_error "No accessions found"
        exit 1
    fi
    print_info "Found $num accessions"
}

download_batch() {
    local batch_file=$1
    local batch_num=$2
    local temp_zip="$OUTPUT_DIR/temp_batch_${batch_num}.zip"
    local temp_dir="$OUTPUT_DIR/temp_batch_${batch_num}"
    local max_retries=3
    local retry_delay=10
    
    # Retry logic for network issues
    local attempt=1
    while [[ $attempt -le $max_retries ]]; do
        print_info "Downloading batch $batch_num (attempt $attempt/$max_retries)..."
        
        if datasets download genome accession \
            --inputfile "$batch_file" \
            --include "$INCLUDE_FORMATS" \
            --filename "$temp_zip" \
            --no-progressbar 2>&1; then
            break
        else
            if [[ $attempt -lt $max_retries ]]; then
                print_warn "Download failed, retrying in ${retry_delay}s..."
                sleep $retry_delay
                attempt=$((attempt + 1))
            else
                print_error "Failed to download batch $batch_num after $max_retries attempts"
                return 1
            fi
        fi
    done
    
    # Extract to temp directory
    unzip -q -o "$temp_zip" -d "$temp_dir"
    rm "$temp_zip"
    
    # Reorganize files and track successful downloads (only if all formats present)
    for accession_dir in "$temp_dir/ncbi_dataset/data"/*; do
        if [[ -d "$accession_dir" ]]; then
            accession=$(basename "$accession_dir")
            target_dir="$OUTPUT_DIR/$accession"
            mkdir -p "$target_dir"
            
            # Copy and rename files
            find "$accession_dir" -name "*_genomic.fna" -exec cp {} "$target_dir/genomic.fna" \;
            find "$accession_dir" -name "*.gff" -exec cp {} "$target_dir/geneannotation.gff" \;
            find "$accession_dir" -name "*protein.faa" -exec cp {} "$target_dir/protein.faa" \;
            
            # Verify all formats are present before marking as successful
            if [[ -f "$target_dir/genomic.fna" && -f "$target_dir/geneannotation.gff" && -f "$target_dir/protein.faa" ]]; then
                echo "$accession" >> "$OUTPUT_DIR/.downloaded_accessions.txt"
            else
                # Remove incomplete genome
                print_warn "Incomplete download for $accession (missing formats), removing directory"
                rm -rf "$target_dir"
            fi
        fi
    done
    
    # Clean up temp directory
    rm -rf "$temp_dir"
    return 0
}

create_metadata() {
    local metadata_file="$OUTPUT_DIR/download_metadata.json"
    local download_date
    download_date=$(date -u +"%d-%m-%YT%H:%M:%SZ")
    
    # Get list of successfully downloaded accessions (compatible with macOS)
    local downloaded_accessions=""
    if [[ -f "$OUTPUT_DIR/.downloaded_accessions.txt" ]]; then
        downloaded_accessions=$(sort -u "$OUTPUT_DIR/.downloaded_accessions.txt" | jq -R . | jq -s .)
    else
        downloaded_accessions="[]"
    fi
    
    # Create JSON metadata
    echo "{
  \"download_date\": \"$download_date\",
  \"source\": \"NCBI\",
  \"formats\": \"$INCLUDE_FORMATS\",
  \"total_accessions\": $(echo "$downloaded_accessions" | jq 'length'),
  \"accessions\": $downloaded_accessions
}" > "$metadata_file"
    
    # Clean up tracking file
    rm -f "$OUTPUT_DIR/.downloaded_accessions.txt"
    
    print_info "Created metadata file: $metadata_file"
}

main() {
    print_info "Starting NCBI genome download"
    
    check_dependencies
    
    if [[ ! -f "$GENOME_ACCESSIONS_FILE" ]]; then
        print_error "File not found: $GENOME_ACCESSIONS_FILE"
        exit 1
    fi
    
    mkdir -p "$OUTPUT_DIR"
    
    # Extract accessions
    local clean_accessions="$OUTPUT_DIR/accessions_clean.txt"
    extract_accessions "$GENOME_ACCESSIONS_FILE" "$clean_accessions"
    
    # Split into batches
    local total=$(wc -l < "$clean_accessions")
    local num_batches=$(( (total + BATCH_SIZE - 1) / BATCH_SIZE ))
    
    print_info "Splitting into $num_batches batches"
    split -l "$BATCH_SIZE" -d -a 3 "$clean_accessions" "$OUTPUT_DIR/batch_"
    
    # Download batches
    local batch_num=0
    local failed=0
    
    for batch_file in "$OUTPUT_DIR"/batch_*; do
        if [[ -f "$batch_file" ]] && [[ ! "$batch_file" =~ \.failed$ ]]; then
            batch_num=$((batch_num + 1))
            
            if ! download_batch "$batch_file" "$batch_num"; then
                failed=$((failed + 1))
                mv "$batch_file" "${batch_file}.failed"
            else
                rm "$batch_file"
            fi
            
            # Longer delay between batches to avoid rate limiting
            sleep 2
        fi
    done
    
    # Create metadata
    create_metadata
    
    # Cleanup
    rm -f "$clean_accessions"
    
    echo ""
    print_info "================================================"
    print_info "Download complete!"
    print_info "Total batches: $num_batches"
    print_info "Failed batches: $failed"
    print_info "Metadata: $OUTPUT_DIR/download_metadata.json"
    print_info "================================================"
    
    if [[ $failed -gt 0 ]]; then
        print_warn "Retry failed batches with:"
        for failed_batch in "$OUTPUT_DIR"/batch_*.failed; do
            if [[ -f "$failed_batch" ]]; then
                echo "  $0 $failed_batch $OUTPUT_DIR"
            fi
        done
    fi
}

main