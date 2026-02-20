#!/usr/bin/env bash
set -euo pipefail

MakeExamples() {
    local mainpath="$1"
    local in_cram_list="$2"
    local ref_file="$3"
    local bed_file="$4"
    local output_path="$5"
    
    # Validate required arguments
    if [[ -z "${mainpath}" ]] || [[ -z "${in_cram_list}" ]] || \
    [[ -z "${ref_file}" ]] || [[ -z "${bed_file}" ]] || [[ -z "${output_path}" ]]; then
        echo "ERROR: Missing required arguments for MakeExamples" >&2
        echo "Usage: MakeExamples <mainpath> <input_cram_list> <input_cram_index_list> <ref_file> <bed_file> <output_path>" >&2
        return 1
    fi

    cram_file=${in_cram_list#$WORKDIR/}
    cram_idx="${cram_file}.crai"
    echo "$cram_file"

    ref_file=${ref_file#$WORKDIR/}
    bed_file=${bed_file#$WORKDIR/}
    output_path=${output_path#$WORKDIR/}

    if [[ ! -f "${mainpath}/${bed_file}" ]]; then
        echo "ERROR: BED file not found: ${mainpath}/${bed_file}" >&2
        return 1
    fi
    if [[ ! -f "${mainpath}/${ref_file}" ]]; then
        echo "ERROR: Reference file not found: ${mainpath}/${ref_file}" >&2
        return 1
    fi
    
    out_number=$(basename "$bed_file" .bed)

    echo "[$(date)] Running MakeExamples:"
    echo "[$(date)] ##################################################"
    echo "[$(date)] Working path: ${mainpath}..."
    echo "[$(date)] Bed file: ${bed_file}..."
    echo "[$(date)] CRAM: ${cram_file}..."
    echo "[$(date)] CRAM index: ${cram_idx}..."
    echo "[$(date)] Reference ${ref_file}..."
    echo "[$(date)] Output path and name ${output_path}/${out_number}..."
    echo "[$(date)] ##################################################"

    # Run docker with relative paths (mainpath is the working directory)
    docker run --rm \
        --volume $mainpath:$mainpath \
        --workdir $mainpath \
        ultimagenomics/make_examples:3.1.9 tool \
        --input $cram_file \
        --cram-index $cram_idx \
        --bed "${bed_file}" \
        --output "${output_path}/${out_number}" \
        --reference "${ref_file}" \
        --min-base-quality 5 \
        --min-mapq 5 \
        --gvcf \
        --p-error 0.005 \
        --cgp-min-count-snps 2 \
        --cgp-min-count-hmer-indels 2 \
        --cgp-min-count-non-hmer-indels 2 \
        --cgp-min-fraction-snps 0.12 \
        --cgp-min-fraction-hmer-indels 0.12 \
        --cgp-min-fraction-non-hmer-indels 0.06 \
        --cgp-min-mapping-quality 5 \
        --max-reads-per-region 1500 \
        --assembly-min-base-quality 0 \
        --optimal-coverages 50 \
        --add-ins-size-channel || {
        echo "ERROR: MakeExamples docker failed with exit code $?" >&2
        return 1
    }
    
    # Verify output was created
    if [[ ! -s "$expected_output" ]]; then
        echo "ERROR: Expected output not found after MakeExamples: $expected_output" >&2
        return 1
    fi
    
    echo "[$(date)] MakeExamples completed successfully for ${samplename} region ${out_number}"
    return 0
}

export -f MakeExamples
