#!/usr/bin/env bash
set -euo pipefail

# Load required parameters for Efficient DV workflow
config_file="$1"

if [[ -z "$config_file" ]]; then
    echo "ERROR: You must provide a config file." >&2
    exit 1
fi

if [[ ! -f "$config_file" ]]; then
    echo "ERROR: Config file not found: $config_file" >&2
    exit 1
fi

# Load all variables
source "$config_file"

# Import functions
source "${WGS_mainpath}/utils/MakeExamples.sh"

# Change working directory to the one provided in 
cd $WGS_mainpath

# Enable nullglob to safely handle glob patterns with no matches
shopt -s nullglob

# 01. SCATTER INTERVAL LIST

# Check intervals files
bed_path="${WGS_mainpath}/00_ScatteredIntervalList"
mkdir -p "$bed_path"

if [ -s "${bed_path}/intervals.txt" ]; then
    echo "[$(date)] Found existing ${bed_path}/intervals.txt"
else
    echo ""
    
    echo "[$(date)] Generating ${bed_path}/intervals.txt"
    find "$bed_path" -maxdepth 1 -type f -name '*.bed' | sort -V > "${bed_path}/intervals.txt"

fi

# Assign number of threads based on number of interval files
threads=$(find "$bed_path" -maxdepth 1 -name "*.bed" | wc -l)
if [ "$threads" -eq 0 ]; then
    echo "ERROR: No BED files found in $bed_path" >&2
    exit 1
fi

echo "[$(date)] Found $threads BED intervals for processing"

# 02. MAKE EXAMPLES

MakeExamplesParallel() {
    
    local bed_file="$1"
    local in_cram="$2"
    local ref_fasta="$3"
    local output_sample_path="$4"
    local wgs_mainpath="$5"

    bed_name=$(basename "$bed_file" .bed)

    output_sample_log_path="${output_sample_path}/log"
    mkdir -p "$output_sample_log_path"
    output_sample_bed_log="${output_sample_log_path}/${bed_name}_MakeExamples.log"

    expected_output="${output_sample_path}/${bed_name}.tfrecord.gz"
    if [ -s "$expected_output" ]; then
        echo "[$(date)] Found existing MakeExamples output: $expected_output"
        return 0
    fi

    MakeExamples \
        "${wgs_mainpath}" \
        "${in_cram}" \
        "${ref_fasta}" \
        "${bed_file}" \
        "${output_sample_path}" \
        > "${output_sample_bed_log}" 2>&1 || {
        echo "ERROR: MakeExamples failed for $bed_name" >&2
        return 1
    }

}

export -f MakeExamplesParallel
export WGS_mainpath
export in_cram
export ref_fasta
export output_sample_path

# Loop through all cram files
work_cram_dir="01_cram"

# Make examples folder
make_examples_path="02_MakeExamples"
mkdir -p "$make_examples_path"

# Call variants folder
call_variants_path="03_CallVariants"
mkdir -p "$call_variants_path"

# Variant postprocessing folder
postproc_path="04_PostProcessing"
mkdir -p "$postproc_path"

# Variant annotation folder
VariantAnnot_path="05_VariantAnnot"
mkdir -p "$VariantAnnot_path"

VariantAnnot_log_path="${VariantAnnot_path}/log"
mkdir -p "$VariantAnnot_log_path"

input_cram_list=($(find ${input_cram_dir} -name '*.cram' | sort -V))

for cram in "${input_cram_list[@]}"; do

    samplename=$(basename "$cram" .cram)

    # Skip if final VCFs already exist
    vcfannot_prefix="${VariantAnnot_path}/${samplename}.pass_annovar"
    vcfannot_out="${vcfannot_prefix}.hg38_multianno.vcf"
    
    if [[ -s "$vcfannot_out" ]]; then
        echo "[$(date)] Skipping ${samplename}: final annotated VCF already exists."
        continue
    else
        echo "[$(date)] Starting variant calling for ${samplename}..."
    fi

    output_sample_path="${make_examples_path}/${samplename}"
    mkdir -p "$output_sample_path"

    # Run in parallel
    echo "[$(date)] Starting MakeExamples for ${in_cram} with $threads parallel jobs..."
    haplotype_failed=0
    find "$bed_path" -maxdepth 1 -name "*.bed" | sort -V | \
        xargs -P "${threads}" -I {} bash -c '
            MakeExamplesParallel "$1" "$2" "$3" "$4" "$5"
            ' _ "{}" "$in_cram" "$ref_fasta" "$output_sample_path" "$WGS_mainpath" _ || haplotype_failed=1

    if [ $haplotype_failed -eq 0 ]; then
        echo "[$(date)] MakeExamples for ${samplename} completed successfully"
        
        # Verify output files were created
        output_count=$(find "$output_sample_path" -name "*.tfrecord.gz" | wc -l)
        if [ "$output_count" -gt 0 ]; then
            echo "[$(date)] Found $output_count tfrecord files. Proceeding to CallVariants..."
        else
            echo "[$(date)] ERROR: No tfrecord output files found in ${output_sample_path}" >&2
            exit 1
        fi
    else
        echo "[$(date)] ERROR: MakeExamples for ${samplename} failed" >&2
        exit 1
    fi

    # 03. CALL VARIANTS

    # Take the params.ini file and add the tfrecord.gz paths
    if [ ! -f "$params_ini_template" ]; then
        echo "ERROR: Template params.ini not found: $params_ini_template" >&2
        exit 1
    fi

    params_ini_path="$call_variants_path/${samplename}_params.ini"
    mkdir -p "$(dirname "$params_ini_path")"
    cp "$params_ini_template" "$params_ini_path"
    chmod 644 "$params_ini_path"

    outputFileName="$call_variants_path/${samplename}"
    mkdir -p $outputFileName

    # Replace existing outputFileName line or append if missing
    if grep -q '^outputFileName' "${params_ini_path}"; then
        tmpfile=$(mktemp)
        sed "s|^outputFileName.*|outputFileName = $samplename|" "${params_ini_path}" > "$tmpfile"
        mv "$tmpfile" "${params_ini_path}"
    else
        echo "outputFileName = $samplename" >> "${params_ini_path}"
    fi

    # Count actual tfrecord files produced by MakeExamples
    actual_tfrecord_count=$(find "$output_sample_path" -maxdepth 1 -name "*scattered.tfrecord.gz" | wc -l)
    echo "[$(date)] Found $actual_tfrecord_count tfrecord files from MakeExamples"

    if [ "$actual_tfrecord_count" -eq 0 ]; then
        echo "[$(date)] ERROR: No tfrecord files found in ${output_sample_path}" >&2
        exit 1
    fi

    # Add tfrecord file paths to params.ini
    for i in $(seq 1 "$actual_tfrecord_count"); do
        tfrecord_path="${make_examples_path}/${samplename}/${i}scattered.tfrecord.gz"
        echo "exampleFile${i} = ${tfrecord_path}" >> "$params_ini_path"
    done

    output_sample_callvariants_log="$call_variants_path/${samplename}_CallVariants.log"

    # Check if output files already exist
    echo "[$(date)] actual_tfrecord_count=$actual_tfrecord_count"
    existing_outputs=0

    for i in $(seq 1 "$actual_tfrecord_count"); do
        expected_output="$call_variants_path/${samplename}/${samplename}.${i}.gz"
        
        if [ -s "$expected_output" ]; then
            echo "[$(date)] Checking for existing output: $expected_output"
            existing_outputs=$((existing_outputs+1))
        fi
    done

    echo "[$(date)] existing_outputs count: $existing_outputs out of $actual_tfrecord_count"

    if [ $existing_outputs -eq "$actual_tfrecord_count" ]; then
        echo "[$(date)] All CallVariants output files for ${samplename} already exist - skipping CallVariants step"
    else
        echo "[$(date)] Running CallVariants for ${samplename}..."
        mkdir -p "$call_variants_path/${samplename}"

        if docker run --rm --gpus all \
            --volume $WGS_mainpath:$WGS_mainpath \
            --workdir $WGS_mainpath \
            ultimagenomics/call_variants:2.2.4 call_variants \
            --parameters-file "${params_ini_path}" \
            > "${output_sample_callvariants_log}" 2>&1; then
            echo "[$(date)] CallVariants completed successfully"
            
            # Move all generated files from call_variants to a specific folder
            mv ${samplename}.*.gz $outputFileName
        else
            echo "[$(date)] ERROR: CallVariants failed (check log: ${output_sample_callvariants_log})" >&2
            exit 1
        fi
    fi

    # 04. POST PROCESSING

    echo "[$(date)] Starting post-processing for ${samplename}..."

    # Collect in comma-separated list the call_variants output files
    call_variants_outputs=""
    for i in $(seq 1 "$actual_tfrecord_count"); do
        call_variants_output="$call_variants_path/${samplename}/${samplename}.${i}.gz"
        if [ -s "$call_variants_output" ]; then
            if [ -z "$call_variants_outputs" ]; then
                call_variants_outputs="${call_variants_output}"
            else
                call_variants_outputs="${call_variants_outputs},${call_variants_output}"
            fi
        else
            echo "[$(date)] ERROR: Missing CallVariants output file: ${call_variants_output}" >&2
            exit 1
        fi
    done

    # Collect in comma-separated list the call_variants .g.vcf.gz files
    call_gvariants_outputs=""
    for i in $(seq 1 "$actual_tfrecord_count"); do
        call_gvariants_output="${make_examples_path}/${samplename}/${i}scattered.gvcf.tfrecord.gz"
        if [ -s "$call_gvariants_output" ]; then
            if [ -z "$call_gvariants_outputs" ]; then
                call_gvariants_outputs="${call_gvariants_output}"
            else
                call_gvariants_outputs="${call_gvariants_outputs},${call_gvariants_output}"
            fi
        else
            echo "[$(date)] WARNING: Missing GVCF output file (optional): ${call_gvariants_output}" >&2
        fi
    done

    postproc_sample_path="${postproc_path}/${samplename}"
    mkdir -p "$postproc_sample_path"

    if [ -s "${postproc_sample_path}/${samplename}.vcf.gz" ] && [ -s "${postproc_sample_path}/${samplename}.g.vcf.gz" ]; then
        echo "[$(date)] Post-processing outputs for ${samplename} already exist - skipping Post-processing step"
    else
        echo "[$(date)] Running Post-processing for ${samplename}..."
        
        # Build docker command with all required parameters
        if docker run --rm --gpus all \
            --volume $WGS_mainpath:$WGS_mainpath \
            --workdir $WGS_mainpath \
            ultimagenomics/make_examples:3.1.9 ug_postproc \
            --infile "${call_variants_outputs}" \
            --ref "${ref_fasta}" \
            --outfile "${postproc_sample_path}"/"${samplename}".vcf.gz \
            --gvcf_outfile "${postproc_sample_path}"/"${samplename}".g.vcf.gz \
            --nonvariant_site_tfrecord_path "$call_gvariants_outputs" \
            --consider_strand_bias \
            --flow_order TGCA \
            --annotate \
            --qual_filter 1 \
            --filter \
            --filters_file filter_params.txt \
            --bed_annotation_files "$exome_twist","$LCR_hs38","$mappability","$hmers_7_and_higher","$ug_hcr" \
            --dbsnp "$dbsnp" \
            > "${postproc_sample_path}/${samplename}_PostProc.log" 2>&1; then
            echo "[$(date)] Post-processing completed successfully"

            bcftools stats "${postproc_sample_path}/${samplename}.vcf.gz" > "${postproc_sample_path}/${samplename}_stats.txt"

            echo "[$(date)] Filtering PASS variants..."
            bcftools view -f PASS -O z "${postproc_sample_path}/${samplename}.vcf.gz" -o "${postproc_sample_path}/${samplename}.pass.vcf.gz"
            bcftools index -t "${postproc_sample_path}/${samplename}.pass.vcf.gz"
            bcftools stats "${postproc_sample_path}/${samplename}.pass.vcf.gz" > "${postproc_sample_path}/${samplename}_stats_pass.txt"

        else
            echo "[$(date)] ERROR: Post-processing failed (check log: ${postproc_sample_path}/${samplename}_PostProc.log)" >&2
            exit 1
        fi
    fi

    # 05. VARIANT ANNOTATION

    # Load annovar path and human database folder
    mkdir -p "$humandb"

    echo "[$(date)] Checking required ANNOVAR databases..."

    databases=(
        refGeneWithVer
        gnomad41_genome
        avsnp151
        clinvar_20250721
        GTEx_v8_eQTL
    )

    for db in "${databases[@]}"; do
        if ! compgen -G "${humandb}/hg38_${db}.*" >/dev/null; then
            echo "[$(date)] Database '$db' not found. Downloading..."
            perl "${annovar_path}/annotate_variation.pl" -buildver hg38 -downdb -webfrom annovar "$db" "$humandb"
        else
            echo "[$(date)] Database '$db' already present."
        fi
    done

    echo "[$(date)] All required databases verified."

    filt_vcf="${postproc_sample_path}/${samplename}.pass.vcf.gz"
    vcfname="$(basename "$filt_vcf" .vcf.gz)"

    # Output prefix
    annovar_log="${VariantAnnot_log_path}/${vcfname}_ANNOVAR.log"

    # Skip if final VCFs already exist
    if [[ -s "$vcfannot_out" ]]; then
        echo "[$(date)] Skipping ${samplename}: final file already exists."
        continue
    fi

    echo "[$(date)] Annotating variants for ${samplename} with ANNOVAR..."

    # Run ANNOVAR annotation
    perl "$annovar_path/table_annovar.pl" "$filt_vcf" "$humandb/" \
        --buildver hg38 \
        --out "$vcfannot_prefix" \
        --remove \
        --protocol refGeneWithVer,gnomad41_genome,avsnp151,clinvar_20250721 \
        --operation g,f,f,f \
        --arg '-hgvs',,, \
        --nastring . \
        --vcfinput \
        --polish \
        --intronhgvs 20 \
        --thread 64 \
        > "${annovar_log}" 2>&1

    echo "[$(date)] Completed annotation for ${samplename}"

done


    
