#!/bin/bash

set -euo pipefail

# load your settings:
source "$(dirname "$0")/variant_call_config.sh"

# ──────────────────────────────────────────────────────────────────────────────
# Pipeline helper functions: download, index, trim, align, & variant‐call.
# ──────────────────────────────────────────────────────────────────────────────


timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
logmsg() { echo "$(timestamp): $*"; }


# create directory function 
create_dir(){
    for dir in "$@"; do
        if [[ -d ${dir} ]] ; then
            logmsg "directory $dir already exists."
        else 
            if mkdir -p ${dir}; then
                logmsg  "directory ${dir} created"
            else
                logmsg  "directory ${dir} not created" ; exit 1
            fi
        fi
    done
}

###############################################################################
# downlad_genome function
#   • Looks up URL in genome_urls[] using the tag (e.g. b73, mo17…)
#   • Downloads the compressed FASTA to $ref_dir if missing
#   • Logs every step
###############################################################################
download_genome() {
  local gbase="$1"
  local url="${genome_urls[$gbase]:-}"               
  local genome="${ref_dir}/${gbase}.fa.gz"
  
  # make parent dirs exist (creates them if missing)
  create_dir "$ref_dir"

  logmsg "Checking genome $gbase : $genome"


    if [[ ${do_download}=="true" ]] ; then 
    logmsg "downloading $gbase required"
        if [[ -s ${genome} ]]; then
            logmsg "Genome ${genome} already exists — skipping download."
        else
            logmsg "Downloading ${gbase} genome..."
            if wget -O ${genome} ${url}; then
            logmsg "Successfully downloaded ${gbase}."
            else
            logmsg "ERROR: download failed or URL missing for ${genox}."
            fi
        fi
    else
        logmsg "download $gbase not required ... skipping"
    fi
}


# helper function to index genomes

index_genome () {
    local gbase=$1
    local genome="${ref_dir}/${gbase}.fa.gz"
    local idx_dir="${ref_dir}/indexs"
    local mmi="${idx_dir}/${gbase}.mmi"
    local link_fa="${idx_dir}/${gbase}.fa.gz"

    create_dir "$idx_dir"

    # check genome
    if [[ ! -s "$genome" ]]; then
        logmsg "ERROR: $genome not found or empty — cannot index."
        return
    fi

    # index the genome using minimap2
    logmsg "Starting Minimap2 indexing for ${genome}"
    if [[ -s "${mmi}" ]]; then
        logmsg "Minimap2 index exists for $gbase — skipping."
    else
        logmsg "Building minimap2 index for $gbase …"
        if minimap2 -t ${threads} -d "${mmi}" ${genome} ; then
            logmsg "Minimap2 index complete: $mmi"
        else
            logmsg "Minimap2 index failed for ${genome}." &&  exit 1 
        fi
    fi

    # index the genome using samtools fadix .fai & .gzi
    if [[ -s "${idx_dir}/${gbase}.fai" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S'): samtools index for ${gbase}.fai already exists. Skipping indexing."
    else
        logmsg "Running samtools faidx on $genome…"
        if samtools faidx ${genome} ; then
            logmsg "Samtools faidx complete for ${gbase}"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S'): samtools index failed for ${gbase}." &&  exit 1
        fi
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S'): ${gbase} genome samtool indexing completed"

    
     # Move auxiliary files to index directory
    for ext in gzi fai; do
        if [[ -f "${genome}.${ext}" ]]; then
            mv -u "${genome}.${ext}" "${idx_dir}" && \
            echo "$(date '+%Y-%m-%d %H:%M:%S'): Moved ${genome}.${ext} to ${idx_dir}" 
        fi
    done

    # create FASTA symlink beside the indexes (handy for downstream)
    if [[ -e "$link_fa" ]]; then
        logmsg "Symlink already exists: $link_fa"
    else
        if ln -s "$genome" "$link_fa" ; then
            logmsg "Symlink created: $link_fa  for $genome"
        else
            logmsg "Error : failed to create symbolic link for genome: $gbase"  &&  exit 1
        fi
    fi
}

trim_reads () {
    local sample=$1 
    local sample_dir="${sample_loc[$sample]}"
    local reads_dir="${sample_dir}/reads"
    local trim_dir="${sample_dir}/minimap_dir/trim_dir"
    local stat_dir="${sample_dir}/minimap_dir/stat_dir"

    create_dir ${trim_dir} ${stat_dir}


    logmsg "${sample_dir}"

    # Initialize the summary file
    local summary="$stat_dir/reads_trim_summary.tsv"
    if [[ -f ${summary} ]]; then 
        echo "$(date '+%Y-%m-%d %H:%M:%S'): ${summary} already exist."
    else
        # Add a header to the summary file
        echo -e "Sample\tOriginal_Reads_R1\tOriginal_Reads_R2\tTrimmed_Reads_R1\tTrimmed_Reads_R2\t\
        Fastp_Version\tTotal_Reads_Before\tTotal_Bases_Before\tQ20_Bases_Before\tQ30_Bases_Before\tGC_Content_Before\t\
        Total_Reads_After\tTotal_Bases_After\tQ20_Bases_After\tQ30_Bases_After\tGC_Content_After\tPassed_Filter_Reads\t\
        Low_Quality_Reads\tToo_Many_N_Reads\tToo_Short_Reads\tDuplication_Rate\tInsert_Size_Peak" > "$summary" && \
        logmsg "summary file created at: $summary"
    fi
    
    # Trim reads
    for R1 in ${reads_dir}/*_R1.fq.gz; do
        local sbase=$(basename ${R1} _R1.fq.gz)
        local R2="${reads_dir}/${sbase}_R2.fq.gz"
        local O1=${trim_dir}/${sbase}_trim_R1.fq.gz
        local O2=${trim_dir}/${sbase}_trim_R2.fq.gz
        local ht=${trim_dir}/${sbase}_fastp_report.html
        local jt=${trim_dir}/${sbase}_fastp_report.json


        # Ensure trim directory exists
        if [[ ! -s ${O1} ]] && [[ ! -s ${O2} ]]; then
            if fastp -i $R1 -I $R2 -o $O1 -O $O2 --detect_adapter_for_pe -h ${ht} -j ${jt}; then
                # Calculate total reads in the output files (divide line count by 4)
                local trR1=$(zcat $R1 | wc -l | awk '{print $1/4}')
                local trR2=$(zcat $R2 | wc -l | awk '{print $1/4}')
                local trO1=$(zcat $O1 | wc -l | awk '{print $1/4}')
                local trO2=$(zcat $O2 | wc -l | awk '{print $1/4}')
            fi

            # Extract metrics from the JSON file and store as a string
            local metrics=$(jq -r '
                [.summary.fastp_version, .summary.before_filtering.total_reads, .summary.before_filtering.total_bases, 
                .summary.before_filtering.q20_bases, .summary.before_filtering.q30_bases, .summary.before_filtering.gc_content, 
                .summary.after_filtering.total_reads, .summary.after_filtering.total_bases, 
                .summary.after_filtering.q20_bases, .summary.after_filtering.q30_bases, .summary.after_filtering.gc_content, 
                .filtering_result.passed_filter_reads, .filtering_result.low_quality_reads, .filtering_result.too_many_N_reads, 
                .filtering_result.too_short_reads, .duplication.rate, .insert_size.peak ] | @tsv' "${jt}")

            # Append detailed metrics to the summary file
            echo -e "${sbase}\t${trR1}\t${trR2}\t${trO1}\t${trO2}\t${metrics}" >> "$summary"

            logmsg "fastp complete and stats added for $sbase"
        else
            logmsg "$(basename ${O1}) and $(basename ${O1}) already exists and trimmed. skipping..."
        fi
    done
}



map_reads (){
    local sample=$1 
    local gbase=$2
    local sample_dir="${sample_loc[$sample]}"
    local trim_dir="${sample_dir}/minimap_dir/trim_dir"
    local stat_dir="${sample_dir}/minimap_dir/stat_dir"
    local align_dir="${sample_dir}/minimap_dir/align_dir"
    local plots_dir="${sample_dir}/minimap_dir/plots_dir"
    local vcf_dir="${sample_dir}/minimap_dir/vcf_dir"
    local snps_dir="${sample_dir}/minimap_dir/snps_dir"
    local idx_dir="${ref_dir}/indexs"

    create_dir ${align_dir} ${plots_dir} ${vcf_dir} ${stat_dir} ${snps_dir}

    # reference paths
    local idx_mmi="$idx_dir/${gbase}.mmi"
    local genome_fa="$idx_dir/${gbase}.fa.gz"

    if [[ -s ${idx_mmi} && -s ${genome_fa} ]]; then
        logmsg "Using index: $idx_mmi and FASTA: $genome_fa"
    else
        logmsg "ERROR – index or FASTA missing for $gbase" &&  exit 1
    fi

    #iterate over trimmed pair
    for O1 in ${trim_dir}/*_trim_R1.fq.gz; do
        local sbase=$(basename "$O1" _trim_R1.fq.gz)
        local O2="$trim_dir/${sbase}_trim_R2.fq.gz"
        local tag="${gbase}_${sbase}"
        local sam="${align_dir}/${tag}.sam"
        local bam="${align_dir}/${tag}.bam"
        local vcf="${vcf_dir}/${tag}.vcf.gz"

        # Run Minimap2 alignment
        logmsg "mapping ${sample} to ${gbase} started"
        if [[ -s "$bam" ]]; then
            logmsg  "${bam} already exists. Skipping alignment."
        else
            if minimap2 -ax sr -t ${threads} "${idx_mmi}" ${O1} ${O2} > ${sam}; then
                logmsg "Alignment completed for ${sbase} & $(basename ${sam}) created"
            else
                logmsg "Aligment failed ${sam}. exiting"  && exit 1
            fi

            # Convert SAM to BAM, sort, and index
            logmsg "samtools viewing and sorting ${bam} started"
            if samtools view -bS ${sam} | samtools sort -o ${bam}; then
                logmsg  "samtools viewing and sorting completeed for ${sbase} & $(basename ${bam}) created ."
            else
                logmsg "creation of ${bam} failed. exiting"  &&  exit 1
            fi

            # Remove the SAM file after successful processing
            logmsg "Removing ${sam}" echo 
            if rm -f ${sam}; then
                logmsg "Removed ${sam}"
            else
                logmsg "Failed to remove ${sam}"
            fi

            # Index the BAM file
            logmsg "samtools indexing ${bam} started"
            if samtools index ${bam}; then
                logmsg "${bam} sucessfully indexed"
            else
                logmsg "${bam} failed. exiting"  &&  exit 1
            fi
        fi


        # Check if overall coverage file exists, create if missing
        local overall_coverage="${stat_dir}/overall_coverage.tsv"

        # create the header only if missing
        # Initialize the merged coverage file (header only)
        if [[ ! -s "${overall_coverage}" ]]; then
            echo -e "Genome\tSample\tCoverage\tTotal_Reads\tMapped_Reads\tPaired_Reads\tAlignment_%" > "${overall_coverage}" \
                && logmsg "Successfully created ${overall_coverage}" || { logmsg "Error creating ${overall_coverage}"; exit 1; }
        fi


        # Run samtools flagstat and extract relevant values
        local stats total_reads mapped_reads properly_paired depth
        logmsg "Flagstat for $bam"
        stats=$(samtools flagstat "${bam}")

        # Extract values using awk
        total_reads=$(echo "$stats" | awk 'NR==1 {print $1}') || total_reads=0
        mapped_reads=$(echo "$stats" | awk 'NR==7 {print $1}') || mapped_reads=0
        properly_paired=$(echo "$stats" | awk 'NR==12 {print $1}') || properly_paired=0

        # Ensure extracted values are valid
        if [[ -z "$total_reads" || -z "$mapped_reads" || "$total_reads" -eq 0 ]]; then
            echo -e "Error:\tCould not retrieve valid read statistics from BAM file."
            exit 1
        fi

        # Compute percentages safely (avoid division by zero)
        local alignment_percentage properly_paired_percentage
        if [[ "${total_reads}" -gt 0 ]]; then
            alignment_percentage=$(awk -v total="${total_reads}" -v mapped="${mapped_reads}" 'BEGIN {printf "%.2f", (mapped/total)*100}')
            properly_paired_percentage=$(awk -v total="${total_reads}" -v paired="${properly_paired}" 'BEGIN {printf "%.2f", (paired/total)*100}')
        else
            alignment_percentage=0
            properly_paired_percentage=0
        fi

        # compute coverage depth using samtoools - depth
        logmsg "computing coverage depth ${bam} started"
        if depth=$(samtools depth -a ${bam} | awk '{sum+=$3} END { print sum/NR }'); then 
            logmsg "Coverage calculated for ${bam}: ${depth}x"
        else
            logmsg "Coverage calculation failed for ${bam}" &&  exit 1
        fi

        # Append structured tab-separated results
        echo -e "${gbase}\t${sbase}\t${depth}\t${total_reads}\t${mapped_reads}\t${properly_paired}\t${alignment_percentage%}" >> "${overall_coverage}"

        # Variant calling with bcftools mpileup + call
        logmsg "checking output vcf ${vcf}"
        logmsg  "mpileup and variant calling ${bam} against ${gbase} started"
        if [[ -f ${vcf} && -s ${vcf} ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S'): ${vcf} already exists. Skipping."
        else
            if ! bcftools mpileup --ignore-RG -f "${genome_fa}" "${bam}" --threads "${threads}" \
                    | bcftools call -m -v -Oz --threads "${threads}" -o "${vcf}"; then
                echo "$(date '+%Y-%m-%d %H:%M:%S'): bcftools mpileup for ${tag} failed" && exit 1
            fi
            echo "$(date '+%Y-%m-%d %H:%M:%S'): mpileup and variant calling for ${tag} completed"

            # index vcf files
            echo "$(date '+%Y-%m-%d %H:%M:%S'): indexing ${vcf} started"
            if ! bcftools index "${vcf}"; then
                echo "$(date '+%Y-%m-%d %H:%M:%S'): bcftools index for ${vcf} failed" && exit 1
            fi
            echo "$(date '+%Y-%m-%d %H:%M:%S'): indexing ${vcf} completed"
        fi

        # bcftool stats
        local bstats="${stat_dir}/${tag}_bcfstats.tsv"
        logmsg "Generating stats for VCF: ${vcf}"
        if [[ -s ${bstats} ]]; then
        logmsg "Stats file already exists: ${bstats}"
        else
            if ! bcftools stats "${vcf}" > "${bstats}"; then
                logmsg "ERROR: bcftools stats failed on ${vcf}" && exit 1
            else 
                logmsg "bcftools stats written to ${bstats}"
            fi     
        fi

        # Statistics and Plotting
        local bplots="${plots_dir}/${tag}_plots"
        logmsg "Plotting stats from ${bstats}"
        if [[ -s ${bplots} ]]; then
            logmsg "Plots directory already exists: ${bplots}"
        else
            if ! plot-vcfstats -t "${tag}" -p "${bplots}" "${bstats}"; then
                logmsg "WARNING: plot-vcfstats failed for ${bstats}"
            else
                logmsg "Plots generated in ${bplots}"
            fi 
        fi

        # decompress .vcf.gz and pipe to bcftools
        local snp_table="${snps_dir}/${tag}_snps.tsv"
        
        # Create readable snp files for downstream analysis
        logmsg "Creating the Final SNP table for ${snp_table} started"

        if [[ -f ${snp_table} && -s ${snp_table} ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S'): ${snp_table} already exists. Skipping."
        else
            echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tFref\tRref\tFalt\tRalt" > ${snp_table}

            if bgzip -d -c ${vcf} | grep -E '^#|^chr[1-9]\b|^chr10\b' | \
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t[%DP4{0}]\t[%DP4{1}]\t[%DP4{2}]\t[%DP4{3}]\n'  >> ${snp_table}; then
                echo "$(date '+%Y-%m-%d %H:%M:%S'): Creating the Final SNP table for ${vcf} completed"
            else
                echo "$(date '+%Y-%m-%d %H:%M:%S'): bcf filter for ${tag} failed" && exit 1
            fi
        fi

        if [[ ${sv_call}="true" ]] ; then
            local manta_dir="${sample_dir}/minimap_dir/manta_dir"
            local manta_run="${manta_dir}/${tag}_manta_sv"
            local sv_vcf="${sample_dir}/minimap_dir/svs_vcf"
            local sv_table="${sample_dir}/minimap_dir/svs_tsv"
            
            # make sure all dirs exist
            create_dir ${manta_dir} ${sv_vcf} ${sv_table} "${manta_run}"

            # Structural variant calling with Manta
            logmsg "Running Manta for: ${tag} started"

            local vcfs=( "$manta_run"/results/variants/*.vcf.gz )

            if [[ -f "${vcfs[0]}" ]]; then
                logmsg "Manta SVs already exist for ${tag}, skipping."
            else
                logmsg "No SV VCFs found — running Manta."
                # Configure Manta and run manta
                if ! configManta.py --bam ${bam} --referenceFasta "${genome_fa}" \
                    --runDir ${manta_run} > "${manta_run}/configManta_${tag}.log" 2>&1; then
                    echo "$(date '+%Y-%m-%d %H:%M:%S'): Manta configuration failed for ${tag}."
                    exit 1
                fi

                if ! ${manta_run}/runWorkflow.py -m local -j ${threads} > "${manta_run}/mantaWorkflow_${tag}.log" 2>&1; then
                    echo "Manta workflow failed. Check log: ${manta_run}/mantaWorkflow_${tag}.log"
                    exit 1
                fi

            fi

            # copy VCFs
            local variant_dir="${manta_run}/results/variants"
            for svf in candidateSmallIndels.vcf.gz candidateSmallIndels.vcf.gz.tbi \
                candidateSV.vcf.gz candidateSV.vcf.gz.tbi diploidSV.vcf.gz diploidSV.vcf.gz.tbi; do  
                if [[ -e $svf ]] ; then  
                    logmsg "${variant_dir}/${svf} already exists—skipping copy."
                else
                    cp "${variant_dir}/${svf}" "${sv_vcf}/${tag}_${svf}" || {
                        echo "$(date '+%Y-%m-%d %H:%M:%S'): Failed to copy ${svf}"; exit 1; }
                fi
            done

            # Extract important SV fields
            for vf in ${sv_vcf}/*.vcf.gz; do
                sv_tsv=${sv_table}/$(basename ${vf} .vcf.gz).tsv
                logmsg "Extracting important SV info from: $(basename ${sv_tsv})"
                # Add headers and extract fields
                if [[ -s ${sv_tsv} ]] ; then
                    logmsg "${sv_tsv} already exists—skipping."
                else
                    echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tSVTYPE\tSVLEN\tEND" > ${sv_tsv}
                    # decompress .vcf.gz and Extract specific fields from the diplodVCF file
                    if bgzip -d -c ${vf} | grep -E '^#|^chr[1-9]\b|^chr10\b' | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\n'  >> ${sv_tsv}; then
                        echo "$(date '+%Y-%m-%d %H:%M:%S'): $(basename ${sv_tsv}) successfully created"
                    else
                        echo "$(date '+%Y-%m-%d %H:%M:%S'): $(basename ${sv_tsv}) failed"
                    fi
                fi
            done
        else
            logmsg "Structural Variant calling for ${tag} not needed. ...skipping."
        fi
    done
}



# before any samples…
for gx in "${goi[@]}"; do
  download_genome "$gx"
  index_genome   "$gx"
done


for sx in "${files[@]}"; do
    logmsg "trimming ${sx}"
    trim_reads "$sx"

    for gx in "${goi[@]}"; do
        map_reads "$sx" "$gx"
    done
done