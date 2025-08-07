#!/bin/bash

#SBATCH --nodes=1   # Number of nodes
#SBATCH --ntasks-per-node=32 # Number of processor cores per node
#SBATCH --output="a619.out"  # Job standard output file
#SBATCH --error="a619.error"  # Job standard error file
#SBATCH -p leiboff_lab   # specify the partition to submit the job to


# 0) Init conda
# Load Conda
source /local/cqls/software/x86_64/miniforge3/etc/profile.d/conda.sh

# Add necessary channels for bioinformatics tools
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels kbchoi


# Create a new Conda environment with specific packages
conda create --name a619 \
  python=3.6 htslib=1.12 sra-tools minimap2 vcftools bwa \
  samtools bcftools plink gatk4 picard gffread \
  matplotlib-base tectonic openjdk g2gtools pysam=0.15.3 -y 


# Activate the environment
conda activate a619

conda env export -n a619 > a619.yaml


# === SETTINGS ===
threads=32
chain="/nfs5/BPP/Leiboff_Lab/Brian/genomes/fake_genomes"
ref_dir="${chain}/ref_genomes"
merge_bam=true
do_download="true"

declare -A sample_loc=(
    [q1]="/nfs5/BPP/Leiboff_Lab/Brian/genomes/fake_genomes/a619_dir"
)

files=(q1)




# Genomes of interest
declare -A genome_urls=(
  [oh43]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0.fa.gz"
)

goi=("oh43")
# only one sample in this run
sra_a619=("SRR12633563" "SRR8907067" "SRR8997919" "SRR10127976" "SRR5725670" "SRR5663982" "SRR5663981")

###############################################################################

############### Hepler functions ###############
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

download_genome() {
  local gbase="$1"
  local url=${genome_urls[${gbase}]}               
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
            logmsg "Samtools faidx complete for $gbase"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S'): samtools index failed for ${gBASE}." &&  exit 1
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

        # merge that genome/sample
        if [[ "${merge_bam}" == true ]]; then
            local merged="${align_dir}/${tag}_merged.bam"
            local mbi="${merged}.bai"

            logmsg "Merging sorted BAMs into $(basename ${merged})"
            if [[ -s "${merged}" && -s "${mbi}" ]]; then
                logmsg "Merged BAM and index already exist. Skipping merge."
            else
                samtools merge -@ "${threads}" "${merged}" "${align_dir}/${tag}_"*.sorted.bam \
                    && logmsg "Created merged BAM: $(basename ${merged})" \
                    || { logmsg "samtools merge failed for ${tag}"; exit 1; }

                samtools index "${merged}" \
                    && logmsg "Created merged index: $(basename ${mbi})" \
                    || { logmsg "Indexing merged BAM failed"; exit 1; }
            fi
        fi
    done
}

# loop through the goi list
for genox in "${goi[@]}"; do
  download_genome "$genox"
  index_genome "$genox"
done

for sx in "${files[@]}"; do
    logmsg "trimming ${sx}"
    trim_reads "$sx"

    for gx in "${goi[@]}"; do
        map_reads "$sx" "$gx"
    done
done


conda deactivate

echo "[$(date '+%F %T')]: [DONE] Script completed successfully."