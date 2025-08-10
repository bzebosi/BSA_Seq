#!/bin/bash

#SBATCH --nodes=1   # Number of nodes
#SBATCH --ntasks-per-node=8 # Number of processor cores per node
#SBATCH --output="variant_call.out"  # Job standard output file
#SBATCH --error="variant_call.error"  # Job standard error file

# Define variables
env_name="variant_call_software"
script=/nfs5/BPP/Leiboff_Lab/Brian/scripts/BSA_Seq/linux_script/run_variant_call.sh


# Load Conda environment management
source /local/cqls/software/x86_64/miniforge3/etc/profile.d/conda.sh

# Ensure necessary Conda channels are present
for channel in defaults bioconda conda-forge; do
    conda config --show channels | grep -q "^${channel}$" || conda config --add channels "${channel}"
done

# Create Conda environment if it doesn't exist
if ! conda env list | grep -q "^${env_name} "; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Creating Conda environment '${env_name}'..."
    conda create --name "${env_name}" minimap2 vcftools samtools bcftools matplotlib-base python tectonic htslib seqkit gatk4 snpeff sift4g -y || {
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Error creating Conda environment '${env_name}'."
        exit 1
    }
else
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Conda environment '${env_name}' already exists."
fi

# Activate the environment
conda activate "${env_name}" || {
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Error: Failed to activate Conda environment '${env_name}'."
    exit 1
}

# Export environment details for reproducibility
conda env export -n "${env_name}" > "${env_name}.yaml"

# Run VCF processing script
if [ -x "${script}" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Running script: ${script}"
    "${script}" && echo "$(date '+%Y-%m-%d %H:%M:%S'): Successfully ran ${script}."
else
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Error: Script not found or not executable: ${script}"
    conda deactivate
    exit 1
fi

# Deactivate and clean up environment
conda deactivate

echo "$(date '+%Y-%m-%d %H:%M:%S'): Interactive SLURM job completed."
