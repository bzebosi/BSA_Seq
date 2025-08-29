#!/bin/bash

#SBATCH --nodes=1   # Number of nodes
#SBATCH --ntasks-per-node=16 # Number of processor cores per node
#SBATCH --output="variant_call.out"  # Job standard output file
#SBATCH --error="variant_call.error"  # Job standard error file
#SBATCH -p leiboff_lab   # specify the partition to submit the job to

# Define variables
env_name="variant_call_packages"
script=/nfs5/BPP/Leiboff_Lab/Brian/scripts/BSA_Seq/linux_script/run_variant_call.sh


# Load Conda environment management
source /local/cqls/software/x86_64/miniforge3/etc/profile.d/conda.sh

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
logmsg() { echo "$(timestamp): $*"; }

# Ensure necessary Conda channels are present
for channel in defaults bioconda conda-forge; do
    conda config --show channels | grep -q "^${channel}$" || conda config --add channels "${channel}"
done

# Create Conda environment if it doesn't exist
if ! conda env list | grep -q "^${env_name} "; then
    logmsg "Creating Conda environment '${env_name}'..."
    conda create --name "${env_name}" minimap2 vcftools samtools bcftools matplotlib-base python tectonic htslib seqkit gatk4 snpeff sift4g -y || {
        logmsg "Error creating Conda environment '${env_name}'."
        exit 1
    }
else
    logmsg "Conda environment '${env_name}' already exists."
fi

# Activate the environment
conda activate "${env_name}" || {
    logmsg "Error: Failed to activate Conda environment '${env_name}'."
    exit 1
}

# Export environment details for reproducibility
conda env export -n "${env_name}" > "${env_name}.yaml"

# Run VCF processing script
if [ -x "${script}" ]; then
    logmsg "Running script: ${script}"
    "${script}" && logmsg "Successfully ran ${script}."
else
    logmsg "Error: Script not found or not executable: ${script}"
    conda deactivate
    exit 1
fi

# Deactivate and clean up environment
conda deactivate

logmsg "SLURM job completed."
