#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the number of samples
#SBATCH --cpus-per-task 1 # This does not use multiple CPUs
#SBATCH --time 02:00:00 # This depends on the number of samples
#SBATCH --job-name RNAseq_R_DESeq2 # Job name that appear in squeue as well as in output and error text files

# This script run DESeq2 analysis


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the paths
# Put in dirPathWithConfigFiles the directory
# where all config files will be written
# this can be on your home
dirPathWithConfigFiles="$PWD/"
# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample (should match the one in 01_RNAseq_XX.sh)
dirPathWithResults="/scratch/ldelisle/HoxBstudy/RNAseq/"
filePathForGTF="${dirPathWithResults}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf"
# This is where should be clone the 2 github repositories:
# https://github.com/lldelisle/toolBoxForMutantAndWTGenomes
# https://github.com/lldelisle/rnaseq_rscripts
dirPathWithDependencies="/home/ldelisle/softwares/"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/8.2.0-2.31.1  
# module load OpenMPI/3.1.3
# module load R/3.6.0

# You can choose to use a conda environment to solve R and associated packages
# You can create it with: conda create -n RNAseq_R_202302 r-base r-colorspace bioconductor-deseq2 r-ggplot2 r-pheatmap r-rcolorbrewer  bioconductor-rtracklayer
# Comment it if you will use module load
condaEnvName=RNAseq_R_202302
######


##################################
####### BEGINING OF SCRIPT #######
##################################

# Check everything is set correctly:
if [ ! -z ${condaEnvName} ]; then
    # Conda environment:
    # This line is to adapt the conda to the shell
    source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
    # We check if the conda environment exists
    exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
    # It if does not exists an error is raised
    if [ $exists -ne 1 ]; then
    echo "conda environment ${condaEnvName} does not exists. Create it before."
    exit 1
    fi
    # Activate the conda environment
    conda activate ${condaEnvName}
fi

# Check all softwares are present and write version to stdout:
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
# Check  are installed:
Rscript -e "library(colorspace);library(DESeq2);library(ggplot2);library(pheatmap);library(RColorBrewer);library(rtracklayer)"
if [ $? -ne 0 ]
then
  echo "Some R packages are missing check rtracklayer, colorspace, DESeq2, ggplot2, pheatmap and RColorBrewer are installed."
  exit 1
fi

# Check the 2 github:
if [ ! -e ${dirPathWithDependencies}rnaseq_rscripts/ ]; then
  echo "${dirPathWithDependencies}rnaseq_rscripts/ does not exists please clone https://github.com/lldelisle/rnaseq_rscripts"
  exit 1
fi
cd ${dirPathWithDependencies}rnaseq_rscripts/
echo "Version of rnaseq_rscripts"
git rev-parse HEAD

if [ ! -e ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/ ]; then
  echo "${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/ does not exists please clone https://github.com/lldelisle/toolBoxForMutantAndWTGenomes"
  exit 1
fi
cd ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/
echo "Version of toolBoxForMutantAndWTGenomes"
git rev-parse HEAD


## START
cd "${dirPathWithConfigFiles}"
# Extract protein coding:
for dir in ${dirPathWithResults}/mergedTables/ ${dirPathWithResults}/mergedTablesE/; do
  if [ ! -e ${dir}/AllHTSeqCounts_subset_onlypc.txt ]; then
    Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/subsetForProteinCoding.R "${filePathForGTF}" ${dir}/AllHTSeqCounts_subset.txt Ens_ID
  fi
  if [ ! -e ${dir}/AllCufflinks_Simplified_subset_onlypc.txt ]; then
    Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/subsetForProteinCoding.R "${filePathForGTF}" ${dir}/AllCufflinks_Simplified_subset.txt gene_id
  fi
done
for configFile in configFileRNAseq_step2multi.R configFileRNAseq_step2multi_E.R; do
    Rscript ${dirPathWithDependencies}rnaseq_rscripts/step2-multi_DESeq2.R $configFile
done
Rscript Get_summary_some_genes.R configFileRNAseq_step2multi.R
