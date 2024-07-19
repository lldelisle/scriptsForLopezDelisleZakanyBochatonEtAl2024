#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the number of samples
#SBATCH --cpus-per-task 1 # This does not use multiple CPUs
#SBATCH --time 02:00:00 # This depends on the number of samples
#SBATCH --job-name RNAseq_R # Job name that appear in squeue as well as in output and error text files

# This script run Rscripts to merge counts and FPKM
# Generate first plots
# And DESeq2 analysis


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
# A samples plan is a tabular table with at least one column named 'sample'
# This should match the table of 01_RNAseq_XX.sh
# The sample names should not start with letter and not contain parenthesis or -
filePathForSamplesPlan="$PWD/samplesplan_embryos.txt"
# Sometimes some chromosomes are excluded from analysis:
chrsToRemove="chrX,chrY,chrM"
# chrsToRemove="chrM"

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

# This script will merge the tables

# First generate tables
if [ ! -e ${dirPathWithResults}/mergedTablesE/AllCufflinks_Simplified_subset.txt ]; then
  configFile="configFileRNAseq_step1E.R"

  echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- \"${dirPathWithResults}/mergedTablesE/\"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- T
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- \"$PWD/genes${chrsToRemove//,/_}.txt\"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- T
chrToRemove <- c(\"${chrsToRemove//,/\",\"}\")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- T
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')).  If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c(\"${chrsToRemove//,/\",\"}\")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
" > ${configFile}

  if [ ! -e genes${chrsToRemove//,/_}.txt ]; then
    Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $filePathForGTF ${chrsToRemove} ./
    mv genesIn${chrsToRemove}from* genes${chrsToRemove//,/_}.txt
  fi
  # Adjust the samplesplan
  if [ ! $(grep "htseq_count_file" $filePathForSamplesPlan) ]; then
    cat $filePathForSamplesPlan | awk -v pa=$dirPathWithResults 'BEGIN{print "htseq_count_file"}NR>1{print pa"/allFinalFiles/counts_FPKM/htseqCount_"$1".txt"}' > htseqCol.txt
    paste -d "\t" $filePathForSamplesPlan htseqCol.txt > ${filePathForSamplesPlan}_withPaths
    rm htseqCol.txt
  fi
  if [ ! $(grep "cufflinks_file" $filePathForSamplesPlan) ]; then
    cat $filePathForSamplesPlan | awk -v pa=$dirPathWithResults 'BEGIN{print "cufflinks_file"}NR>1{print pa"/allFinalFiles/counts_FPKM/FPKM_"$1".txt"}' > CuffCol.txt
    if [ -e ${filePathForSamplesPlan}_withPaths ];then
      mv ${filePathForSamplesPlan}_withPaths tmp
      paste -d "\t" tmp CuffCol.txt > ${filePathForSamplesPlan}_withPaths
      rm tmp
    else
      paste -d "\t" $filePathForSamplesPlan CuffCol.txt > ${filePathForSamplesPlan}_withPaths
    fi
    rm CuffCol.txt
  fi
  if [ -e ${filePathForSamplesPlan}_withPaths ]; then
    cat $configFile | sed "s#$filePathForSamplesPlan#${filePathForSamplesPlan}_withPaths#" > ${configFile}_withPaths
    configFile=${configFile}_withPaths
  fi

  Rscript ${dirPathWithDependencies}rnaseq_rscripts/step1-generateTables.R $configFile
  # copy to GEO:
  mkdir -p ${dirPathWithResults}/toGEO
  cp ${dirPathWithResults}/mergedTablesE/AllCufflinks_Simplified.txt ${dirPathWithResults}/toGEO/AllCufflinks_Simplified_E.txt
  cp ${dirPathWithResults}/mergedTablesE/AllHTSeqCounts.txt ${dirPathWithResults}/toGEO/AllHTSeqCounts_E.txt
  mkdir -p ${dirPathWithConfigFiles}/outputs
  gzip -c ${dirPathWithResults}/mergedTablesE/AllCufflinks_Simplified.txt > ${dirPathWithConfigFiles}/outputs/AllCufflinks_Simplified_E.txt.gz
  gzip -c ${dirPathWithResults}/mergedTablesE/AllHTSeqCounts.txt > ${dirPathWithConfigFiles}/outputs/AllHTSeqCounts_E.txt.gz
fi
