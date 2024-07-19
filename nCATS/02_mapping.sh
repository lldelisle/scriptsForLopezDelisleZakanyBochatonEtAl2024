#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 10
#SBATCH --time 24:00:00
#SBATCH --array=0-1
#SBATCH --job-name nCATS_mapping
#SBATCH --chdir /scratch/ldelisle/HoxBstudy/nCATS/

gitHubDirectory=$1
genomeDir=$2

genomes=('mm10' 'mm10_HoxB_deli9-13insCBS5-10_del')

path="$PWD/"
sample=delBins5-10
nThreads=10

genome=${genomes[${SLURM_ARRAY_TASK_ID}]}
pathForFasta="${genomeDir}/fasta/${genome}.fa"
pathForMinimap2Index="${genomeDir}/minimap2/${genome}.mmi"

mkdir -p $(dirname $pathForMinimap2Index)

condaEnvName=atac202209
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

pathResults=${path}/${sample}_${genome}/

mkdir -p $pathResults
echo $sample
cd $pathResults

mkdir -p combinedFastq
fastq="${pathResults}/combinedFastq/${sample}.fastq.gz"
if [ -e $fastq ]; then
  echo "The fastq already exists"
else
  cat ${path}/guppy_output/fastq_*.fastq | gzip > $fastq
fi

if [ ! -e $pathForMinimap2Index ]; then
    ~/softwares/minimap2-2.28_x64-linux/minimap2 -t $nThreads -d $pathForMinimap2Index $pathForFasta
fi

if [ ! -e ${sample}_algn_splice.sam ]; then
  ~/softwares/minimap2-2.28_x64-linux//minimap2 -t $nThreads -ax splice $pathForMinimap2Index $fastq > ${sample}_algn_splice.sam
fi

# Due to a design flaw, BAM does not work with CIGAR strings with >65535 operations (SAM and CRAM work). 
# So this may not work
# Here it worked
if [ ! -e ${sample}_mapped_splice_sorted.bam ]; then
  samtools sort --threads $nThreads -o ${sample}_mapped_splice_sorted.bam ${sample}_algn_splice.sam
fi

# We remove the non primary alignments
if [ ! -e  ${sample}_mapped_splice_sorted_filtered2.bam.bai ]; then
  samtools view -b ${sample}_mapped_splice_sorted.bam --threads $nThreads -F0x104 > ${sample}_mapped_splice_sorted_filtered2.bam
  samtools index ${sample}_mapped_splice_sorted_filtered2.bam
fi

# We generate coverage
if [ ! -e ${sample}_splice_filtered2_cov.bedGraph.gz ]; then
  bedtools genomecov -ibam ${sample}_mapped_splice_sorted_filtered2.bam -bg -split | LC_ALL=C sort -k1,1 -k2,2n > ${sample}_splice_filtered2_cov.bedGraph
  bedGraphToBigWig ${sample}_splice_filtered2_cov.bedGraph ${pathForFasta}.fai ${sample}_splice_filtered2_cov.bw
  gzip ${sample}_splice_filtered2_cov.bedGraph
fi

# Generate bed12
if [ ! -e ${sample}_mapped_splice_sorted_filtered2.bed12.gz ]; then
  bedtools bamtobed -bed12 -i ${sample}_mapped_splice_sorted_filtered2.bam | gzip > ${sample}_mapped_splice_sorted_filtered2.bed12.gz
fi

# Select bed12
if [ ! -e ${sample}_HoxBdel.bed ]; then
  zcat ${sample}_mapped_splice_sorted_filtered2.bed12.gz | awk '
NR==FNR{
  reads[$0] = 1
}
NR!=FNR{
  if ($4 in reads) {
    print $0
  }
}' ${gitHubDirectory}/nCATS/selected_HoxBdel.txt - > ${sample}_HoxBdel.bed
fi
if [ ! -e ${sample}_cassette.bed ]; then
  zcat ${sample}_mapped_splice_sorted_filtered2.bed12.gz | awk '
NR==FNR{
  reads[$0] = 1
}
NR!=FNR{
  if ($4 in reads) {
    print $0
  }
}' ${gitHubDirectory}/nCATS/selected_cassette.txt - > ${sample}_cassette.bed
fi

mkdir -p ../toGEO
cp ${sample}_splice_filtered2_cov.bw ../toGEO/nCATS_${sample}_on${genome}.bw
mkdir -p ${gitHubDirectory}/nCATS/output/${genome}
cp ${sample}*bed ${gitHubDirectory}/nCATS/output/${genome}/
