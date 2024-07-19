mkdir -p /scratch/ldelisle/HoxBstudy/genomes/fasta/
cd /scratch/ldelisle/HoxBstudy/genomes/fasta/

# Get mm10
wget "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" -O mm10.fa.gz -nc
gunzip -k mm10.fa.gz

# Get the mutant genome mm10_HoxB_deli9-13insCBS5-10_del
# TODO: put on zenodo
aws s3 --profile perso --endpoint-url=https://s3.epfl.ch/ cp s3://11705-388fd8245175782087c769d3c1f8dabd/custom_genomes/mm10_deli9-13insCTCF_delB.fa ./mm10_HoxB_deli9-13insCBS5-10_del.fa

# Once zenodo published:
# # Get chr11_delB
# # Get chr11 mutant
# for i in {1..10} {12..19} X Y M; do
#   echo "chr${i}" >> listOfChrs.txt
# done
# seqtk seq -U mm10.fa.gz | seqtk subseq -l 60 - listOfChrs.txt | gzip > allChrsExceptchr11.fa.gz

# cat chr11_HoxBDel.fa.gz 'chr11_HoxBDel(i9-13):Ins(CBS5-10).fa.gz' allChrsExceptchr11.fa.gz > mm10_HoxB_deli9-13insCBS5-10_del.fa.gz
# gunzip -k mm10_HoxB_deli9-13insCBS5-10_del.fa.gz

mkdir -p /scratch/ldelisle/HoxBstudy/RNAseq/
cd /scratch/ldelisle/HoxBstudy/RNAseq/
wget "https://zenodo.org/records/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1" -O mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz -nc
gunzip -k mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
