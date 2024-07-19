#! /bin/bash
# Using pygenometracks version 3.9
gitHubDirectory=$PWD/
condaEnvName=lastVersion
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

pathWithAnnotations="${gitHubDirectory}/annotations"
publishedDataDirectory="${gitHubDirectory}/publishedData"

mkdir -p ${pathWithAnnotations}/outputs/

# Generate additional annotations:
cd ${pathWithAnnotations}/outputs/

# Custom function gtf2bed4
gtf2bed4() {
awk -F "\t" -v OFS="\t" '
BEGIN{
    start[""] = 0
}
{
    gene_name = ""
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    chrom[gene_name] = $1
    if (!(gene_name in start)) {
        start[gene_name] = $4
        end[gene_name] = $5
    }
    if (end[gene_name] < $5) {
        end[gene_name] = $5
    }
    if (start[gene_name] > $4) {
        start[gene_name] = $4
    }
}
END{
    for (gene_name in chrom) {
        if (gene_name != "") {
            print chrom[gene_name], start[gene_name] - 1, end[gene_name], gene_name
        }
    }
}' $1 | sort -k1,1 -k2,3n > $2
}

mkdir -p $publishedDataDirectory

# Get the gene annotation used in RNAseq:
wget "https://zenodo.org/records/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1" -O ${publishedDataDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz -nc

# Get protein codings in HoxB/HoxD region (removing "Hoxd3-203")
zcat ${publishedDataDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz | grep protein_coding | awk '($1=="chr11" && $4 < 97100000 && $5 > 95500000) || ($1=="chr2" && $4 < 75800000 && $5 > 73640001) {print}' | grep -v "Hoxd3-203" > protein_coding_around_HoxBD.gtf

# Adjust label on the left
sort -k1,1 -k4,4n protein_coding_around_HoxBD.gtf | awk -F "\t" '
BEGIN{
    HoxToShort="a5 b8 b7 b4 c6 c5 d10"
    split(HoxToShort,HoxArray, " ")
    for (i in HoxArray){
        HoxNamesToShort["Hox"HoxArray[i]] = 1
    }
}
{
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    if (gene_name in HoxNamesToShort) {
        gsub("Hox", "", gene_name)
    }
    if (! (gene_name in printed) ) {
        print $1"\t"$4"\t"$4+1"\t"gene_name
        printed[gene_name] = "Y"
    }
}' > protein_coding_around_HoxBD_left.bed

# Convert it to bed4
gtf2bed4 protein_coding_around_HoxBD.gtf protein_coding_around_HoxBD.bed4

# Get non coding around HoxB
zcat ${publishedDataDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz | awk '$1=="chr11" && $4 < 97100000 && $5 > 95500000 {print}' | grep -v protein_coding > non_protein_coding_around_HoxB.gtf

sort -k1,1 -k4,4n non_protein_coding_around_HoxB.gtf | awk -F "\t" '
{
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    if (! (gene_name in printed) ) {
        print $1"\t"$4"\t"$4+1"\t"gene_name
        printed[gene_name] = "Y"
    }
}' > non_protein_coding_around_HoxB_left.bed

# Select the 2 non coding genes of interest
grep "Gm53" non_protein_coding_around_HoxB.gtf > gm53.gtf
grep "Mir19" non_protein_coding_around_HoxB.gtf  > mir.gtf

echo -e "chr11\t122082543" > mm10.chr11.txt

# Enlarge mir 5x to be able to see it:
bedtools slop -b 2 -pct -i mir.gtf -g mm10.chr11.txt > mir_extended.gtf

# Get both nc into a single file
cat gm53.gtf mir_extended.gtf > selected_nc.gtf

# Convert all gtf to bed4
gtf2bed4 mir_extended.gtf mir_extended.bed4
gtf2bed4 gm53.gtf gm53.bed4
gtf2bed4 selected_nc.gtf selected_nc.bed4

grep "Gm53" non_protein_coding_around_HoxB_left.bed > gm53_left.bed
grep "Mir19" non_protein_coding_around_HoxB_left.bed > mir_left.bed

# To make prettier plots we remove 'Hox' and we remove '196a-1'.
cat protein_coding_around_HoxBD_left.bed gm53_left.bed mir_left.bed | sed "s/Hox//g" | sed "s/Mir196a-1/Mir/g" > selected_left.bed

cat selected_left.bed | awk -v OFS="\t" '
{
    if($4=="b13") {
        $4="Hoxb13"
    } else if($4=="Gm53") {
        $2-=13000
        $3-=13000
    } else if($4=="Mir") {
        $2-=5000
        $3-=5000
    } else if($4~/^b[1-8]/) {
        gsub("^b", "", $4)
    }
    print
}' > selected_left_5A.bed

cat selected_left.bed | awk -v OFS="\t" '
{
    if ($4=="Gm53") {
        $2-=13000
        $3-=13000
        print
    } else if($4=="Mir") {
        $2-=5000
        $3-=5000
        print
    } else if($4=="Evx2") {
        $2-=10000
        $3-=10000
        print
    }
}' > selected_left_nonHox.bed

cat selected_left.bed | awk -v OFS="\t" '
{
    if($4=="b13") {
        $4="Hoxb13"
        print
    } else if($4=="b9") {
        print
    } else if($4~/^b[1-8]/) {
        gsub("^b", "", $4)
        print
    } else if($4 == "d13") {
        $2-=7000
        $3-=7000
        print
    } else if($4~/^d/) {
        gsub("^d", "", $4)
        if ($4~/1[0-2]/) {
            $2-=2000
            $3-=2000
        }
        print
    }
}' > selected_left_Hox.bed

# Get genes present in delHoxd10-12
bedtools intersect -a protein_coding_around_HoxBD.bed4 -b ${pathWithAnnotations}/mice_HoxD_deletion.bed -v -wa > protein_coding_around_HoxBD_del1012.bed4

touch empty.bed

# Get a simplified version of selected_left.bed for zoom out:
cat selected_left.bed | sort -k1,1 -k2,2n | awk -v OFS="\t" -v shift=20000 '
{
    if ($4 == "b13") {
        $2 -= shift
        $3 -= shift
        print
    } else if ($4 ~ "^b[1-9]") {
        if ($4 == "b9") {
            $4="b9-b1"
            print
        }
    }
}' > selected_left_grouped.bed

# For zoomed figures
cat ${pathWithAnnotations}/CTCF_inHox_colored.bed | awk -v shift=6000 -v OFS="\t" '{split($4, a, "_");cshift = shift; if (a[1] != "CBS10") {gsub("CBS", "", a[1]); cshift = shift / 3}; print $1,$2-cshift,$3-cshift,a[1]}' > CTCF_inHox_label_zoom.bed
cat ${pathWithAnnotations}/CTCF_inHox_colored.bed | awk -v shift=3000 -v OFS="\t" '{split($4, a, "_");cshift = shift; gsub("CBS", "", a[1]); print $1,$2-cshift,$3-cshift,a[1]}' > CTCF_inHox_label_zoom2.bed

# Get all mus gtfs:
mkdir -p mus_gtf
cd mus_gtf
wget https://ftp.ensembl.org/pub/release-102/gtf/ -O ${publishedDataDirectory}/102_gtf -nc
musDirs=$(awk -F "\"" '{for(i = 1;i<=NF;i++){if ($i~/^mus_/){print $i}}}' ${publishedDataDirectory}/102_gtf)
for musDir in $musDirs; do
    name=$(basename $musDir)
    wget https://ftp.ensembl.org/pub/release-102/gtf/$musDir -O ${publishedDataDirectory}/102_gtf_${name} -nc
    gtfName=$(awk -F "\"" '{for(i = 1;i<=NF;i++){if ($i~/.chr.gtf.gz$/){print $i}}}' ${publishedDataDirectory}/102_gtf_${name})
    wget https://ftp.ensembl.org/pub/release-102/gtf/${musDir}$gtfName -P ${publishedDataDirectory}/ -nc
done
# Generate fake gtfs where Hoxb9 is centered at 0.5e6
echo -e "11\t122082543" > mm10.11.txt
for gtf in ${publishedDataDirectory}/Mus*chr.gtf.gz; do
    echo $gtf
    name=$(basename $gtf | awk -F "\\." '{print $1}')
    if [ ! -e ${name}_Hoxb9_centered.gtf ]; then
        hoxb9mid=$(zcat $gtf | grep Hoxb9 | awk '$3=="gene"{printf("%d\n", ($4+$5)/2)}')
        if [ $name = "Mus_pahari" ]; then
            echo -e "14\t$((hoxb9mid - 5000000))\t$((hoxb9mid + 5000000))" > temp.bed
            bedtools intersect -a $gtf -b temp.bed -wa | awk -F "\t" -v OFS="\t" '{$1=11;print}' | bedtools shift -i - -g mm10.11.txt -s $((5000000 - hoxb9mid)) > ${name}_Hoxb9_centered.gtf
        else
            echo -e "11\t$((hoxb9mid - 5000000))\t$((hoxb9mid + 5000000))" > temp.bed
            bedtools intersect -a $gtf -b temp.bed -wa | bedtools shift -i - -g mm10.11.txt -s $((5000000 - hoxb9mid)) > ${name}_Hoxb9_centered.gtf
        fi
        rm temp.bed
    fi
    grep protein_coding ${name}_Hoxb9_centered.gtf > ${name}_Hoxb9_centered_pc.gtf
    grep Gm53 ${name}_Hoxb9_centered.gtf > ${name}_Hoxb9_centered_npc_selected.gtf
    grep "Mir19" ${name}_Hoxb9_centered.gtf >> ${name}_Hoxb9_centered_npc_selected.gtf
    gtf2bed4 ${name}_Hoxb9_centered_pc.gtf ${name}_Hoxb9_centered_pc.bed4
    gtf2bed4 ${name}_Hoxb9_centered_npc_selected.gtf ${name}_Hoxb9_centered_npc_selected.bed4
done

# Get the order:
for f in *_Hoxb9_centered_pc.gtf; do
    hoxb13start=$(grep Hoxb13 $f | head -n 1 | cut -f 4)
    echo $(basename $f _Hoxb9_centered_pc.gtf) $hoxb13start
done | sort -k2,2n > Hoxb13_order.txt

# Get labels:
first_name=$(head -n 1 Hoxb13_order.txt | cut -f 1 -d " ")
cat ${first_name}_Hoxb9_centered_pc.gtf ${first_name}_Hoxb9_centered_npc_selected.gtf | sort -k1,1 -k4,4n | awk -F "\t" '
BEGIN{
    HoxToShort="b8 b7 b4"
    split(HoxToShort,HoxArray, " ")
    for (i in HoxArray){
        HoxNamesToShort["Hox"HoxArray[i]] = 1
    }
}
{
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    if (gene_name in HoxNamesToShort) {
        gsub("Hox", "", gene_name)
    }
    if (gene_name~/Mir/) {
        gene_name = "Mir"
    }
    if (! (gene_name in printed) ) {
        print $1"\t"$4"\t"$4+1"\t"gene_name
        printed[gene_name] = "Y"
    }
}' > selected_left_Hoxb9_centered.bed
#### HERE ####
name=Mus_musculus
cat ${name}_Hoxb9_centered_pc.gtf ${name}_Hoxb9_centered_npc_selected.gtf | sort -k1,1 -k4,4n | awk -F "\t" '
BEGIN{
    HoxToShort="b8 b7 b4"
    split(HoxToShort,HoxArray, " ")
    for (i in HoxArray){
        HoxNamesToShort["Hox"HoxArray[i]] = 1
    }
}
{
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    if (gene_name in HoxNamesToShort) {
        gsub("Hox", "", gene_name)
    }
    if (gene_name~/Mir/) {
        gene_name = "Mir"
    }
    if (! (gene_name in printed) ) {
        print $1"\t"$4"\t"$4+1"\t"gene_name
        printed[gene_name] = "Y"
    }
}' > ${name}_selected_left_Hoxb9_centered.bed
# wget https://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz -O ${publishedDataDirectory}/MGI.gff3.gz -nc

# These ones are mm39
# mm10 chr11:96183548-96376004 -> mm39 chr11:96074374-96266830
for i in {1..4}; do
    wget https://ftp.ensembl.org/pub/release-111/maf/ensembl-compara/multiple_alignments/21_murinae.epo/21_murinae.epo.11_${i}.maf.gz -nc -P ${publishedDataDirectory}/
done
gunzip -k ${publishedDataDirectory}/21_murinae.epo.11_*.maf.gz

# Try to reduce to 11:96074374-96266830
echo -e "11\t96074374\t96266830" > mm39_region.bed

# First index all:
for f in ${publishedDataDirectory}/21_murinae.epo.11_*.maf; do
    maf_build_index.py -s mus_musculus $f
done
# Then extract
for f in ${publishedDataDirectory}/21_murinae.epo.11_*.maf; do
    maf_extract_ranges_indexed.py -c -p mus_musculus. $f < mm39_region.bed > $(basename $f .maf)_HoxB.maf
done
# Only the number 3 has a block
maf_build_index.py -s mus_musculus 21_murinae.epo.11_3_HoxB.maf

# Cut them into chunk
# Remove ancestral sequences
# and put them back together
mkdir -p chunk
chunksize=1000
current_start=$(cut -f 2 mm39_region.bed)
current_end=$((current_start + chunksize))
end=$(cut -f 3 mm39_region.bed)
i=1
while [ $current_end -lt $end ]; do
    echo $current_start $current_end
    echo -e "11\t${current_start}\t${current_end}" > temp.bed
    maf_extract_ranges_indexed.py -c -p mus_musculus. 21_murinae.epo.11_3_HoxB.maf < temp.bed > chunk/chunk_${i}.maf
    if [ $i = "1" ]; then
        grep -v "ancestral_sequences" chunk/chunk_${i}.maf > 21_murinae.epo.11_3_HoxB_chunked.maf
    else
        grep -v "^#" chunk/chunk_${i}.maf | grep -v "ancestral_sequences" >> 21_murinae.epo.11_3_HoxB_chunked.maf
    fi
    current_start=$((current_start + chunksize))
    current_end=$((current_start + chunksize))
    i=$((i + 1))
done
if [ $current_start != $end ]; then
    echo -e "11\t${current_start}\t${end}" > temp.bed
    maf_extract_ranges_indexed.py -c -p mus_musculus. 21_murinae.epo.11_3_HoxB.maf < temp.bed > chunk/chunk_${i}.maf
    grep -v "^#" chunk/chunk_${i}.maf | grep -v "ancestral_sequences" >> 21_murinae.epo.11_3_HoxB_chunked.maf
fi
rm temp.bed
rm -r chunk
# Remove the intermediate mafs
rm 21_murinae.epo.11_[1-4]_HoxB.maf*


# I need annotations for this mm39:
wget "https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.chr.gtf.gz" -P ${publishedDataDirectory}/
bedtools intersect -a ${publishedDataDirectory}/Mus_musculus.GRCm39.111.chr.gtf.gz -b mm39_region.bed -wa > 21_murinae.epo.11_3_HoxB_mus_musculus.gtf
grep protein_coding 21_murinae.epo.11_3_HoxB_mus_musculus.gtf > 21_murinae.epo.11_3_HoxB_mus_musculus_pc.gtf
grep Gm53 21_murinae.epo.11_3_HoxB_mus_musculus.gtf > 21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.gtf
grep "Mir19" 21_murinae.epo.11_3_HoxB_mus_musculus.gtf >> 21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.gtf

cat 21_murinae.epo.11_3_HoxB_mus_musculus_pc.gtf 21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.gtf | sort -k1,1 -k4,4n | awk -F "\t" '
BEGIN{
    HoxToShort="b8 b7 b4"
    split(HoxToShort,HoxArray, " ")
    for (i in HoxArray){
        HoxNamesToShort["Hox"HoxArray[i]] = 1
    }
}
{
    split($9, a, "\"")
    for (i=1;i<=length(a);i++) {
        if (a[i]~/gene_name/) {
            gene_name=a[i + 1]
        }
    }
    if (gene_name in HoxNamesToShort) {
        gsub("Hox", "", gene_name)
    }
    if (gene_name~/Mir/) {
        gene_name = "Mir"
    }
    if (! (gene_name in printed) ) {
        print $1"\t"$4"\t"$4+1"\t"gene_name
        printed[gene_name] = "Y"
    }
}' > 21_murinae.epo.11_3_HoxB_mus_musculus_selected_left.bed

gtf2bed4 21_murinae.epo.11_3_HoxB_mus_musculus_pc.gtf 21_murinae.epo.11_3_HoxB_mus_musculus_pc.bed4
gtf2bed4 21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.gtf 21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.bed4

# The maf for 60 species (mm10)
wget "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/multiz60way/maf/chr11.maf.gz" -O ${publishedDataDirectory}/chr11.maf.gz -nc
gunzip -k ${publishedDataDirectory}/chr11.maf.gz
maf_build_index.py -s mm10 ${publishedDataDirectory}/chr11.maf

# Prepare annotations for mutants
inputCTCF=${pathWithAnnotations}/CTCF_inHox_colored.bed
inputGenesPC=protein_coding_around_HoxBD.gtf
inputSelectedNPC=selected_nc.gtf
inputLabels=selected_left.bed
inputCTCFcassette=${pathWithAnnotations}/CTCF_cassette.bed
## First for the Deli9-13
br_file=${gitHubDirectory}/prepare/BR-file-deli9-13.txt
genomeName=$(cat $br_file | cut -f 1 | tail -n 1)
genome=${genomeName}

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_compatibleInv_splitIfOV_chr11.R $br_file ${inputCTCF} 1 2 3 6 TRUE ${genome}_CTCF.bed
mkdir -p ${pathWithAnnotations}/${genome}/
mv ${gitHubDirectory}/prepare/${genomeName}/${genomeName}_vSplit_${genome}_CTCF.bed ${pathWithAnnotations}/${genome}/${genome}_CTCF.bed
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" -v colorOther="0,0,0" '
{
    if ($6 == "+"){
        color = colorPos
    } else if ($6 == "-" ){
        color = colorNeg
    } else {
        color = colorOther
    }
    print $1, $2, $3, $4, $5, $6, $2, $3, color
}' ${pathWithAnnotations}/${genome}/${genome}_CTCF.bed > ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_compatibleInv_splitIfOV_chr11.R $br_file ${inputGenesPC} 1 4 5 7 FALSE ${genome}_protein_coding_around_HoxBD.gtf
mv ${gitHubDirectory}/prepare/${genomeName}/${genomeName}_vSplit_${genome}_protein_coding_around_HoxBD.gtf ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.gtf

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_compatibleInv_splitIfOV_chr11.R $br_file ${inputSelectedNPC} 1 4 5 7 FALSE ${genome}_selected_nc.gtf
mv ${gitHubDirectory}/prepare/${genomeName}/${genomeName}_vSplit_${genome}_selected_nc.gtf ${pathWithAnnotations}/${genome}/${genome}_selected_nc.gtf

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_compatibleInv_splitIfOV_chr11.R $br_file ${inputLabels} 1 2 3 6 TRUE ${genome}_selected_left.bed

cat ${gitHubDirectory}/prepare/${genomeName}/${genomeName}_vSplit_${genome}_selected_left.bed | awk -v OFS="\t" '
{
    if($4=="b13") {
        $4="13"
        $2-=4000
        $3-=4000
    } else if($4=="Gm53") {
        $2-=13000
        $3-=13000
    } else if($4=="Mir") {
        $2-=4000
        $3-=4000
    } else if($4~/^b[1-8]/) {
        gsub("^b", "", $4)
    }
    print
}' > ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed

rm ${gitHubDirectory}/prepare/${genomeName}/${genomeName}_vSplit_${genome}_selected_left.bed

gtf2bed4 ${pathWithAnnotations}/${genome}/${genome}_selected_nc.gtf ${pathWithAnnotations}/${genome}/${genome}_selected_nc.bed4
gtf2bed4 ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.gtf ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4

# Shift annotations
br_file=${gitHubDirectory}/prepare/BR-file-deli9-13insCBS5-10.txt
genomeName=$(cat $br_file | cut -f 1 | tail -n 1)
genome=$genomeName
mkdir -p ${pathWithAnnotations}/${genome}/
Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R ${inputCTCF} ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 2 3 6 TRUE ${pathWithAnnotations}/${genome}/${genome}_CTCF.bed
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" -v colorOther="0,0,0" '
{
    if ($6 == "+"){
        color = colorPos
    } else if ($6 == "-" ){
        color = colorNeg
    } else {
        color = colorOther
    }
    print $1, $2, $3, $4, $5, $6, $2, $3, color
}' ${pathWithAnnotations}/${genome}/${genome}_CTCF.bed > ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R ${inputCTCFcassette} ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 2 3 6 TRUE ${genome}_CTCF_cassette.bed
# Fix the thickStart and thickEnd
awk -v OFS="\t" '{
    $7 = $2
    $8 = $3
    print
}' ${genome}_CTCF_cassette.bed > ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette.bed
rm ${genome}_CTCF_cassette.bed

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R  ${inputGenesPC} ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 4 5 7 FALSE ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.gtf

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R  ${inputSelectedNPC} ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 4 5 7 FALSE ${pathWithAnnotations}/${genome}/${genome}_selected_nc.gtf

# Remove Gm53
grep -v "Gm53" ${pathWithAnnotations}/${genome}/${genome}_selected_nc.gtf > ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.gtf

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R ${inputLabels} ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 2 3 0 TRUE ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed

# Put the full b names
sed -i 's/\tb/\tHoxb/g' ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed
# Remove Gm and Mir
grep -v "Gm" ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed | grep -v "Mir" > ${pathWithAnnotations}/${genome}/${genome}_selected_left_pc.bed

gtf2bed4 ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.gtf ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.bed4
gtf2bed4 ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.gtf ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4

Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R ${pathWithAnnotations}/CTCF_cassette_breakpoints.bed ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 2 3 0 TRUE ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette_breakpoints.bed
cat ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette_breakpoints.bed | awk -v OFS="\t" -v prefix=${pathWithAnnotations}/${genome}/${genome}_ '
$4~/_start/{
  gsub("_start", "", $4)
  start[$4] = $2
}
$4~/_end/{
  gsub("_end", "", $4)
  $2 = start[$4]
  print > prefix $4 ".bed"
}'
# Get left CTCF:
for f in ${pathWithAnnotations}/${genome}/${genome}_*_del.bed; do
    bedtools intersect -a ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed \
        -b $f -v > ${f/.bed/_CTCF_left.bed}
done

# For zoomed figures
cat ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed | awk '{split($4, a, "_");cshift = 1700; if (a[1] != "CBS10") {gsub("CBS", "", a[1]); cshift = 100}; printf("%s\t%d\t%d\t%s\n",$1,$2-cshift,$3-cshift,a[1])}' > ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed

# Get CRISPR positions for nCATS
Rscript ${gitHubDirectory}/prepare/scripts/shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R ${pathWithAnnotations}/CRISPR_nCATS.bed ${gitHubDirectory}/prepare/${genome}.txt $br_file 1 2 3 0 TRUE ${pathWithAnnotations}/${genome}/${genome}_CRISPR_nCATS.bed


# Get repeat mask
# mm10_rmsk.bed.gz
# which can be obtained through UCSC website
# in tools > table browser > variations and repeat.
# And put it into $publishedDataDirectory
