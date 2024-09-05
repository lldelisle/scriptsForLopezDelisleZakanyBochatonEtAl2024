#! /bin/bash
# Using pygenometracks version 3.9
gitHubDirectory=$PWD/
GEODirectory=~/mnt/scratch/HoxBstudy/toGEO/

condaEnvName=lastVersion # 3.9
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

pathWithAnnotations="${gitHubDirectory}/annotations"
pathWithCR="${GEODirectory}/CUTandRUN"
pathWithChIPM="${GEODirectory}/ChIPM"
pathWithnCATS="${GEODirectory}/nCATS"
plottingDirectory="${gitHubDirectory}/plots/outputs"
publishedDataDirectory="${gitHubDirectory}/publishedData"

mkdir -p $publishedDataDirectory

mkdir -p $plottingDirectory

cd $plottingDirectory

# Figure 1A:
ini_file="Fig1A.ini"
echo "[scalebar]
file_type = scalebar
height = 0.14
fontsize = 6

[CTCF label]
file = ${pathWithAnnotations}/outputs/CTCF_inHox_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
arrowhead_fraction = 0.01
" > ${ini_file}
for f in ${pathWithAnnotations}/outputs/empty.bed ${pathWithAnnotations}/mice_deletion.bed ${pathWithAnnotations}/mice_rescue.bed ${pathWithAnnotations}/mice_b13hd.bed; do
    echo "[$(basename $f)]
file = $f
display = deletions
color = red
labels = false
height = 0.5

[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
overlay_previous = share-y
color = black
display = collapsed
color_utr = black
border_color = none
labels = false

[miR]
file = ${pathWithAnnotations}/outputs/mir_extended.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
labels = false
" >> ${ini_file}
    if [[ $(basename $f) == "empty.bed" || $(basename $f) = "mice_b13hd.bed" ]]; then
        echo "[gm]
file = ${pathWithAnnotations}/outputs/gm53.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false
" >> ${ini_file}
    fi
done
echo "[genes label Hox]
file = ${pathWithAnnotations}/outputs/selected_left_Hox.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 5
height = 0.1

[genes label non-Hox]
file = ${pathWithAnnotations}/outputs/selected_left_nonHox.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 4.5
overlay_previous = share-y
" >> ${ini_file}
pgt --region chr11:96190000-96370000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 5.5


# Figure 1B:
ini_file="Fig1B.ini"
echo "[scalebar]
file_type = scalebar
height = 0.14
fontsize = 6

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
arrowhead_fraction = 0.01
" > ${ini_file}
for f in ${pathWithAnnotations}/outputs/empty.bed ${pathWithAnnotations}/mice_HoxD_deletion.bed ${pathWithAnnotations}/outputs/empty.bed; do
    echo "[$(basename $f)]
file = $f
display = deletions
color = red
labels = false
height = 0.5
" >> ${ini_file}
    if [[ $f = *"deletion.bed" ]]; then
        echo "[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD_del1012.bed4
overlay_previous = share-y
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
" >> ${ini_file}
    else
        echo "[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
overlay_previous = share-y
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
" >> ${ini_file}
    fi
done
echo "[genes label Hox]
file = ${pathWithAnnotations}/outputs/selected_left_Hox.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 5
height = 0.1

[genes label non-Hox]
file = ${pathWithAnnotations}/outputs/selected_left_nonHox.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 4.5
overlay_previous = share-y
" >> ${ini_file}
# Same scale as 1A: 180kb
pgt --region chr2:74,630,000-74,810,000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 5.5

# Get all published data needed for Fig5A
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226299&format=file&file=GSM6226299%5Fwt%5F120h%5FH3K27ac%5FHiChIP%5F10kb%5FRaw%5Fmm10%2Ecool%2Egz" -O "${publishedDataDirectory}/wt_120h_H3K27ac_HiChIP_10kb_Raw_mm10.cool.gz" -nc
gunzip -k ${publishedDataDirectory}/wt_120h_H3K27ac_HiChIP_10kb_Raw_mm10.cool.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226246&format=file&file=GSM6226246%5Fwt%5F120h%5FH3K27ac%5Freptc%5FNormalized%2Ebigwig" -O "${publishedDataDirectory}/wt_120h_H3K27ac_reptc_Normalized.bigwig" -nc

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226247&format=file&file=GSM6226247%5Fwt%5F132h%5FH3K27ac%5Freptc%5FNormalized%2Ebigwig" -O "${publishedDataDirectory}/wt_132h_H3K27ac_reptc_Normalized.bigwig" -nc

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226271&format=file&file=GSM6226271%5Fwt%5F132h%5FNIPBL%5FNormalized%2Ebigwig" -O "${publishedDataDirectory}/wt_132h_NIPBL_Normalized.bigwig" -nc

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226294&format=file&file=GSM6226294%5Fwt%5F132h%5FRAD21%5Frep1%5FNormalized%2Ebigwig" -O "${publishedDataDirectory}/wt_132h_RAD21_rep1_Normalized.bigwig" -nc

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6226240&format=file&file=GSM6226240%5Fwt%5F168h%5FCTCF%2Ebigwig" -O "${publishedDataDirectory}/wt_168h_CTCF.bigwig" -nc

# Figure 5A
ini_file=Fig5A.ini
echo "[scalebar]
file_type = scalebar
height = 0.3
fontsize = 6

[wt_132h_H3K27ac]
file = ${publishedDataDirectory}/wt_132h_H3K27ac_reptc_Normalized.bigwig
title = H3K27ac (132h)
height = 0.7
color = #55cf21
min_value = 0
max_value = 350
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[NIPBL]
file = ${publishedDataDirectory}/wt_132h_NIPBL_Normalized.bigwig
title = NIPBL (132h)
height = 0.7
color = #204D48
min_value = 0
max_value = 300
number_of_bins = 2000
show_data_range = true
file_type = bigwig

[RAD21]
file = ${publishedDataDirectory}/wt_132h_RAD21_rep1_Normalized.bigwig
title = RAD21 (132h)
height = 0.7
color = #d167a8
min_value = 0
max_value = 200
number_of_bins = 2000
show_data_range = true
file_type = bigwig

[wt_168h_CTCF]
file = ${publishedDataDirectory}/wt_168h_CTCF.bigwig
title = CTCF (168h)
height = 0.7
color = #fc8403
min_value = 0
max_value = 120
number_of_bins = 2000
show_data_range = true
file_type = bigwig

[CTCF label]
file = ${pathWithAnnotations}/outputs/CTCF_inHox_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.01
height = 0.3

[genes label]
file = ${pathWithAnnotations}/outputs/selected_left_5A.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.3

[nc]
file = ${pathWithAnnotations}/outputs/selected_nc.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false" > ${ini_file}
pyGenomeTracks --tracks ${ini_file} --region chr11:96183548-96376004 --outFileName ${ini_file/.ini/.pdf} --dpi 250 --plotWidth 10.1 --fontSize 6 --trackLabelFraction 0.2

# Figure 5B
ini_file=Fig5B.ini
echo "[scalebar]
file_type = scalebar
height = 0.3
fontsize = 6

[wt_120h_HiChIP]
file = ${publishedDataDirectory}/wt_120h_H3K27ac_HiChIP_10kb_Raw_mm10.cool
title = 120h HiChIP H3K27ac
colormap = hot
depth = 1000000
show_masked_bins = false
min_value = 0
max_value = 50

[spacer]
height = 0.1

[wt_120h_H3K27ac]
file = ${publishedDataDirectory}/wt_120h_H3K27ac_reptc_Normalized.bigwig
title = H3K27ac (120h)
height = 0.7
color = #55cf21
min_value = 0
max_value = 350
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.3

[nc]
file = ${pathWithAnnotations}/outputs/selected_nc.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[labels]
file = ${pathWithAnnotations}/outputs/selected_left_grouped.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

" > ${ini_file}
pyGenomeTracks --tracks ${ini_file} --region chr11:95,634,641-97,036,681 --outFileName ${ini_file/.ini/.pdf} --dpi 250 --plotWidth 10.1 --fontSize 6 --trackLabelFraction 0.3


# Figure 6:
ini_file=Fig6A.ini
echo "[CTCF label]
file = ${pathWithAnnotations}/outputs/CTCF_inHox_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.01
height = 0.2

[genes label]
file = ${pathWithAnnotations}/outputs/selected_left_5A.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2
" > ${ini_file}
for f in ${pathWithAnnotations}/outputs/empty.bed ${pathWithAnnotations}/gastruloid_deletion_HoxB.bed ${pathWithAnnotations}/gastruloid_deletion.bed; do
    echo "[$(basename $f)]
file = $f
display = deletions
color = red
labels = false
height = 0.5
" >> ${ini_file}
    if [[ $(basename $f) != "gastruloid_deletion_HoxB.bed" ]]; then
        echo "[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
overlay_previous = share-y
color = black
display = collapsed
color_utr = black
border_color = none
labels = false

[miR]
file = ${pathWithAnnotations}/outputs/mir_extended.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false
" >> ${ini_file}
    else
        echo "[spacer]
height = 0.1
" >> ${ini_file}
    fi
    if [[ $(basename $f) == "empty.bed" ]]; then
        echo "[gm]
file = ${pathWithAnnotations}/outputs/gm53.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[CTCF cassette]
file = ${pathWithAnnotations}/CTCF_cassette.bed
display = collapsed
overlay_previous = share-y
color = bed_rgb
color_utr = bed_rgb
border_color = none
labels = false

[scalebar]
file_type = scalebar
scalebar_start_position = 96197447
# Should be: 96271457
scalebar_end_position = 96271447
file_type = scalebar
where = right
height = 0.2
fontsize = 6
" >> ${ini_file}
    fi
done
genome=deli9-13
echo "[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.01
height = 0.2

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.35

[nc]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_nc.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[scalebar]
file_type = scalebar
scalebar_start_position = 96197447
# Should be 96204046
scalebar_end_position = 96204047
file_type = scalebar
where = right
height = 0.2
fontsize = 6
" >> ${ini_file}
genome=deli9-13insCBS5-10
echo "[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.01
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.35

[nc]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[CTCF cassette]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette.bed
display = collapsed
overlay_previous = share-y
color = bed_rgb
color_utr = bed_rgb
border_color = none
labels = false

[scalebar]
file_type = scalebar
scalebar_start_position = 96197447
# Should be 96206976
scalebar_end_position = 96206947
file_type = scalebar
where = right
height = 0.2
fontsize = 6" >> ${ini_file}
pgt --region chr11:96190000-96370000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 6.5


# Figure 6AB:
ini_file=Fig6AB.ini
echo "" > ${ini_file}
genome=deli9-13insCBS5-10
echo "[CTCF label]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.01
height = 0.2

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
height = 0.35
color = black
display = collapsed
color_utr = black
border_color = none
labels = false

[nc]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[CTCF cassette]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette.bed
display = collapsed
overlay_previous = share-y
color = bed_rgb
color_utr = bed_rgb
border_color = none
labels = false

[scalebar]
file_type = scalebar
scalebar_start_position = 96197447
# Should be 96206976
scalebar_end_position = 96206947
file_type = scalebar
where = right
height = 0.2
fontsize = 6
" >> ${ini_file}
pgt --region chr11:96192000-96212200 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 6.5

# Figure 6B
ini_file=Fig6B.ini
echo "[scalebar]
file_type = scalebar
scalebar_start_position = 96192500
size = 5000
where = right
height = 0.2
fontsize = 6
" > ${ini_file}
time=96h
for sample in delBins5-10_${time} wt_${time}; do
    echo "[$sample]
file = ${pathWithChIPM}/ChIPM_CTCF_${sample}_onmm10_HoxB_deli9-13insCBS5-10_del.bw
color = #fc8403
min_value = 0
max_value = 1.5
height = 1
title = ChIP-M $sample

[$sample Rad21]
file = ${pathWithChIPM}/ChIPM_RAD21_${sample}_onmm10_HoxB_deli9-13insCBS5-10_del.bw
color = #d167a8
min_value = 0
max_value = 2
height = 1
title = Rad21 ChIP-M $sample
" >> ${ini_file}
done
genome=deli9-13insCBS5-10
echo "[CTCF label]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left_pc.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
overlay_previous = share-y

[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
overlay_previous = share-y
" >> ${ini_file}

pgt --region chr11:96192000-96210000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 6.2 --fontSize 6 --trackLabelFraction 0.3


# Figure 7A:
ini_file=Fig7A.ini
genome=deli9-13insCBS5-10
echo "[scalebar]
file_type = scalebar
scalebar_start_position = 96192500
size = 5000
where = right
height = 0.2
fontsize = 6
" > ${ini_file}
for clone in delBins5-10 delBins5-10del7-10; do
    chipmfile=${pathWithChIPM}/ChIPM_CTCF_${clone}_96h_onmm10_HoxB_deli9-13insCBS5-10_del.bw
    if [ "$clone" = "delBins5-10" ]; then
        name="full cassette"
    elif [ "$clone" = "delBins5-10del7-10" ]; then
        name="Del(CBS7-10)"
    fi
    echo "[chipm]
file = $chipmfile
color = #fc8403
min_value = 0
height = 1
title = $name 96h ChIPM
" >> $ini_file
    if [ "$clone" != "delBins5-10" ]; then
        echo "[del $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del.bed
display = deletions
color = red
labels = false
height = 0.35

[CTCF left $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del_CTCF_left.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
overlay_previous = share-y
" >> $ini_file
    else
        echo "[empty]
file = ${pathWithAnnotations}/outputs/empty.bed
display = deletions
color = red
labels = false
height = 0.35

[CTCF $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
overlay_previous = share-y
" >> $ini_file
    fi
done
for clone in delBins5-10 delBins5-10del8-10; do
    cnr=${pathWithCR}/CUTnRUN_CTCF_${clone}_mESC_onmm10_HoxB_deli9-13insCBS5-10_del.bw
    if [ "$clone" = "delBins5-10" ]; then
        name="full cassette"
    elif [ "$clone" = "delBins5-10del8-10" ]; then
        name="Del(CBS8-10)"
    fi
    echo "[cut and run]
file = $cnr
color = #fc8403
min_value = 0
height = 1
title = $name mES CUT&RUN
" >> $ini_file
    if [ "$clone" != "delBins5-10" ]; then
        echo "[del $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del.bed
display = deletions
color = red
labels = false
height = 0.35

[CTCF left $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del_CTCF_left.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
overlay_previous = share-y
" >> $ini_file
    else
        echo "[empty]
file = ${pathWithAnnotations}/outputs/empty.bed
display = deletions
color = red
labels = false
height = 0.35

[CTCF $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
overlay_previous = share-y
" >> $ini_file
    fi
done
echo "[CTCF label]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.25

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left_pc.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
overlay_previous = share-y

[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.015
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
height = 0.35
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
overlay_previous = share-y
" >> ${ini_file}

pgt --region chr11:96192000-96210000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 6.3 --fontSize 6 --trackLabelFraction 0.3

# Figure 7C:
ini_file=Fig7C.ini
genome=deli9-13insCBS5-10
echo "[empty]
file = ${pathWithAnnotations}/outputs/empty.bed
display = deletions
color = red
labels = false
height = 0.5

[CTCF left $genome]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.1
overlay_previous = share-y
" > ${ini_file}
for clone in delBins5-10del8-10 delBins5-10del7-10; do
    echo "[del $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del.bed
display = deletions
color = red
labels = false
height = 0.5

[CTCF left $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del_CTCF_left.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5
arrowhead_fraction = 0.1
overlay_previous = share-y
" >> $ini_file
done

pgt --region chr11:96198000-96202000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 1.1


# Figure S1B:
ini_file=FigS1B.ini
echo "[scalebar]
file_type = scalebar
fontsize = 6
height = 0.2

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2

[genes label]
file = ${pathWithAnnotations}/outputs/selected_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.3

[miR]
file = ${pathWithAnnotations}/outputs/selected_nc.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[spacer]
height = 0.2

[maf]
file = ${publishedDataDirectory}/chr11.maf
file_type = maf
reference = mm10
color_identical = black
color_mismatch = grey
color_gap = lightgrey
species_order = rn5 hg19 otoGar3 felCat5 sarHar1 galGal4 chrPic1 latCha1 danRer7 ornAna1
species_labels = Rat Human Bushbaby Cat Tasmanian_devil Chicken Painted_turtle Coelacanth Zebrafish Platypus
species_order_only = true
height = 3
rasterize = true

[spacer]
height = 0.2

[repeatmask]
file = ${publishedDataDirectory}/mm10_rmsk.bed.gz
title = Repeat Mask
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.2
arrowhead_fraction = 0
" > ${ini_file}
# Maf is 6757 iterations about 5 minutes
pgt --region chr11:96183548-96376004 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --dpi 250 --plotWidth 14.5 --fontSize 6

# Figure S1C:
ini_file=FigS1C.ini
species_order=$(cut -f 1 -d " " ${pathWithAnnotations}/outputs/mus_gtf/Hoxb13_order.txt | grep -vw "Mus_musculus" | tr "\n" " " | sed "s/M/m/g")
echo "[x-axis]
where = top
fontsize = 6
height = 1

[scalebar]
file_type = scalebar
fontsize = 6
height = 0.2

[genes label]
file = ${pathWithAnnotations}/outputs/mus_gtf/21_murinae.epo.11_3_HoxB_mus_musculus_selected_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[pc]
file = ${pathWithAnnotations}/outputs/mus_gtf/21_murinae.epo.11_3_HoxB_mus_musculus_pc.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.3

[mus_musculus npc]
file = ${pathWithAnnotations}/outputs/mus_gtf/21_murinae.epo.11_3_HoxB_mus_musculus_npc_selected.bed4
overlay_previous = share-y
color = grey
border_color = white
display = collapsed
labels = false

[maf]
file = ${pathWithAnnotations}/outputs/mus_gtf/21_murinae.epo.11_3_HoxB_chunked.maf
file_type = maf
reference = mus_musculus
height = 3.6
species_order = $species_order
species_order_only = true
rasterize = true
" > ${ini_file}
# Maf is 193 iterations about 6 minutes
pgt --region chr11:96074374-96266830 --tracks ${ini_file} -o ${ini_file/.ini/.pdf}  --dpi 250 --plotWidth 13.5 --fontSize 6

# Figure S1D:
ini_file=FigS1D.ini
name=Mus_musculus
echo "[scalebar]
file_type = scalebar
size = 50000
fontsize = 6
height = 0.2

[labels]
file = ${pathWithAnnotations}/outputs/mus_gtf/Mus_musculus_selected_left_Hoxb9_centered.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[$name pc]
file = ${pathWithAnnotations}/outputs/mus_gtf/${name}_Hoxb9_centered_pc.bed4
title = ${name/Mus/mus}
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.2

[$name npc]
file = ${pathWithAnnotations}/outputs/mus_gtf/${name}_Hoxb9_centered_npc_selected.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
labels = false

[genes label]
file = ${pathWithAnnotations}/outputs/mus_gtf/selected_left_Hoxb9_centered.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2
" > $ini_file
for name in $(cut -f 1 -d " " ${pathWithAnnotations}/outputs/mus_gtf/Hoxb13_order.txt | grep -vw "Mus_musculus"); do
    echo "[$name pc]
file = ${pathWithAnnotations}/outputs/mus_gtf/${name}_Hoxb9_centered_pc.bed4
title = ${name/Mus/mus}
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
height = 0.2

[$name npc]
file = ${pathWithAnnotations}/outputs/mus_gtf/${name}_Hoxb9_centered_npc_selected.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
labels = false
" >> $ini_file
done
pgt --region chr11:4890000-5110000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 13.5 --fontSize 5

# Figure S8A:
ini_file=FigS8A.ini
genome=mm10
echo "[scalebar]
file_type = scalebar
size = 50000
fontsize = 6
height = 0.2

[nCATS coverage]
file = ${pathWithnCATS}/nCATS_delBins5-10_on${genome}.bw
height = 1
title = coverage
min_value = 0
max_value = 60

[CRISPR]
file = ${pathWithAnnotations}/CRISPR_nCATS.bed
labels = false
display = collapsed
color = black
title = CRISPR
height = 0.2

[spacer]
height = 0.1

[selected cassette]
file = ${gitHubDirectory}/nCATS/output/${genome}/delBins5-10_cassette.bed
labels = false
title = individual reads

[selected HoxBdel]
file = ${gitHubDirectory}/nCATS/output/${genome}/delBins5-10_HoxBdel.bed
labels = false
height = 0.25

[del HoxB]
file = ${pathWithAnnotations}/gastruloid_deletion_HoxB.bed
display = deletions
color = red
labels = false
height = 0.5

[CTCF label]
file = ${pathWithAnnotations}/outputs/CTCF_inHox_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.3

[CTCF]
file = ${pathWithAnnotations}/CTCF_inHox_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.3
arrowhead_fraction = 0.01

[genes label]
file = ${pathWithAnnotations}/outputs/selected_left_5A.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/outputs/protein_coding_around_HoxBD.bed4
color = black
display = collapsed
color_utr = black
border_color = none
labels = false

[nc]
file = ${pathWithAnnotations}/outputs/selected_nc.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[CTCF cassette]
file = ${pathWithAnnotations}/CTCF_cassette.bed
display = collapsed
overlay_previous = share-y
color = bed_rgb
color_utr = bed_rgb
border_color = none
labels = false
" > ${ini_file}

pgt --region chr11:96180000-96380000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 8.5 --fontSize 6

# Figure S8B:
ini_file=FigS8B.ini
genome=mm10_HoxB_deli9-13insCBS5-10_del
echo "[scalebar]
file_type = scalebar
fontsize = 6
height = 0.2

[nCATS coverage]
file = ${pathWithnCATS}/nCATS_delBins5-10_on${genome}.bw
height = 1
title = coverage
min_value = 0
max_value = 60

[CRISPR]
file = ${pathWithAnnotations}/deli9-13insCBS5-10/deli9-13insCBS5-10_CRISPR_nCATS.bed
labels = false
display = collapsed
color = black
title = CRISPR
height = 0.2

[spacer]
height = 0.1

[selected cassette]
file = ${gitHubDirectory}/nCATS/output/${genome}/delBins5-10_cassette.bed
labels = false
title = individual reads
" > ${ini_file}
genome=deli9-13insCBS5-10
echo "[CTCF label]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.3

[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.3
arrowhead_fraction = 0.01

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
height = 0.5
color = black
display = collapsed
color_utr = black
border_color = none
labels = false

[nc]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_nc_noGm53.bed4
overlay_previous = share-y
color = grey
border_color = none
display = collapsed
color_utr = black
labels = false

[CTCF cassette]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_cassette.bed
display = collapsed
overlay_previous = share-y
color = bed_rgb
color_utr = bed_rgb
border_color = none
labels = false
" >> ${ini_file}
pgt --region chr11:96192000-96213000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 8.5 --fontSize 6


# Figure S9A:
ini_file=FigS9A.ini
genome=deli9-13insCBS5-10
echo "[scalebar]
file_type = scalebar
scalebar_start_position = 96192500
size = 5000
fontsize = 6
height = 0.2
where = right
" > ${ini_file}
for clone in delBins5-10 delBinv5-10 delBinv5-10del7-10; do
    cnr=${pathWithCR}/CUTnRUN_CTCF_${clone}_mESC_onmm10_HoxB_deli9-13insCBS5-10_del.bw
    if [ "$clone" = "delBins5-10" ]; then
        name="full cassette"
    elif [ "$clone" = "delBinv5-10" ]; then
        name="inverted cassette"
    elif [ "$clone" = "delBinv5-10del7-10" ]; then
        name="inverted Del(CBS7-10)"
    fi
    echo "[cut and run]
file = $cnr
color = #fc8403
min_value = 0
height = 1.5
title = $name mES CUT&RUN
" >> $ini_file
    if [ "$clone" == "delBinv5-10" ]; then
        echo "[inv $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_inv.bed
display = inversions
color = black
labels = false
height = 1.5
arrowhead_fraction = 0.007

[empty]
file = ${pathWithAnnotations}/outputs/empty.bed
display = deletions
color = red
labels = false
overlay_previous = share-y

[CTCF deli9-13insCBS5-10]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.02
overlay_previous = share-y
" >> $ini_file
    elif [ "$clone" == "delBins5-10" ]; then
        echo "[empty]
file = ${pathWithAnnotations}/outputs/empty.bed
display = deletions
color = red
labels = false
height = 0.75

[CTCF $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5
arrowhead_fraction = 0.02
overlay_previous = share-y
" >> $ini_file
    else
        echo "[inv $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_inv5-10_inv.bed
display = inversions
color = black
labels = false
height = 1.5
arrowhead_fraction = 0.007

[del $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del.bed
display = deletions
color = red
labels = false
overlay_previous = share-y

[CTCF left $clone]
file = ${pathWithAnnotations}/${genome}/${genome}_${clone/delB/}_del_CTCF_left.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.5
arrowhead_fraction = 0.02
overlay_previous = share-y
" >> ${ini_file}
    fi
done
echo "[CTCF label]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_label_zoom.bed
color = none
border_color = none
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.3

[CTCF]
file = ${pathWithAnnotations}/${genome}/${genome}_CTCF_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
arrowhead_fraction = 0.02
height = 0.5

[genes label]
file = ${pathWithAnnotations}/${genome}/${genome}_selected_left_pc.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
arrowhead_fraction = 0
fontsize = 6
height = 0.2

[genes]
file = ${pathWithAnnotations}/${genome}/${genome}_protein_coding_around_HoxBD.bed4
height = 0.5
color = black
display = collapsed
color_utr = black
border_color = none
labels = false
" >> ${ini_file}
pgt --region chr11:96192000-96210000 --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --plotWidth 8.5 --fontSize 6
