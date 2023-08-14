#!/bin/bash
# Input: allQCBAFilterVariantMap.txt.gz
# Output: pLOF.txt, npLOF.txt, other.txt, Missense.txt, Synonymous.txt, AnnoMerged_prefiltered.txt, genotypes_annotations.csv.gz filtered_genos.csv

regenie_folder=$1

top=$regenie_folder
dir="$top/Anno"

# Merges the separate var-geno files into a single anno file used by regenie for burden testing
## used as --anno-file in step 2
## you need to adjust mask.mask file to include the specific variant effects you want

mkdir $dir
touch ${dir}/AnnoMerged_prefiltered.txt

echo "Beginning Variant Annotation file creation"

# This command searches for lines in the file ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz (which is uncompressed using zcat) where the fourth column contains the pattern "synonymous_variant". It then extracts the first, fifth, and fourth columns from those lines and saves the results in the "temp" file.
awk '{if($4 ~ /synonymous_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) > temp
awk '{if($4 ~ /start_retained/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /stop_retained_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
# In summary, the command reads the file "temp" and replaces the contents of the third column in each line with the value "synonymous". It then saves the modified output, containing columns 1, 2, and the modified column 3, to a file named "Synonymous.txt" in the specified directory.
awk '{ $3="synonymous" ; print $1,$2,$3}' temp > ${dir}/Synonymous.txt
rm temp

# The command searches for lines in the uncompressed file ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz where the fourth column contains the pattern "missense_variant". It then extracts the first, fifth, and fourth columns from those lines and saves the results in the "temp" file.
awk '{if($4 ~ /missense_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) > temp
awk '{if($4 ~ /missense_variant+splice_region_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
# The command reads the file "temp" and replaces the contents of the third column in each line with the value "missense". It then saves the modified output, containing columns 1, 2, and the modified column 3, to a file named "Missense.txt" in the specified directory.
awk '{ $3="missense" ; print $1,$2,$3}' temp > ${dir}/Missense.txt
rm temp

awk '{if($4 ~ /frameshift_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) > temp
awk '{if($4 ~ /stop_gained/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /splice_donor_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /splice_acceptor_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp

awk '{ $3="pLOF" ; print $1,$2,$3}' temp > ${dir}/pLOF.txt
rm temp


# The command searches for lines in the uncompressed file ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz where the fourth column contains the pattern "conservative_inframe_insertion". It then extracts the first, fifth, and fourth columns from those lines and appends the results to the "temp" file.
awk '{if($4 ~ /conservative_inframe_insertion/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /conservative_inframe_deletion/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /rare_amino_acid_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /coding_sequence_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{ $3="other" ; print $1,$2,$3}' temp > ${dir}/other.txt
rm temp

awk '{if($4 ~ /start_loss/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) > temp
awk '{if($4 ~ /disruptive_inframe_deletion/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /disruptive_inframe_insertion/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp
awk '{if($4 ~ /exon_loss_variant/) {print $1,$5,$4} }' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) >>temp

awk '{ $3="npLOF" ; print $1,$2,$3}' temp > ${dir}/npLOF.txt
rm temp

cat ${dir}/pLOF.txt ${dir}/npLOF.txt ${dir}/other.txt ${dir}/Missense.txt ${dir}/Synonymous.txt > ${dir}/AnnoMerged_prefiltered.txt

# awk 1 ${dir}/Missense.txt ${dir}/pLOF.txt > ${dir}/temp
# awk 1 ${dir}/temp ${dir}/npLOF.txt > ${dir}/temp
# awk 1 ${dir}/temp ${dir}/other.txt > ${dir}/temp
# awk 1 ${dir}/temp ${dir}/Synonymous.txt > ${dir}/AnnoMerged_prefiltered.txt
# The command removes duplicate lines from the file "AnnoMerged.txt" based on the values in the first column. It saves the output to a temporary file named "tmp" and then replaces the original file with the updated contents.
awk '!visited[$1]++' ${dir}/AnnoMerged_prefiltered.txt > tmp && mv -f tmp ${dir}/AnnoMerged_prefiltered.txt
rm ${dir}/temp

#                      V1     V2         V3
#      1:    1-861349-C-T SAMD11   missense
#      2:    1-865595-A-G SAMD11   missense
#      3:    1-865602-G-A SAMD11   missense
#      4:    1-865628-G-A SAMD11   missense
#      5:    1-865633-C-A SAMD11   missense
#     ---                                  
# 640962: X-155127859-A-G  VAMP7 synonymous
# 640963: X-155127898-C-G  VAMP7 synonymous
# 640964: X-155234224-C-T   IL9R synonymous
# 640965: X-155235041-C-G   IL9R synonymous
# 640966: X-155235086-T-C   IL9R synonymous

echo "AnnoMerged.txt created"

##uses the all varant map created from Gundula's R script and the gene_boundary file. Output is listFile.txt
## This gets used as --set-list in step2

gene_boundaries="/nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.ensembl87.txt"
# The command reads the contents of the file specified by $gene_boundaries. It modifies the third column by replacing ( with .., splits the modified third column into an array, and extracts specific values from the input file. The extracted values, along with the values from the first and second columns, are saved to a new file named "listFile.txt" in the specified directory.
awk -F " " '{ gsub("\\(","..",$3);  split($3,arr1,"\\.."); print $1,$2,arr1[2]}' $gene_boundaries > ${dir}/listFile.txt

awk '{ print $5,$1}' <(zcat ${top}/allQCBAFilter/allQCBAFilterVariantMap.txt.gz) > ${dir}/temp.txt

awk '{
      if(NR!=1){a[$1]=$2","a[$1]}
      else print $0}
    END{
      n = asorti(a, b);
      for (n in b) {
      print b[n],a[b[n]]
      }
    }' ${dir}/temp.txt > tmp.txt && mv -f tmp.txt ${dir}/temp.txt
sed 's/,$//' ${dir}/temp.txt > tmp.txt && mv -f tmp.txt ${dir}/temp.txt
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' ${dir}/temp.txt ${dir}/listFile.txt > tmp.txt && mv -f tmp.txt ${dir}/listFile.txt

rm ${dir}/temp.txt

awk -F " "  '$4!="" {print $0}' ${dir}/listFile.txt > tmp.txt && mv -f tmp.txt ${dir}/listFile.txt
awk -F " "  '$2!="Y" {print $0}' ${dir}/listFile.txt > tmp.txt && mv -f tmp.txt ${dir}/listFile.txt

echo "listFile.txt created"

# The command writes the formatted text to a file named "mask.mask" in the specified directory. The contents of the file will include lines of text with specific tab-separated values, representing different categories and associated values.
printf "Synonymous\tsynonymous\nLOF\tpLOF,npLOF\npLOF\tpLOF\nMissense\tMPC0,MPC2,MPC3\nAll\tpLOF,npLOF,MPC0,MPC2,MPC3,other\nMPC0\tMPC0\nMPC2\tMPC2,MPC3\nDeleterious\tMPC2,MPC3,pLOF" > ${dir}/mask.mask
# printf "Synonymous\tsynonymous\nLOF\tLOF\nMissense\tmissense,MTR50,LIMBR50,MultiPath\nAll\tLOF,missense,LIMBR50,MTR50,MultiPath\nMTR50\tMTR50\nMultiPath\tMultiPath\nLIMBR50\tLIMBR50" > ${dir}/mask.mask

echo "mask.mask created, now sorting genofile. This takes a few hours"

# The command creates a new CSV file named "filtered_genos.csv" in the ${top}/allQCBAFilter/genotypes/ directory. The file contains the header from any CSV file in the directory, followed by the unique and sorted lines from the remaining CSV files, excluding their headers.
(head -n 1 ${top}/allQCBAFilter/genotypes/*genotypes_annotations.csv && ( tail -n +2 ${top}/allQCBAFilter/genotypes/*genotypes_annotations.csv | sort -k1,1 -t "," -u )) > ${top}/allQCBAFilter/genotypes/filtered_genos.csv

#           Variant ID Variant Type Ref Allele Alt Allele   Rs Number   Impact                                 Effect
     # 1: 10-100003875-G-A          snv          G          A rs148923672 MODERATE                       missense_variant
     # 2: 10-100010822-C-A          snv          C          A        <NA> MODERATE missense_variant+splice_region_variant
     # 3: 10-100010922-G-T          snv          G          T rs192574171 MODERATE                       missense_variant
     # 4: 10-100012115-G-A          snv          G          A        <NA> MODERATE                       missense_variant
     # 5: 10-100012198-G-A          snv          G          A rs372115263      LOW                     synonymous_variant
     #    Canonical Transcript Effect Gene Name Transcript Stable Id Has CCDS Transcript
     # 1:                                                              missense_variant 'R3HCC1L'      ENST00000298999                TRUE
     # 2: missense_variant+splice_region_variant|missense_variant|splice_region_variant   'LOXL4'      ENST00000260702                TRUE
     # 3:                                                              missense_variant   'LOXL4'      ENST00000260702                TRUE
     # 4:                                                              missense_variant   'LOXL4'      ENST00000260702                TRUE
     # # 5:                                                            synonymous_variant   'LOXL4'      ENST00000260702                TRUE           
     # HGVS_c      HGVS_p Polyphen Humdiv Score Polyphen Humdiv Prediction Polyphen Humvar Score Polyphen Humvar Prediction
     # 1: c.2297G>A p.Arg766Gln                 1.000                   probably                 0.998                   probably
     # 2: c.2200G>T p.Gly734Trp                 0.815                   possibly                 0.656                   possibly
     # 3: c.2100C>A p.Asn700Lys                 0.997                   probably                 0.973                   probably
     # 4: c.1946C>T p.Pro649Leu                 0.727                   possibly                 0.734                   possibly
     # 5: c.1863C>T p.Thr621Thr                    NA                       <NA>                    NA                       <NA>
# and more
echo "sort complete, filtered_genos.csv created, gzipping genotypes.csv"

annotation_file="${top}/allQCBAFilter/genotypes/*genotypes_annotations.csv"
gzip -d $annotation_file
# rm ${top}/allQCBAFilter/genotypes/*genotypes_annotations.csv

echo "gzipping complete"

exit
