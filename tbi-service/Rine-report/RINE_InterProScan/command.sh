
head -n1 ../genes.de.xls > genes.de.headers

head -n1 ../genes.de.xls | cut -f 2,17,18,19,20 | tr '\t' '\n' > genes.de.headers.Ga_8-DG_4
head -n1 ../genes.de.xls | cut -f 2,22,23,24,25 | tr '\t' '\n' > genes.de.headers.Ga_8-O_10
head -n1 ../genes.de.xls | cut -f 2,27,28,29,30 | tr '\t' '\n' > genes.de.headers.Ga_8-Ma_5
head -n1 ../genes.de.xls | cut -f 2,32,33,34,35 | tr '\t' '\n' > genes.de.headers.Ga_8-Gi_8
head -n1 ../genes.de.xls | cut -f 2,37,38,39,40 | tr '\t' '\n' > genes.de.headers.Ga_8-Heart_5
head -n1 ../genes.de.xls | cut -f 2,42,43,44,45 | tr '\t' '\n' > genes.de.headers.Ga_8-Muscle_9
head -n1 ../genes.de.xls | cut -f 2,47,48,49,50 | tr '\t' '\n' > genes.de.headers.Ga_8-PPG_1
head -n1 ../genes.de.xls | cut -f 2,52,53,54,55 | tr '\t' '\n' > genes.de.headers.Ga_8-BG_1
head -n1 ../genes.de.xls | cut -f 2,57,58,59,60 | tr '\t' '\n' > genes.de.headers.Ga_8-T_2


deg_tags=("Ga_8-DG_4" "Ga_8-O_10" "Ga_8-Ma_5" "Ga_8-Gi_8" "Ga_8-Heart_5" "Ga_8-Muscle_9" "Ga_8-PPG_1" "Ga_8-BG_1" "Ga_8-T_2")

for deg_tag in ${deg_tags[@]}
do
    /usr/bin/python ~/YoosScripts/tbi-service/Rine-report/genesColumnSelector.py -ig ../genes.de.xls -hl genes.de.headers.$deg_tag -hi GeneId -ofn GenLst.$deg_tag.tmp
    head -n1 GenLst.$deg_tag.tmp > GenLst.$deg_tag.txt
    cat GenLst.$deg_tag.tmp | awk '$2 == "Y"' >> GenLst.$deg_tag.txt
    rm GenLst.$deg_tag.tmp
done

/usr/bin/python bin/mkIPRS4DEG.py ../InterproScanAnnotation.xls

rm fileLst.txt
touch fileLst.txt
for deg_tag in ${deg_tags[@]}
do
    echo -e "InterproScanAnnotation.$deg_tag.xls\t$deg_tag" >> fileLst.txt
done

for deg_tag in ${deg_tags[@]}
do
    /usr/bin/python bin/InterProScanGeneToGOTable.py InterproScanAnnotation.$deg_tag.xls Gene2GO.$deg_tag.xls
done

/usr/bin/python bin/InterProScanGOlevelToTable.py fileLst.txt GO2GENE

/usr/bin/python bin/InterProScanKEGGPathwayToTable.py fileLst.txt KEGG4DEG
