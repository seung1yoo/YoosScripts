#!/bin/bash

vcftools="/BiO/BioTools/vcftools/vcftools_0.1.11/bin/vcftools"
significantSNP_finder="./TBD171076_var-filter.py"
vcf_ifn="multisample.vcf"
vcf_ofp="multisample"
chrs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")
minGQ=15
minDP=10

vcf_spliteds=( )
for chr in ${chrs[@]};do
    #echo ${chr}
    echo "${vcftools} --vcf ${vcf_ifn} --out ${vcf_ofp}.${chr} --chr ${chr} --recode --recode-INFO-all"
    #echo "${vcftools} --vcf ${vcf_ifn} --out ${vcf_ofp}.${chr} --chr ${chr} --recode --recode-INFO-all" | /bin/bash
    vcf_spliteds=( "${vcf_spliteds[@]}" "${vcf_ofp}.${chr}.recode.vcf" )
done;

for vcf in ${vcf_spliteds[@]};do
    #echo ${vcf}
    arrays=($(echo $vcf | tr '.' ' ')) 
    outprefix="${vcf_ofp}.sig_var.${arrays[1]}"
    echo "/usr/bin/python ${significantSNP_finder} --vcf ${vcf} --vcftools ${vcftools} --plink plink --outprefix ${outprefix}"
    echo "/usr/bin/python ${significantSNP_finder} --vcf ${vcf} --vcftools ${vcftools} --plink plink --outprefix ${outprefix}" | /bin/bash
done;


