#!/bin/bash

vcftools="vcftools"
significantSNP_finder="./TBD200193_var-filter.py"
vcf_ifn="multisample.haedeul.vcf"
vcf_ofp="multisample.haedeul"
chrs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12")
minGQ=15
minDP=10
sample1="gopum156"
sample2="gangwon4ho"

vcf_spliteds=( )
for chr in ${chrs[@]};do
    echo ${chr}
    echo "${vcftools} --vcf ${vcf_ifn} --out ${vcf_ofp}.${chr} --chr ${chr} --recode --recode-INFO-all"
    echo "${vcftools} --vcf ${vcf_ifn} --out ${vcf_ofp}.${chr} --chr ${chr} --recode --recode-INFO-all" | /bin/bash
    vcf_spliteds=( "${vcf_spliteds[@]}" "${vcf_ofp}.${chr}.recode.vcf" )
done;

for vcf in ${vcf_spliteds[@]};do
    echo ${vcf}
    arrays=($(echo $vcf | tr '.' ' ')) 
    outprefix="${vcf_ofp}.sig_var.${arrays[2]}"
    echo "/usr/bin/python2 ${significantSNP_finder} --vcf ${vcf} --vcftools ${vcftools} --plink plink --outprefix ${outprefix} --sample1 ${sample1} --sample2 ${sample2}"
    echo "/usr/bin/python2 ${significantSNP_finder} --vcf ${vcf} --vcftools ${vcftools} --plink plink --outprefix ${outprefix} --sample1 ${sample1} --sample2 ${sample2}" | /bin/bash
done;


