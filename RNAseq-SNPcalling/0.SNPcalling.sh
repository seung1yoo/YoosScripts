#!/bin/bash

### SNP calling From RNA-seq data (tophat bam files)
### JUST SIMPLE !

### variations
samtools_mpileup_path="/BiO/BioTools/Rine/Tools/samtools/samtools-1.2/samtools mpileup"
bcftools_call_path="/BiO/BioPeople/siyoo/BioTools/bcftools/bcftools-1.2/bcftools call"
java_path="/BiO/BioPeople/shinig/BioTools/GATK4_ServerAnalysis/Tools/Java/jdk1.8.0_40/bin/java"
snpeff_path="/BiO/BioPeople/shinig/BioTools/GATK4_ServerAnalysis/Tools/SNPEFF/snpEff.jar"
snpeff_config="/BiO/BioPeople/shinig/BioTools/snpeff/v4.1g/snpEff.config"
snpeff_OnePerLine="/BiO/BioPeople/shinig/BioTools/GATK4_ServerAnalysis/Tools/SNPEFF/scripts/vcfEffOnePerLine.pl"
snpeff_sift="/BiO/BioPeople/shinig/BioTools/GATK4_ServerAnalysis/Tools/SNPEFF/SnpSift.jar"
script_path="1.split_snpeff.xls.py"

project_path="/BiO/BioPeople/siyoo/TBD180045-gntech-Human-RNAref-SNPcalling"
snpeff_v="H_sapiens_ENS_72"
ref="/BiO/BioResources/References/Homo_sapiens/H_sapiens_ENS_72/H_sapiens_ENS_72.chr.fa"
chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
in_bams="bam_tophat/S0001.bam bam_tophat/S0002.bam bam_tophat/S0003.bam"

### make "multisample.vcf"
MPILEUP="$samtools_mpileup_path -d 10000 -L 10000 -I -B -t DP,DPR,DV,DP4,SP -f $ref $in_bams -go multisample.bcf"
echo $MPILEUP
echo $MPILEUP | /bin/bash
CALL="$bcftools_call_path -Avm -O v -f GQ -o multisample.vcf multisample.bcf"
echo $CALL
echo $CALL | /bin/bash

### SNPEFF
SNPEFF="$java_path -Xmx24g -Djava.io.tmpdir=$project_path/tmp/ -jar $snpeff_path -canon -geneId -c $snpeff_config -v $snpeff_v -s $project_path/multisample.snpeff.html -o vcf $project_path/multisample.vcf > $project_path/multisample.snpeff.vcf"
echo ${SNPEFF}
echo ${SNPEFF} | /bin/bash

### SNPSIFT
sample_ID_list=$(cat $project_path/multisample.vcf | grep '#CHROM' | awk -F '\t' '{for (i=10;i<=NF;i=i+1) print $i}' | tr '\n' '\t' | sed 's/\t$//g')
SNPSIFT="cat $project_path/multisample.snpeff.vcf | $snpeff_OnePerLine | $java_path -Xmx24g -Djava.io.tmpdir=$project_path/tmp/ -jar $snpeff_sift extractFields - -e '.' CHROM POS REF ALT QUAL DP GEN[*].GT 'EFF[*].EFFECT' 'EFF[*].IMPACT' 'EFF[*].CODON' 'EFF[*].AA' 'EFF[*].GENE' 'EFF[*].BIOTYPE' 'EFF[*].TRID' | sed \"s/GEN\[\*\].GT/${sample_ID_list}/g\" > $project_path/multisample.snpeff.xls"
echo ${SNPSIFT}
echo ${SNPSIFT} | /bin/bash

### split by CHRM & filter out DP10
SPLIT="python $script_path --chrs $chrs --xls multisample.snpeff.xls --prefix multisample.snpeff --depth 10"
echo ${SPLIT}
echo ${SPLIT} | /bin/bash

