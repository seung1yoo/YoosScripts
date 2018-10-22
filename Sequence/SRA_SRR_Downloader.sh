#!/bin/bash
srapath="./sratoolkit.2.9.2-ubuntu64/bin/srapath"
fastqdump="./sratoolkit.2.9.2-ubuntu64/bin/fastq-dump"

srrs=("SRR7132515" "SRR7132516" "SRR7132517" "SRR7132518" "SRR7132519" "SRR7132520" "SRR7132521" "SRR7132522")

for srr in ${srrs[@]}
do
    srr_path=$($srapath $srr)
    wget_cmd="wget $srr_path"
    echo $wget_cmd 
    echo $wget_cmd | /bin/bash
    #
    fqdump_cmd="$fastqdump --split-files ./${srr}"
    echo $fqdump_cmd
    echo $fqdump_cmd | /bin/bash
    #
    gzip_cmd="gzip -c ${srr}_1.fastq > ${srr}_1.fq.gz"
    echo $gzip_cmd
    echo $gzip_cmd | /bin/bash
    #
    gzip_cmd="gzip -c ${srr}_2.fastq > ${srr}_2.fq.gz"
    echo $gzip_cmd
    echo $gzip_cmd | /bin/bash
    #
    rm_cmd="rm -f ${srr}_1.fastq ${srr}_2.fastq"
    echo $rm_cmd
    echo $rm_cmd | /bin/bash
done


