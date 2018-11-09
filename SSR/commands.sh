## TBI SSR Finder v1
## 

unigene="TBI.unigene.fa"
outprefix="TBI.unigene.SSR"

## Manual

## make reference.stat
cmd="clc_sequence_info -r -n $unigene > $unigene.SeqInfo"
echo $cmd
# make reference.stat MANUAL
# #Sequences	TotalLength	MaxLength	MinLength	N50	N90
# 2221	516564219	4955513	10030	649252	115747

## Auto-mation

## make SSR Result
cmd="perl bin/ssr.pl $unigene > $outprefix"
echo $cmd
cmd="/usr/bin/python bin/SSRFilter.py $unigene $unigene.xls"
echo $cmd


## SSR Parsing
cmd="/usr/bin/python bin/SSRParsing.py $outprefix.xls reference.stat $outprefix.Freq.Mer.xls $outprefix.Freq.Motif.xls"
echo $cmd

## SSR Statistics
cmd="/usr/bin/python bin/SSRStatistics.py $outprefix.Freq.Motif.xls $outprefix.Freq.Motif.di-tri.xls"
echo $cmd
