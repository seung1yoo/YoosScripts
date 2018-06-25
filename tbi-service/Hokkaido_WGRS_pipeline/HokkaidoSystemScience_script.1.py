
import os

def programPathMaker():
    #titan
    binDic = {'cutadapt':'/BiO/BioTools/cutadapt/cutadapt-1.1/bin/cutadapt',
              'trimmomatic':'/BiO/BioTools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar',
              'bwa':'/BiO/BioTools/Rine/Tools/bwa/bwa-0.7.10/bwa',
              'samtools':'/BiO/BioTools/Rine/Tools/samtools/samtools-1.2/samtools',
              'gatk':'/BiO/BioTools/GATK/gatk-lite-2.3.0-0/jar/GenomeAnalysisTKLite.jar',
              'picard':'/BiO/BioPeople/siyoo/BioTools/picard/picard-tools-1.133/picard.jar',
              'bcftools':'/BiO/BioPeople/siyoo/BioTools/bcftools/bcftools-1.2/bcftools',
              'ref':'/BiO/BioProjects/TBO180061-Hokkaido-Cow-WGRS-20180418/Ref/Bos_taurus.UMD3.1.dna.toplevel.fa'}
    return binDic

def main(args):
    binDic = programPathMaker()
    #print(args)
    #print(binDic)

    # cutadapt 1.1
    ### ORIGINAL
    cmds = ['{0} --match-read-wildcards -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {1}_1.nonadapt.fq.gz {2}'.format(binDic['cutadapt'], args.sample, args.forward),
            '{0} --match-read-wildcards -O 1 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {1}_2.nonadapt.fq.gz {2}'.format(binDic['cutadapt'], args.sample, args.reverse)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # trimmomatic 0.32
    cmds = ['java -jar {0} PE -threads 10 -phred33 {1}_1.nonadapt.fq.gz {1}_2.nonadapt.fq.gz {1}_1.clean.paired.fq.gz {1}_1.clean.unpaired.fq.gz {1}_2.clean.paired.fq.gz {1}_2.clean.unpaired.fq.gz LEADING:0 TRAILING:0 SLIDINGWINDOW:20:20 MINLEN:70'.format(binDic['trimmomatic'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # BWA 0.7.10
    cmds = ["{0} mem -t 10 -aM -R '@RG\\tID:{1}\\tSM:{1}\\tLB:{1}\\tPL:ILLUMINA' {2} {1}_1.clean.paired.fq.gz {1}_2.clean.paired.fq.gz > {1}.mapped.sam".format(binDic['bwa'], args.sample, binDic['ref'])]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # samtools 1.2
    cmds = ['{0} fixmate -O bam {1}.mapped.sam {1}.fixmate.bam'.format(binDic['samtools'], args.sample),
            '{0} sort -@ 10 -O bam -o {1}.sorted.bam -T {1}.temp {1}.fixmate.bam'.format(binDic['samtools'], args.sample),
            '{0} index {1}.sorted.bam'.format(binDic['samtools'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # Extract Unmap Reads
    cmds = ['{0} view -bh -f 0x4 -f 0x8 {1}.sorted.bam | {0} bam2fq -O - > {1}.unmap_pair.fastq'.format(binDic['samtools'], args.sample),
            '''cat %s.unmap_pair.fastq  |  paste - - - - - - - -  |  awk  ' BEGIN{ FS = "\\t"; OFS = "\\n"; }{ if( $1 != gensub(/\/2$/, "/1", "", $5) ){ print $1, $2, $3, $4, $5, $6, $7, $8  |  "gzip  -c  >  %s.unmapped_pair_PE.fq.gz "; }else{ print $1, $2, $3, $4  |  "gzip  -c  >  %s.unmapped_pair_R1.fq.gz ";  print $5, $6, $7, $8  |  "gzip  -c  >  %s.unmapped_pair_R2.fq.gz"; } }' '''% (args.sample, args.sample, args.sample, args.sample),
            '''{0} view  -bh  -f  0x4  -F  0x8  {1}.sorted.bam  |  {0} bam2fq  -O  -  >  {1}.unmap_single.fastq'''.format(binDic['samtools'], args.sample),
            '''cat  %s.unmap_single.fastq  |  paste - - - -  |  awk  ' BEGIN{ FS = "\\t"; OFS = "\\n"; }{ if( $1 ~ /\/[12]$/ ){ rid = "R" substr($1, length($1), 1) }else{ rid = "SR"; }; print $1, $2, $3, $4  |  "gzip  -c  >  %s.unmapped_single_"rid".fq.gz"; }' ''' % (args.sample, args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # GenomeAnalysisTK:Lite-2.3.0-0
    #cmds = ['java -jar {0} -T RealignerTargetCreator --num_threads 10 -R {1} -I {2}.sorted.bam -o {2}.intervals'.format(binDic['gatk'], binDic['ref'], args.sample),
    #        'java -jar {0} -T IndelRealigner -R {1} -I {2}.sorted.bam -targetIntervals {2}.intervals  -rf NotPrimaryAlignment  -o {2}.realigned.bam'.format(binDic['gatk'], binDic['ref'], args.sample)]
    cmds = ['java -jar {0} -T IndelRealigner -R {1} -I {2}.sorted.bam -targetIntervals {2}.intervals  -rf NotPrimaryAlignment  -o {2}.realigned.bam'.format(binDic['gatk'], binDic['ref'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # picard 1.133
    cmds = ['java -jar {0} MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT={1}.realigned.bam  OUTPUT={1}.mrd.bam  M={1}.metrics'.format(binDic['picard'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # samtools 1.2
    cmds = ['{0} index {1}.mrd.bam'.format(binDic['samtools'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # GenomeAnalysisTK:Lite-2.3.0-0
    cmds = ['java -jar {0} -T RealignerTargetCreator --num_threads 10 -R {1} -I {2}.mrd.bam -o {2}.2nd.intervals'.format(binDic['gatk'], binDic['ref'], args.sample),
            'java -jar {0} -T IndelRealigner -R {1} -I {2}.mrd.bam -targetIntervals {2}.2nd.intervals -rf NotPrimaryAlignment -o {2}.final.bam'.format(binDic['gatk'], binDic['ref'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # Mapping QC (samtools 1.2)
    cmds = ['{0} flagstat {1}.final.bam > {1}.flagstat.txt'.format(binDic['samtools'], args.sample)]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    # picard 1.133
    cmds = ['java -jar {0} CollectMultipleMetrics VALIDATION_STRINGENCY=LENIENT INPUT={1}.final.bam OUTPUT={1}.cmm.metrics'.format(binDic['picard'], args.sample),
            'java -jar {0} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT INPUT={1}.final.bam OUTPUT={1}.cwm.metrics REFERENCE_SEQUENCE={2}'.format(binDic['picard'], args.sample, binDic['ref'])]
    for cmd in cmds:
        #print cmd
        #os.system(cmd)
        pass

    print 'DONE'

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--forward', help='forward reads')
    parser.add_argument('-r', '--reverse', help='reverse reads')
    parser.add_argument('-s', '--sample', help='sample name')
    args = parser.parse_args()
    main(args)
