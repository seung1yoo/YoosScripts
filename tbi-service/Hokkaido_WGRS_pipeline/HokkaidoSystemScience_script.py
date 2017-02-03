
import os

def programPathMaker():
    binDic = {'cutadapt':'/BiO/BioPeople/siyoo/TBO160406-Hokkaido-WGRS-20170112/11.Tools/cutadapt-1.1/bin/cutadapt',
              'trimmomatic':'/BiO/BioTools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar',
              'bwa':'/BiO/BioTools/bwa/bwa-0.7.10/bwa',
              'ref':'/BiO/BioPeople/siyoo/TBO160406-Hokkaido-WGRS-20170112/10.Ref/GRCm38_68.fa',
              'samtools':'/BiO/BioTools/samtools/1.2/bin/samtools',
              'gatk':'/BiO/BioPeople/siyoo/TBO160406-Hokkaido-WGRS-20170112/11.Tools/gatk-lite-2.3.0-0/jar/GenomeAnalysisTKLite.jar',
              'picard':'/BiO/BioPeople/siyoo/TBO160406-Hokkaido-WGRS-20170112/11.Tools/picard-tools-1.133/picard.jar',
              'bcftools':'/BiO/BioTools/bcftools/1.2/bcftools'}
    return binDic

def main(args):
    binDic = programPathMaker()
    print(args)
    print(binDic)

    # cutadapt 1.6
    cmds = ['{0} --match-read-wildcards -O 1 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {1}_1.nonadapt.fq.gz {2}'.format(binDic['cutadapt'], args.sample, args.forward),
            '{0} --match-read-wildcards -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {1}_2.nonadapt.fq.gz {2}'.format(binDic['cutadapt'], args.sample, args.reverse)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # trimmomatic 0.32
    cmds = ['java -jar {0} PE -phred33 {1}_1.nonadapt.fq.gz {1}_2.nonadapt.fq.gz {1}_1.clean.paired.fq.gz {1}_1.clean.unpaired.fq.gz {1}_2.clean.paired.fq.gz {1}_2.clean.unpaired.fq.gz '\
            'LEADING:0 TRAILING:0 SLIDINGWINDOW:20:20 MINLEN:70'.format(binDic['trimmomatic'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # BWA 0.7.10
    cmds = ["{0} mem -aM -R '@RG\tID:{1}\tSM:{1}\tPL:ILLUMINA' {2} {1}_1.clean.paired.fq.gz {1}_2.clean.paired.fq.gz | grep -v '@PG' > {1}.mapped.sam".format(binDic['bwa'], args.sample, binDic['ref'])]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # samtools 1.2
    cmds = ['{0} fixmate -O bam {1}.mapped.sam {1}.mapped.bam'.format(binDic['samtools'], args.sample),
            '{0} sort -O bam -o {1}.mapped.sorted.bam -T {1}.mapped.sorted.temp {1}.mapped.bam'.format(binDic['samtools'], args.sample),
            '{0} index {1}.mapped.sorted.bam'.format(binDic['samtools'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # GenomeAnalysisTK:Lite-2.3.0-0
    cmds = ['java -jar {0} -T RealignerTargetCreator -R {1} -I {2}.mapped.sorted.bam -o {2}.intervals'.format(binDic['gatk'], binDic['ref'], args.sample),
            'java -jar {0} -T IndelRealigner -R {1} -I {2}.mapped.sorted.bam -targetIntervals {2}.intervals  -rf NotPrimaryAlignment  -o {2}.realigned.bam'.format(binDic['gatk'], binDic['ref'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # picard 1.133
    cmds = ['java -jar {0} MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT={1}.realigned.bam  OUTPUT={1}.mrd.bam  M={1}.metrics'.format(binDic['picard'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # samtools 1.2
    cmds = ['{0} index {1}.mrd.bam'.format(binDic['samtools'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # GenomeAnalysisTK:Lite-2.3.0-0
    cmds = ['java -jar {0} -T RealignerTargetCreator -R {1} -I {2}.mrd.bam -o {2}.2nd.intervals'.format(binDic['gatk'], binDic['ref'], args.sample),
            'java -jar {0} -T IndelRealigner -R {1} -I {2}.mrd.bam -targetIntervals {2}.2nd.intervals -rf NotPrimaryAlignment -o {2}.final.bam'.format(binDic['gatk'], binDic['ref'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # samtools 1.2
    cmds = ['{0} mpileup -d 10000 -L 10000 -B -t DP,DPR,DV,DP4,SP -f {1} {2}.final.bam -go {2}.raw.bcf'.format(binDic['samtools'], binDic['ref'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # bcftools 1.2
    cmds = ['{0} call -Avm -O z -f GQ -o {1}.var.raw.vcf.gz {1}.raw.bcf'.format(binDic['bcftools'], args.sample)]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    # picard 1.133
    cmds = ['java -jar {0} CollectMultipleMetrics VALIDATION_STRINGENCY=LENIENT INPUT={1}.final.bam OUTPUT={1}.cmm.metrics'.format(binDic['picard'], args.sample),
            'java -jar {0} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT INPUT={1}.final.bam OUTPUT={1}.cwm.metrics REFERENCE_SEQUENCE={2}'.format(binDic['picard'], args.sample, binDic['ref'])]

    for cmd in cmds:
        print cmd
        os.system(cmd)

    print 'DONE'

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--forward', help='forward reads')
    parser.add_argument('-r', '--reverse', help='reverse reads')
    parser.add_argument('-s', '--sample', help='sample name')
    args = parser.parse_args()
    main(args)
