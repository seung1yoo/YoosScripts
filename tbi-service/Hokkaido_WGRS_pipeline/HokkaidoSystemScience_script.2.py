
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
    ###################
    # samtools 1.2
    cmds = ['{0} mpileup -d 10000 -L 10000 -B -t DP,DPR,DV,DP4,SP -f {1} {2} -go multisamples.raw.bcf'.format(binDic['samtools'], binDic['ref'], ' '.join(args.finalBams))]
    for cmd in cmds:
        print cmd
        os.system(cmd)
        pass

    # bcftools 1.2
    cmds = ['{0} call -Avm -O z -f GQ -o multisamples.var.raw.vcf.gz multisamples.raw.bcf'.format(binDic['bcftools'])]
    for cmd in cmds:
        print cmd
        os.system(cmd)
        pass
    ###################

    print 'DONE'

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fbs', '--finalBams', nargs='+')
    args = parser.parse_args()
    main(args)


