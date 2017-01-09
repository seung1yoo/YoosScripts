

def main(args):
    print args

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fq1', '--fastq1', help='input forward fastq')
    parser.add_argument('-fq2', '--fastq2', help='input reverse fastq')
    parser.add_argument('-idx', '--bwaidx', help='bwa index name')
    parser.add_argument('-b', '--bwa', help='bwa full path')
    parser.add_argument('-s', '--samtools', help='samtools full path')
    parser.add_argument('-bf', '--bam2fastq', help='bam2fastq full path')
    parser.add_argument('-o', '--outdir', help='outdir full path')
    parser.add_argument('-s', '--sample', help='sample name')
    args = parser.parse_args()
    main(args)

