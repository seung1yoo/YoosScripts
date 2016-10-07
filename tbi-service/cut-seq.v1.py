

def main(args):
    cmd = 'cutadapt --cut {0} -o {1} {2}'.format(args.cutLen, args.outputFastq, args.inputFastq)
    print cmd
    os.system(cmd)


if __name__=='__main__':
    #from Bio import SeqIO
    #import glob
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputFastq',
            default='/home/siyoo/BioProject/TBD160661/TN1609R0051/TN1609R0051--ATCACGAT-1_S24_L004_R2_001.fastq.gz')
    parser.add_argument('-o', '--outputFastq',
            default='/home/siyoo/BioProject/TBD160661/TN1609R0051/TN1609R0051--ATCACGAT-1_S24_L004_R2_001.101.fastq.gz')
    parser.add_argument('-c', '--cutLen',
            default=50)
    args = parser.parse_args()
    main(args)

