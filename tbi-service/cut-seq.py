

def main(args):
    print args

if __name__=='__main__':
    #from Bio import SeqIO
    #import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputFastq',
            default='/home/siyoo/BioProject/TBD160661/TN1609R0051/TN1609R0051--ATCACGAT-1_S24_L004_R2_001.fastq.gz')
    parser.add_argument('-l', '--length',
            default=101)
    args = parser.parse_args()
    main(args)

