


def main(args):
    print(args)

if __name__=='__main__':
    import glob
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', help='The PATH for file search',
            default='/BiO/BioProjects/TBD160661/rawdata/1')
    parser.add_argument('-e', '--extension', help='File extension to search',
            default='fastq.gz')
    parser.add_argument('-s', '--script', help='FasterFastqStatistics',
            default='/home/shsong/work/Pipeline/dnaseq/util/FasterFastqStatistics')
    parser.add_argument('-o', '--outFile', help='Summary file',
            default='SequenceInfo.tbi.py.summary.xls')
    args = parser.parse_args()
    main(args)

