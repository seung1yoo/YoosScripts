
def fileFinder(path, extension, firstTag, secondTag):
    cmd = "find {0} -iname '*.{1}'".format(path, extension)
    files = subprocess.check_output(cmd, shell=True)
    for afile in files.split('\n'):
        if not afile:
            continue
        sample_name = afile.split('/')[-1]
        print sample_name
        if firstTag in afile:
            print afile
        elif secondTag in afile:
            print afile
        else:
            print "check the firstTag, secondTag arguments"
            import sys
            sys.exit()

def main(args):
    print(args)
    fileDic = fileFinder(args.path, args.extension, args.firstTag, args.secondTag)

if __name__=='__main__':
    import glob
    import subprocess
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', help='The PATH for file search',
            default='/BiO/BioProjects/TBD160661/rawdata/1')
    parser.add_argument('-e', '--extension', help='File extension to search',
            default='fastq.gz')
    parser.add_argument('-s', '--script', help='FasterFastqStatistics',
            default='/home/siyoo/YoosScripts/tbi-service/SequenceInfo/FasterFastqStatistics')
    parser.add_argument('-o', '--outFile', help='Summary file',
            default='SequenceInfo.tbi.py.summary.xls')
    parser.add_argument('-1', '--firstTag', help='first read tag keyword',
            default='R1')
    parser.add_argument('-2', '--secondTag', help='second read tag keyword',
            default='R2')
    args = parser.parse_args()
    main(args)
