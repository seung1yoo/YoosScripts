
def cutadapt_exe(cutLen, outputFastq, inputFastq):
    cmd = 'cutadapt --cut {0} -o {1} {2}'.format(cutLen, outputFastq, inputFastq)
    print cmd
    os.system(cmd)

def fileFinder(path, extension):
    cmd = 'find {0} -iname *.{1} > ./fileFinder.result'.format(path, extension)
    print cmd
    os.system(cmd)
    files = []
    for line in open('./fileFinder.result'):
        files.append(line.strip())
    return files

def main(args):
    files = fileFinder(args.path, args.extension)
    for file in files:
        inputFastq = file
        outputFastq = '{0}{1}.{2}'.format(file.split(args.extension)[0],
                str(int(args.oriLen)-int(args.cutLen)), args.extension)
        #print inputFastq, outputFastq
        cutadapt_exe(args.cutLen, outputFastq, inputFastq)


if __name__=='__main__':
    #from Bio import SeqIO
    import glob
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path',
            default='/BiO/BioProjects/TBD160661/rawdata/1')
    parser.add_argument('-e', '--extension',
            default='fastq.gz')
    parser.add_argument('-c', '--cutLen',
            default=50)
    parser.add_argument('-o', '--oriLen',
            default=151)
    args = parser.parse_args()
    main(args)
