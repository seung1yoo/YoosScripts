def makeThreshold(sizes, cutoff):
    modi_sizes = []
    for size in sizes:
        if int(size) > cutoff:
            modi_sizes.append(cutoff)
        else:
            modi_sizes.append(size)
    return modi_sizes

def listmaker(file, col):
    units = [float(line.strip().split('\t')[col]) for line in open(file)]
    print 'average', sum(units)/len(units)
    sizes = makeThreshold(units, 4000) ### of histogram
    pylab.hist(sizes, bins=50, color=['crimson'])
    pylab.title('{0} sequences\nDepth {1} to {2}'.format(len(units), min(units), max(units)))
    pylab.xlabel('Depth')
    pylab.ylabel('Count')
    pylab.show()

def main(args):
    print args
    listmaker(args.file, col=2) # index.

if __name__=='__main__':
    import pylab ## black server OK, ngs2 server not OK
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help = 'specify coverage file as the output of genomeCovFromNonModel.py ',
            default='/tier3/Codes/niab/4.Centipede_2015/7.IncoMapping/1.ContamSeq_vs_Genome/ContamSeq_vs_Genome.sort.bam.depth.cov')
    args = parser.parse_args()
    main(args)
