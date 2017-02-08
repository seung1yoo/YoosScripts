def sizes_modifier(sizes, cutoff):
    modi_sizes = []
    for size in sizes:
        if int(size) > cutoff:
            modi_sizes.append(cutoff)
        else:
            modi_sizes.append(size)
    return modi_sizes

def lengthDistribution(fasta):
    sizes = [len(record) for record in SeqIO.parse(fasta, 'fasta')]
    print (sum(sizes)/len(sizes))
    max_size = sorted(sizes)[-1]
    sizes = sizes_modifier(sizes, 1000)
    pylab.hist(sizes, bins=50, color=['crimson'])
    pylab.title('{0} sequences\nLength {1} to {2}'.format(len(sizes), min(sizes), max_size))
    pylab.xlabel('Sequence length (bp)')
    pylab.ylabel('Count')
    pylab.show()

def gcContents(fasta):
    gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse(fasta, 'fasta'))
    #pylab.plot(gc_values)
    print (sum(gc_values)/len(gc_values))
    pylab.hist(gc_values, bins=60)
    pylab.title('{0} sequences\nGC% {1} to {2}'.format(len(gc_values), round(min(gc_values),2), round(max(gc_values),2)))
    pylab.xlabel('GC%')
    pylab.ylabel('Genes')
    pylab.show()

def main(args):
    print(args)
    lengthDistribution(args.fasta)
    gcContents(args.fasta)

if __name__=='__main__':
    import pylab ## black server OK, ngs2 server not OK
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help = 'specify fasta file',
            default='/tier1/Codes/nfrdi/2.HaliotisDiscus_2015/16.IncoGDB/Haliotis.IncoGDB.CDS.fna')
    args = parser.parse_args()
    main(args)
