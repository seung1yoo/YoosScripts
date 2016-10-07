
def main(args):
    print args
    lenDic = dict()
    for line in open(args.fia):
        items = line.strip().split('\t')
        lenDic.setdefault(items[0], int(items[1]))

    depthDic = dict()
    for line in open(args.depth):
        items = line.strip().split('\t')
        contig, position, depth = items
        depthDic.setdefault(contig, {}).setdefault('d', 0)
        depthDic[contig]['d'] += int(depth)
        depthDic.setdefault(contig, {}).setdefault('c', 0)
        depthDic[contig]['c'] += 1

    out = open('{0}.cov'.format(args.depth), 'w')
    for contig, length in lenDic.iteritems():
        if depthDic.has_key(contig):
            avg_depth = depthDic[contig]['d']/float(length)*100
            avg_cov = depthDic[contig]['c']/float(length)*100
            out.write('{0}\n'.format('\t'.join([str(x) for x in [contig, length, avg_depth, avg_cov]])))
        else:
            out.write('{0}\n'.format('\t'.join([str(x) for x in [contig, length, 0, 0]])))
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fia', default='Rockfish.polished_assembly.fasta.fai')
    parser.add_argument('-d', '--depth', default='Rockfish.polished_assembly.fasta_vs_MiSeqs.sort.bam.depth')
    args = parser.parse_args()
    main(args)
