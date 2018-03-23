

def main(args):
    out = open(args.outGenesFn, 'w')
    for line in open(args.inGenesFn):
        if line.startswith('Order'):
            out.write(line)
            continue
        items = line.rstrip('\n').split('\t')
        if 'Y' in [items[args.idx]]:
            out.write(line)
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ig', '--inGenesFn', default='genes.kegg.13_1d_vs_20_1d.xls')
    parser.add_argument('-og', '--outGenesFn', default='genes.kegg.13_1d_vs_20_1d.Y.xls')
    parser.add_argument('-idx', '--idx', type=int, default=17)
    args = parser.parse_args()
    main(args)
