
def main(args):
    print(args)
    out = open('{0}.wego'.format(args.xls), 'w')
    out.write('!WGOP_samp_inpu_1\n')
    for line in open(args.xls):
        items = line.rstrip('\n').split('\t')
        if line.startswith(args.title):
            go_idx = items.index(args.go_title)
            id_idx = items.index(args.id_title)
            print(go_idx)
            continue
        trId = items[id_idx]
        gos = items[go_idx].strip('"').split(',')
        if gos == ['-']:
            continue
        else:
            out.write('{0}\t{1}\n'.format(trId, '\t'.join(gos)))
    out.close()
    print("DONE")


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xls', help='xls file from RINE')
    parser.add_argument('-t', '--title', help='what is start chr of title line')
    parser.add_argument('-g', '--go_title', help='go col title')
    parser.add_argument('-i', '--id_title', help='id col title')
    args = parser.parse_args()
    main(args)
