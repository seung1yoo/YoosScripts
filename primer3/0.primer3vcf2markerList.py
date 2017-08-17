

def main(args):
    out = open(args.outFile, 'w')

    for line in open(args.inVcf):
        if line.startswith('#'):
            continue
        #
        items = line.rstrip('\n').split('\t')
        #
        chrom = items[0]
        pos = int(items[1])
        qual = items[5]
        ref_n = len(items[3])
        alt_n = len(items[4])
        ref = items[3]
        alt = items[4]
        start = pos
        end = pos + ref_n - 1
        #
        depth = '0'
        for info in items[7].split(';'): # INFO
            if info.startswith('DP'):
                units = info.split('=')
                key = units[0]
                value = units[1]
                if key in ['DP']:
                    depth = value
                else:
                    pass
            else:
                pass
        #
        new_items = [chrom, str(start), str(end), qual, depth, ref, alt]
        out.write('{0}\n'.format('\t'.join(new_items)))
    out.close()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-iv', '--inVcf')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()
    main(args)
