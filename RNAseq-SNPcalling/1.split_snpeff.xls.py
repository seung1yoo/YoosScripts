
def make_outfile_handle_dic(prefix, chrs, depth):
    ofh_dic = dict()
    for achr in chrs:
        outfn = '{0}.DP{1}.{2}.xls'.format(prefix, depth, achr)
        outfh = open(outfn, 'w')
        ofh_dic.setdefault(achr, outfh)
    return ofh_dic

def close_outfile_handle(ofh_dic):
    for achr, outfh in ofh_dic.iteritems():
        outfh.close()

def main(args):
    ofh_dic = make_outfile_handle_dic(args.prefix, args.chrs, args.depth)
    #
    for line in open(args.xls):
        if line.startswith('#CHROM'):
            for achr, outfh in ofh_dic.iteritems():
                outfh.write(line)
            continue
        items = line.rstrip('\n').split('\t')
        achr = items[0]
        dp = items[5]
        if int(dp) < int(args.depth):
            continue
        if achr in args.chrs:
            ofh_dic[achr].write(line)
        elif not achr in args.chrs:
            continue
    #
    close_outfile_handle(ofh_dic)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrs', nargs='+',
        default=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','X','Y','MT'])
    parser.add_argument('--xls',
        default='multisample.snpeff.xls')
    parser.add_argument('--prefix',
        default='multisample.snpeff')
    parser.add_argument('--depth',
        default='10')
    args = parser.parse_args()
    main(args)
