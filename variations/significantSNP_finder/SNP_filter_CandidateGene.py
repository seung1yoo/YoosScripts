
def main(args):
    geneDic = dict()
    for gene in args.genes:
        geneDic.setdefault(gene, {})

    for line in open(args.invcf):
        if line.startswith('#CHROM'):
            title = line
            continue
        elif line.startswith('#'):
            continue
        #
        items = line.rstrip('\n').split('\t')
        #
        chrom = items[0]
        pos = items[1]
        ref_allele = items[3]
        alt_allele = items[4]
        primary_key = '{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele)
        #
        infos = items[7].split(';')
        for info in infos:
            key = info.split('=')[0]
            value = info.split('=')[-1]
            if key in ['ANN']:
                units = value.split('|')
                for unit in units:
                    if unit in geneDic:
                        geneDic[unit].setdefault(primary_key, line)
        #

    out = open('5.Candidate_genes.xls', 'w')
    out.write('Gene\t{0}'.format(title))
    for gene, keyDic in geneDic.iteritems():
        for key, line in keyDic.iteritems():
            out.write('{0}\t{1}'.format(gene, line))
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--invcf', default='2.brachygnathia_A-specific.vcf')
    parser.add_argument('-g', '--genes', nargs='+', default=[
        'ENSCAFG00000005970', 'ENSCAFG00000012374', 'ENSCAFG00000014993', 'ENSCAFG00000002506', 'ENSCAFG00000016731',
        'ENSCAFG00000002528', 'ENSCAFG00000011628', 'ENSCAFG00000008894', 'ENSCAFG00000019985', 'ENSCAFG00000002589',
        'ENSCAFG00000002997', 'ENSCAFG00000009059', 'ENSCAFG00000006755', 'ENSCAFG00000001872', 'ENSCAFG00000014548',
        'ENSCAFG00000018864', 'ENSCAFG00000018857', 'ENSCAFG00000009872', 'ENSCAFG00000002459', 'ENSCAFG00000003535',
        'ENSCAFG00000008889', 'ENSCAFG00000002191', 'ENSCAFG00000002187'])
    args = parser.parse_args()
    main(args)
