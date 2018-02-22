
def biotype_parser(gtf):
    bt2gDic = dict()
    g2btDic = dict()
    for line in open(gtf):
        if line.startswith('#'):
            continue
        items = line.rstrip('\n').split('\t')
        biotype = items[1] ### Checking point for biotype. in some case, biotype will be in 9th column.
        atts = items[8].rstrip(';').split(';')
        attDic = dict()
        for att in atts:
            att = att.strip()
            key, value = att.split(' ')
            value = value.strip('"')
            attDic.setdefault(key, value)
        bt2gDic.setdefault(biotype, {}).setdefault(attDic['gene_id'], attDic['gene_name'])
        g2btDic.setdefault(attDic['gene_id'], biotype)
    return bt2gDic, g2btDic

def outHeaderLineFinder(xls):
    idxDic = dict()
    for line in open(xls):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['Order', 'GeneId']:
            outHeaderLine = line
            for idx, item in enumerate(items):
                idxDic.setdefault(item, idx)
            break
        #
    return outHeaderLine, idxDic

def gene_classificator(xls, outHeaderLine, idxDic, bt2gDic, g2btDic, outprefix, extend):
    outDic = dict()
    for biotype, geneDic in bt2gDic.iteritems():
        out = open('{0}.{1}.{2}'.format(args.outprefix, biotype, extend), 'w')
        out.write(outHeaderLine)
        outDic.setdefault(biotype, out)
    #
    for line in open(xls):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['Order', 'GeneId']:
            continue
        outDic[g2btDic[items[idxDic['GeneAcc']]]].write(line)
    #
    for biotype, out in outDic.iteritems():
        out.close()

def main(args):
    print args
    bt2gDic, g2btDic = biotype_parser(args.gtf)
    #
    outHeaderLine, idxDic = outHeaderLineFinder(args.gene_xls)
    gene_classificator(args.gene_xls, outHeaderLine, idxDic, bt2gDic, g2btDic, args.outprefix, 'genes.xls')
    #
    outHeaderLine, idxDic = outHeaderLineFinder(args.de_xls)
    gene_classificator(args.de_xls, outHeaderLine, idxDic, bt2gDic, g2btDic, args.outprefix, 'genes.de.xls')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', default='/BiO/BioResources/References/Drosophila_melanogaster/D_melanogaster_ENS_23/D_melanogaster_ENS_23.chr.gtf')
    parser.add_argument('--gene-xls', default='/BiO/BioProjects/TBD170539-HUMC-Drosophila-RNAref-20180104/Rine_Quant/report/Assembly/genes.xls')
    parser.add_argument('--de-xls', default='/BiO/BioProjects/TBD170539-HUMC-Drosophila-RNAref-20180104/Rine_Quant/report/DEG/gene.de.xls')
    parser.add_argument('-op', '--outprefix', default='TBD170539-HUMC-Drosophila-RNAref-BioType')
    args = parser.parse_args()
    main(args)

