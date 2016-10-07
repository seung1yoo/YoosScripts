def makeHisto(file, sample):
    covDic = dict()
    for line in open(file):
        items = line.strip().split('\t')
        if not len(items) in [16]:
            print items
            sys.exit()
        ref = items[0]
        depth, cov_base, len_gene, cov_ratio = items[-4:]

        covDic.setdefault(ref, {}).setdefault('num_gene', 0)
        covDic[ref]['num_gene'] += 1
        covDic.setdefault(ref, {}).setdefault('len_gene', 0)
        covDic[ref]['len_gene'] += int(len_gene)
        covDic.setdefault(ref, {}).setdefault('cov_base', 0)
        covDic[ref]['cov_base'] += int(cov_base)

    for ref, infoDic in sorted(covDic.iteritems()):
        print sample, ref, infoDic['num_gene'], infoDic['len_gene'], infoDic['cov_base'], infoDic['cov_base']/float(infoDic['len_gene'])*100

def main():
    files = glob.glob('*/*.TAIR10.bed.coverage')
    for n, file in enumerate(files):
        sample = file.split('/')[0]
        makeHisto(file, sample)


if __name__=='__main__':
    import glob
    import sys
    main()

