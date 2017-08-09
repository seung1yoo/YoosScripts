import os

class VCF_REDUCER():
    def __init__(self, invcf, outvcf, distance):
        inDic = self.inDicMaker(invcf)
        self.dicStatistics(inDic, 'origin')
        #
        outDic = self.outDicMaker(inDic, distance)
        self.dicStatistics(outDic, 'reduced')
        #
        self.vcf_writer(outDic, invcf, outvcf)

    def inDicMaker(self, invcf):
        inDic = dict()
        for line in open(invcf):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            pos = int(items[1])
            inDic.setdefault(chrom, {}).setdefault(pos, '')
        return inDic

    def outDicMaker(self, inDic, distance):
        outDic = dict()
        for chrom, posDic in inDic.iteritems():
            ## SELECT ONLY CHR
            if not chrom.startswith('chr'):
                continue
            ##
            pre_pos = sorted(posDic.keys())[0]
            for pos, temp in sorted(posDic.iteritems()):
                dist = pos - pre_pos
                if dist in [0]:
                    outDic.setdefault(chrom, {}).setdefault(pos, '')
                elif dist >= distance:
                    outDic.setdefault(chrom, {}).setdefault(pos, '')
                    pre_pos = pos
                else:
                    pass
            #
        return outDic

    def dicStatistics(self, aDic, tag):
        total_pos = 0
        for chrom, posDic in sorted(aDic.iteritems()):
            pos_n = len(posDic.keys())
            print '#{0} ==> {1} / {2}'.format(tag, chrom, pos_n)
            total_pos += pos_n
        print '#{0} ==> total / {1}'.format(tag, total_pos)

    def vcf_writer(self, outDic, invcf, outvcf):
        out = open(outvcf, 'w')
        for line in open(invcf):
            if line.startswith('#'):
                out.write(line)
                continue
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            pos = int(items[1])
            if chrom in outDic:
                if pos in outDic[chrom]:
                    out.write(line)
                else:
                    pass
            else:
                pass
        out.close()


def main(args):
    print args
    reducer = VCF_REDUCER(args.invcf, args.outvcf, args.distance)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-iv', '--invcf',
            default='../multisample.snpeff.dev.vcf')
    parser.add_argument('-ov', '--outvcf',
            default='multisample.snpeff.selected.dev.vcf')
    parser.add_argument('-d', '--distance', type=int,
            default=10000)
    args = parser.parse_args()
    main(args)

