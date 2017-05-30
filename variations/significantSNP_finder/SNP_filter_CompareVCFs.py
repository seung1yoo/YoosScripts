
class MUTANTCOMPARE():
    def __init__(self):
        self.varDic = dict()
        self.vcfDic = dict()
        pass

    def variation_update(self, invcf, name):
        self.vcfDic.setdefault(name, invcf)
        #
        for n, line in enumerate(open(invcf)):
            if line.startswith('#CHROM'):
                self.vcftitle = line
                continue
            elif line.startswith('#'):
                continue
            #
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            pos = items[1]
            ref_allele = items[3]
            alt_allele = items[4]
            primary_key = '{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele)
            self.varDic.setdefault(primary_key, {}).setdefault(name, line)
            #
            if n in [100]:
                #break
                pass

    def variation_stats(self, names):
        self.countDic = dict()
        #
        for primary_key, nameDic in self.varDic.iteritems():
            #
            name_key = '-vs-'.join(sorted(nameDic.keys()))
            self.countDic.setdefault(name_key, {}).setdefault(primary_key, '')
        #
        total = 0
        for name_key, primary_keyDic in self.countDic.iteritems():
            total += len(primary_keyDic)
            print '#STATS : {0} vcf files has {1} variations.'.format(name_key, len(primary_keyDic))
        print '#STATS : Total {0} variations.'.format(total)

    def variation_outmake(self, outprefix):
        outDic = dict()
        for name_key, aDic in self.countDic.iteritems():
            outDic.setdefault(name_key, open('{0}_{1}.vcf'.format(outprefix, name_key), 'w'))
            outDic[name_key].write(self.vcftitle)
        #
        for primary_key, nameDic in self.varDic.iteritems():
            line = nameDic[nameDic.keys()[0]]
            name_key = '-vs-'.join(sorted(nameDic.keys()))
            outDic[name_key].write(line)
        #
        for name_key, out in outDic.iteritems():
            out.close()





def main(args):
    print args
    mc = MUTANTCOMPARE()
    mc.variation_update(args.vcf1, args.name1)
    mc.variation_update(args.vcf2, args.name2)
    mc.variation_stats([args.name1, args.name2])
    mc.variation_outmake(args.outprefix)
    #print mc.varDic


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v1', '--vcf1', default='2.brachygnathia_A-specific.vcf')
    parser.add_argument('-v2', '--vcf2', default='3.donor-cell_D.vcf')
    parser.add_argument('-n1', '--name1', default='brachygnathia_A-specific')
    parser.add_argument('-n2', '--name2', default='donor-cell_D')
    parser.add_argument('-o', '--outprefix', default='SNP_filter_result_5')
    args = parser.parse_args()
    main(args)
