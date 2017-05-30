

class MUTANTCOMPARE():
    def __init__(self):
        pass

    def compare_variant(self, invcf, sample1, sample2):
        self.sample1 = sample1
        self.sample2 = sample2
        self.commonDic = dict()
        self.specificDic = dict()
        self.missingDic = dict()
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('#CHROM'):
                n_col = len(items)
                print len(items)
                print items
                samples = items[9:]
                samIdxDic = dict()
                for idx, sample in enumerate(samples):
                    if sample in [self.sample1, self.sample2]:
                        print idx, sample
                        samIdxDic.setdefault(sample, idx)
                continue
            if not len(items) in [n_col]:
                print '#ERROR'
                print len(items)
                print items
                import sys
                sys.exit()
            #print samIdxDic
            #
            formatTag = items[8]
            tagIdxDic = dict()
            for idx, atag in enumerate(formatTag.split(':')):
                tagIdxDic.setdefault(atag, idx)
            #print tagIdxDic
            #
            chrom = items[0]
            pos = items[1]
            ref_allele = items[3]
            alt_allele = items[4]
            samples = items[9:]
            muDic = dict()
            for sample, idx in samIdxDic.iteritems():
                sampleTags = samples[samIdxDic[sample]].split(':')
                #
                if sampleTags[tagIdxDic['GT']] in ['./.']:
                    muDic.setdefault('GT', []).append(sampleTags[tagIdxDic['GT']])
                    continue
                else:
                    muDic.setdefault('GT', []).append(sampleTags[tagIdxDic['GT']])
                ###
                genotypes = sampleTags[tagIdxDic['GT']].split('/')
                ref_gt = genotypes[0]
                alt_gt = genotypes[1]
                #
                depths = sampleTags[tagIdxDic['AD']].split(',')
                ref_depth = depths[int(ref_gt)]
                alt_depth = depths[int(alt_gt)]
                total_depth = sampleTags[tagIdxDic['DP']]
                ###
            #
            #print muDic
            if './.' in muDic['GT']:
                self.missingDic.setdefault('{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele), '')
            elif muDic['GT'][0] in [muDic['GT'][1]]:
                self.commonDic.setdefault('{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele), '')
            elif not muDic['GT'][0] in [muDic['GT'][1]]:
                self.specificDic.setdefault('{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele), '')
            else:
                print '#ERROR'
                print muDic
                import sys
                sys.exit()



    def vcf_maker(self, invcf, outprefix):
        print '{0} versus {1} : Common : {2}'.format(self.sample1, self.sample2, len(self.commonDic))
        print '{0} versus {1} : Specific : {2}'.format(self.sample1, self.sample2, len(self.specificDic))
        print '{0} versus {1} : Missing : {2}'.format(self.sample1, self.sample2, len(self.missingDic))
        #
        out_common = open('{0}_{1}-versus-{2}_common.vcf'.format(outprefix, self.sample1, self.sample2), 'w')
        out_specific = open('{0}_{1}-versus-{2}_specific.vcf'.format(outprefix, self.sample1, self.sample2), 'w')
        out_missing = open('{0}_{1}-versus-{2}_missing.vcf'.format(outprefix, self.sample1, self.sample2), 'w')
        #
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('#CHROM'):
                out_common.write(line)
                out_specific.write(line)
                out_missing.write(line)
                continue
            chrom = items[0]
            pos = items[1]
            ref_allele = items[3]
            alt_allele = items[4]
            primary_key = '{0}:{1}:{2}:{3}'.format(chrom, pos, ref_allele, alt_allele)
            if primary_key in self.commonDic:
                out_common.write(line)
            elif primary_key in self.specificDic:
                out_specific.write(line)
            elif primary_key in self.missingDic:
                out_missing.write(line)
            else:
                print '#ERROR'
                print line
                import sys
                sys.exit()
        out_common.close()
        out_specific.close()

def main(args):
    print args
    mc = MUTANTCOMPARE()
    mc.compare_variant(args.invcf, args.sample1, args.sample2)
    mc.vcf_maker(args.invcf, args.outprefix)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--invcf', default='SNP_filter_result.bi.vcf')
    parser.add_argument('-o', '--outprefix', default='SNP_filter_result')
    parser.add_argument('-s1', '--sample1', default='NT_1_SEQ')
    parser.add_argument('-s2', '--sample2', default='NT_3')
    args = parser.parse_args()
    main(args)
