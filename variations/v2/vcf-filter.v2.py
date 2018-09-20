#!/usr/bin/python

import os



def select_snp(invcf, vcftools, prefix):
    prefix = '{0}.SNP'.format(prefix)
    outvcf = '{0}.recode.vcf'.format(prefix)
    if os.path.isfile(outvcf):
        print '{0} OK'.format(outvcf)
    else:
        cmd = '{0} --vcf {1} --remove-indels --recode --recode-INFO-all --out {2}'.format(vcftools, invcf, prefix)
        print cmd
        os.system(cmd)

    return outvcf, prefix

def select_biallele(invcf, prefix):
    prefix = '{0}.bi'.format(prefix)
    outvcf = '{0}.vcf'.format(prefix)
    if os.path.isfile(outvcf):
        print '{0} OK'.format(outvcf)
    else:
        out = open(outvcf, 'w')
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('#'):
                out.write(line)
                continue
            ref = items[3]
            alt = items[4]
            if ',' in alt:
                continue
            #if not len(alt) in [1] or not len(ref) in [1]:
            #    continue
            out.write(line)
        out.close()
    return outvcf, prefix

def select_goodvar(invcf, prefix):
    prefix = '{0}.good'.format(prefix)
    outvcf = '{0}.vcf'.format(prefix)
    min_DP = 10.0 # default 5.0
    min_GQ = 60.0 # default 60.0
    min_rate = 20.0
    max_rate = 80.0
    homo_min_rate = 0.0 # default 10.0
    if os.path.isfile(outvcf):
        print '{0} OK'.format(outvcf)
    else:
        out = open(outvcf, 'w')
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('#'):
                out.write(line)
                continue
            #
            formatTag = items[8]
            idxDic = dict()
            for idx, atag in enumerate(formatTag.split(':')):
                idxDic.setdefault(atag, idx)
            #
            samples = items[9:]
            new_samples = []
            for sample in samples:
                sampleTags = sample.split(':')
                #
                genotypes = sampleTags[idxDic['GT']].split('/')
                ref_gt = genotypes[0]
                alt_gt = genotypes[1]
                if ref_gt in ['.'] and alt_gt in ['.']:
                    new_samples.append('./.:.:.:.:.')
                    continue
                #
                total_depth = sampleTags[idxDic['DP']]
                if int(total_depth) <= min_DP:
                    new_samples.append('./.:.:.:.:.')
                    continue
                #
                genotype_quality = sampleTags[idxDic['GQ']]
                if int(genotype_quality) <= min_GQ:
                    new_samples.append('./.:.:.:.:.')
                    continue
                #

                if ref_gt in [alt_gt]: ##This is homo
                    depths = sampleTags[idxDic['AD']].split(',')
                    homo_depth = depths[int(ref_gt)]
                    total_depth = sampleTags[idxDic['DP']]
                    homo_rate = int(homo_depth) / float(total_depth) * 100.0
                    if homo_rate >= homo_min_rate:
                        #new_samples.append(sample)
                        #GT:AD:DP:GQ:PL
                        new_samples.append('{0}'.format(':'.join([sampleTags[idxDic['GT']],
                                                                  sampleTags[idxDic['AD']],
                                                                  sampleTags[idxDic['DP']],
                                                                  sampleTags[idxDic['GQ']],
                                                                  sampleTags[idxDic['PL']]])))
                    else:
                        new_samples.append('./.:.:.:.:.')
                else: ##This is hetoro
                    depths = sampleTags[idxDic['AD']].split(',')
                    ref_depth = depths[int(ref_gt)]
                    alt_depth = depths[int(alt_gt)]
                    #depth_rate = int(ref_depth) / float(alt_depth) * 100.0
                    depth_rate = int(alt_depth) / float(int(ref_depth)+int(alt_depth)) * 100.0
                    #
                    if min_rate <= depth_rate <= max_rate:
                        #new_samples.append(sample)
                        new_samples.append('{0}'.format(':'.join([sampleTags[idxDic['GT']],
                                                                  sampleTags[idxDic['AD']],
                                                                  sampleTags[idxDic['DP']],
                                                                  sampleTags[idxDic['GQ']],
                                                                  sampleTags[idxDic['PL']]])))
                    else:
                        new_samples.append('./.:.:.:.:.')
            #
            new_items = items[:8]
            new_items.append('GT:AD:DP:GQ:PL')
            if not len(samples) in [len(new_samples)]:
                print 'ERROR : check the samples and new_samples'
                import sys
                sys.exit()
            new_items.extend(new_samples)
            #
            if not './.:.:.:.:.' in new_items:
                out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
            else:
                continue
        out.close()

    return outvcf, prefix

def filter_diffhomo(invcf, prefix, sample_1, sample_2):
    prefix = '{0}.diffhomo.{1}_vs_{2}'.format(prefix, sample_1, sample_2)
    outvcf = '{0}.vcf'.format(prefix)
    if os.path.isfile(outvcf):
        print '{0} OK'.format(outvcf)
    else:
        out = open(outvcf, 'w')
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('#CHROM'):
                out.write(line)
                s_idxDic = dict()
                for idx, s_name in enumerate(items[9:]):
                    s_idx = idx
                    s_idxDic.setdefault(s_name, s_idx)
                continue
            elif line.startswith('#'):
                out.write(line)
                continue
            #
            formatTag = items[8]
            idxDic = dict()
            for idx, atag in enumerate(formatTag.split(':')):
                idxDic.setdefault(atag, idx)
            #
            samples = items[9:]
            gt_1 = samples[s_idxDic[sample_1]].split(':')[idxDic['GT']]
            gt_1_units = list(set(gt_1.split('/')))
            gt_2 = samples[s_idxDic[sample_2]].split(':')[idxDic['GT']]
            gt_2_units = list(set(gt_2.split('/')))
            if not len(gt_1_units) in [1] or not len(gt_2_units) in [1]:
                continue # only homo allele in both samples
            if gt_1 in [gt_2]:
                continue # only diff-homo alleles between samples
            #
            out.write('{0}\n'.format('\t'.join(items)))
        out.close()

    return outvcf, prefix

def vcf2xls(invcf, prefix, sample_1, sample_2):
    prefix = '{0}'.format(prefix)
    outxls = '{0}.xls'.format(prefix)
    if os.path.isfile(outxls):
        print '{0} OK'.format(outvcf)
    else:
        out = open(outxls, 'w')
        for line in open(invcf):
            items = line.strip('\n').split('\t')
            if line.startswith('##'):
                continue
            new_items = list()
            if line.startswith('#CHROM'):
                new_items.append(items[0]) #chrom
                new_items.append(items[1]) #pos
                new_items.append(items[3]) #ref
                new_items.append(items[4]) #alt
                new_items.extend(['SN={0}'.format(x) for x in items[9:]]) #samples for GT
                new_items.extend(['SN={0}'.format(x) for x in items[9:]]) #samples for AD
                new_items.extend(['SN={0}'.format(x) for x in items[9:]]) #samples for Origin
                out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
                #
                s_idxDic = dict()
                for idx, s_name in enumerate(items[9:]):
                    s_idx = idx
                    s_idxDic.setdefault(s_name, s_idx)
                #
                continue
            #
            formatTag = items[8]
            idxDic = dict()
            for idx, atag in enumerate(formatTag.split(':')):
                idxDic.setdefault(atag, idx)
            #
            samples = items[9:]
            new_samples = []
            for sample in samples:
                sampleTags = sample.split(':')
                new_samples.append('GT={0}'.format(sampleTags[idxDic['GT']]))
            #
            for sample in samples:
                sampleTags = sample.split(':')
                new_samples.append('AD={0}'.format(sampleTags[idxDic['AD']]))
            #
            gts = []
            for sample in samples:
                sampleTags = sample.split(':')
                gts.append(sampleTags[idxDic['GT']])
            for gt in gts:
                if gt in [samples[s_idxDic[sample_1]].split(':')[idxDic['GT']]]:
                    new_samples.append('ORIGIN={0}'.format(sample_1))
                elif gt in [samples[s_idxDic[sample_2]].split(':')[idxDic['GT']]]:
                    new_samples.append('ORIGIN={0}'.format(sample_2))
                else:
                    new_samples.append('ORIGIN=OTHER')
            #
            new_items.append(items[0]) #chrom
            new_items.append(items[1]) #pos
            new_items.append(items[3]) #ref
            new_items.append(items[4]) #alt
            new_items.extend(new_samples) #samples for GT and AD
            out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
        out.close()

    return outxls, prefix

class VarComparing:
    def __init__(self, vcf_fn):
        self.fixed_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        self.vcf_fn = vcf_fn
        self.idx_sample_dic = self.extract_sample_idx()
        self.idx_format_col = self.extract_format_idx()

    def extract_format_idx(self):
        aIdx = ''
        for line in open(self.vcf_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#CHROM']:
                for idx, item in enumerate(items):
                    if item in ['FORMAT']:
                        aIdx = idx
                break
        return aIdx

    def extract_sample_idx(self):
        aDic = dict()
        for line in open(self.vcf_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#CHROM']:
                for idx, item in enumerate(items):
                    if not item in self.fixed_cols:
                        aDic.setdefault(item, idx)
                break
        return aDic

    def format_indexing(self, items):
        format_col = items[self.idx_format_col]
        #
        aDic = dict()
        for idx, unit in enumerate(format_col.split(':')):
            aDic.setdefault(unit, idx)
        return aDic

    def call_gt(self, items, samples, idx_format_dic):
        gts = list()
        for sample in samples:
            info = items[self.idx_sample_dic[sample]]
            infos = info.split(':')
            gt = infos[idx_format_dic['GT']]
            gts.append(gt)
        return gts

    def is_hetero(self, gts):
        for gt in gts:
            ref, alt = gt.split('/')
            #
            if ref in ['.']:
                return 0
            if alt in ['.']:
                return 0
            #
            if not ref in [alt]:
                return 1
        return 0

    def is_homo(self, gts):
        for gt in gts:
            ref, alt = gt.split('/')
            #
            if ref in ['.']:
                return 0
            if alt in ['.']:
                return 0
            #
            if ref in [alt]:
                return 1
        return 0

    def have_same_gt(self, gts):
        if len(list(set(gts))) in [1]:
            return 1
        return 0

    def is_diff_a_b(self, a_gts, b_gts):
        if not a_gts[0] in [b_gts[0]]:
            return 1
        return 0

    def run(self, a_homo, b_homo, hete):
        a_homos = a_homo.split(',')
        b_homos = b_homo.split(',')
        hetes = hete.split(',')
        print a_homos, b_homos, hetes
        #
        out = open('{0}.VarComparing'.format(self.vcf_fn),'w')
        for line in open(self.vcf_fn):
            if line.startswith('#'):
                out.write(line)
                continue
            items = line.rstrip('\n').split('\t')
            idx_format_dic = self.format_indexing(items)
            #
            a_homo_gts = self.call_gt(items, a_homos, idx_format_dic)
            b_homo_gts = self.call_gt(items, b_homos, idx_format_dic)
            hete_gts = self.call_gt(items, hetes, idx_format_dic)
            #
            if not self.have_same_gt(a_homo_gts):
                continue
            if not self.have_same_gt(b_homo_gts):
                continue
            if not self.is_hetero(hete_gts):
                continue
            if not self.is_homo(a_homo_gts):
                continue
            if not self.is_homo(b_homo_gts):
                continue
            if not self.is_diff_a_b(a_homo_gts, b_homo_gts):
                continue
            out.write(line)
        out.close()


def main(args):
    print args
    avcf = args.vcf
    prefix = args.outprefix
    #
    avcf, prefix = select_snp(avcf, args.vcftools, prefix)
    avcf, prefix = select_biallele(avcf, prefix)
    avcf, prefix = select_goodvar(avcf, prefix)
    #
    ## below steps were made for TBD171076 project
    #avcf, prefix = filter_diffhomo(avcf, prefix, "1-HWANGGEUM", "2-DAEPUNG")
    #outxls, prefix = vcf2xls(avcf, prefix, "1-HWANGGEUM", "2-DAEPUNG")
    ## below steps were made for L_edodes WGRS project
    vc = VarComparing(avcf)
    print vc.idx_sample_dic
    print vc.idx_format_col
    vc.run('47','46,48,49','39')



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='specify input vcf file name',
            default='multisample.snpeff.vcf')
    parser.add_argument('-vtp', '--vcftools', help='vcftools path',
            default="/BiO/BioTools/vcftools/vcftools_0.1.13/bin/vcftools")
    parser.add_argument('-pp', '--plink', help='plink v1.9 path',
            default='plink')
    parser.add_argument('-o', '--outprefix', help='output file name prefix',
            default='L_edodes_newgeneset')
    args = parser.parse_args()
    main(args)
