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

def filter_badvar(invcf, prefix):
    prefix = '{0}.good'.format(prefix)
    outvcf = '{0}.vcf'.format(prefix)
    min_DP = 10.0
    min_GQ = 60.0
    min_rate = 40.0
    max_rate = 60.0
    homo_min_rate = 90.0
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
                if int(total_depth) < min_DP:
                    new_samples.append('./.:.:.:.:.')
                    continue
                #
                genotype_quality = sampleTags[idxDic['GQ']]
                if int(genotype_quality) < min_GQ:
                    new_samples.append('./.:.:.:.:.')
                    continue
                #

                if ref_gt in [alt_gt]: ##This is homo
                    depths = sampleTags[idxDic['AD']].split(',')
                    homo_depth = depths[int(ref_gt)]
                    total_depth = sampleTags[idxDic['DP']]
                    homo_rate = int(homo_depth) / float(total_depth) * 100.0
                    if homo_rate >= homo_min_rate:
                        new_samples.append(sample)
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
                        new_samples.append(sample)
                    else:
                        new_samples.append('./.:.:.:.:.')
            #
            new_items = items[:9]
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
        print '{0} OK'.format(outxls)
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


def main(args):
    print args
    avcf = args.vcf
    prefix = args.outprefix
    #
    avcf, prefix = select_snp(avcf, args.vcftools, prefix)
    avcf, prefix = select_biallele(avcf, prefix)
    avcf, prefix = filter_badvar(avcf, prefix)
    avcf, prefix = filter_diffhomo(avcf, prefix, args.sample1, args.sample2)
    #
    outxls, prefix = vcf2xls(avcf, prefix, args.sample1, args.sample2)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='specify input vcf file name',
            default='multisample.1.recode.vcf')
    parser.add_argument('-vtp', '--vcftools', help='vcftools path',
            default="/BiO/BioTools/vcftools/vcftools_0.1.11/bin/vcftools")
    parser.add_argument('-pp', '--plink', help='plink v1.9 path',
            default='plink')
    parser.add_argument('-o', '--outprefix', help='output file name prefix',
            default='SNP_filter_result')
    parser.add_argument('-s1', '--sample1')
    parser.add_argument('-s2', '--sample2')
    args = parser.parse_args()
    main(args)
