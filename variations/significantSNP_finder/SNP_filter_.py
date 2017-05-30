import os

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

def select_depth(invcf, depth, vcftools, prefix):
    prefix = '{0}.DP{1}'.format(prefix, depth)
    outvcf = '{0}.recode.vcf'.format(prefix)
    if os.path.isfile(outvcf):
        print '{0} OK'.format(outvcf)
    else:
        cmd = '{0} --vcf {1} --depth --out {2}.depth'.format(vcftools, invcf, prefix.split('.DP')[0])
        print cmd
        os.system(cmd)
        #
        cmd = '{0} --vcf {1} --min-meanDP {2} --recode --recode-INFO-all --out {3}'.format(vcftools, invcf, depth, prefix)
        print cmd
        os.system(cmd)

    return outvcf, prefix

def make_plinkInput(invcf, vcftools, prefix):
    plink_fileset = '{0}.p'.format(prefix)
    #
    #invcf4plink = '{0}.vcf'.format(plink_fileset)
    #out = open(invcf4plink, 'w')
    #chrMapDic = dict()
    #chrMapId = 0
    #for line in open(invcf):
    #    if line.startswith('#'):
    #        out.write(line)
    #        continue
    #    items = line.rstrip('\n').split('\t')
    #    chrom = items[0]
    #    if chrom in chrMapDic:
    #        items[0] = chrMapDic[chrom]
    #        out.write('{0}\n'.format('\t'.join(items)))
    #    else:
    #        chrMapId += 1
    #        chrMapDic.setdefault(chrom, str(chrMapId))
    #        items[0] = chrMapDic[chrom]
    #        out.write('{0}\n'.format('\t'.join(items)))
    #out.close()
    #
    #chrMap4plink = '{0}.chrMap'.format(plink_fileset)
    #out = open(chrMap4plink, 'w')
    #for chrom, chrMapId in sorted(chrMapDic.iteritems()):
    #    out.write('{0}\n'.format('\t'.join([chrom, chrMapId])))
    #out.close()
    #
    pedFile = '{0}.ped'.format(plink_fileset)
    mapFile = '{0}.map'.format(plink_fileset)
    #
    if os.path.isfile(pedFile) and os.path.isfile(mapFile):
        print '{0} OK'.format(pedFile)
        print '{0} OK'.format(mapFile)
    else:
        cmd = '{0} --vcf {1} --plink --out {2}'.format(vcftools, invcf, plink_fileset)
        print cmd
        os.system(cmd)
    #
    mapTempFile = '{0}.temp'.format(mapFile)
    os.system('mv {0} {1}'.format(mapFile, mapTempFile))
    out = open(mapFile, 'w')
    for line in open(mapTempFile):
        items = line.rstrip('\n').split('\t')
        snp_id, position = items[1].split(':')
        items[0] = snp_id
        items[1].replace(':', '_')
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()

    return pedFile, mapFile, plink_fileset

def make_bed(plink_fileset, plink):
    binary_fileset = '{0}b'.format(plink_fileset)
    bedFile = '{0}.bed'.format(binary_fileset)
    bimFile = '{0}.bim'.format(binary_fileset)
    famFile = '{0}.fam'.format(binary_fileset)
    #
    if os.path.isfile(bedFile) and os.path.isfile(bimFile) and os.path.isfile(famFile):
        print '{0} OK'.format(bedFile)
        print '{0} OK'.format(bimFile)
        print '{0} OK'.format(famFile)
    else:
        cmd = '{0} --file {1} --make-bed --allow-extra-chr --out {2}'.format(plink, plink_fileset, binary_fileset)
        print cmd
        os.system(cmd)

    return bedFile, bimFile, famFile, binary_fileset

def filter_geno(binary_fileset, plink):
    geno_fileset = '{0}.geno'.format(binary_fileset)
    bedFile = '{0}.bed'.format(geno_fileset)
    bimFile = '{0}.bim'.format(geno_fileset)
    famFile = '{0}.fam'.format(geno_fileset)
    #
    if os.path.isfile(bedFile) and os.path.isfile(bimFile) and os.path.isfile(famFile):
        print '{0} OK'.format(bedFile)
        print '{0} OK'.format(bimFile)
        print '{0} OK'.format(famFile)
    else:
        cmd = '{0} --bfile {1} --geno 0.2 --make-bed --allow-extra-chr --out {2}'.format(plink, binary_fileset, geno_fileset)
        print cmd
        os.system(cmd)

    return bedFile, bimFile, famFile, geno_fileset

def filter_maf(geno_fileset, plink):
    maf_fileset = '{0}.maf'.format(geno_fileset)
    bedFile = '{0}.bed'.format(maf_fileset)
    bimFile = '{0}.bim'.format(maf_fileset)
    famFile = '{0}.fam'.format(maf_fileset)
    #
    if os.path.isfile(bedFile) and os.path.isfile(bimFile) and os.path.isfile(famFile):
        print '{0} OK'.format(bedFile)
        print '{0} OK'.format(bimFile)
        print '{0} OK'.format(famFile)
    else:
        cmd = '{0} --bfile {1} --maf 0.05 --make-bed --allow-extra-chr --out {2}'.format(plink, geno_fileset, maf_fileset)
        print cmd
        os.system(cmd)

    return bedFile, bimFile, famFile, maf_fileset

def filter_hwe(maf_fileset, plink):
    hwe_fileset = '{0}.hwe'.format(maf_fileset)
    bedFile = '{0}.bed'.format(hwe_fileset)
    bimFile = '{0}.bim'.format(hwe_fileset)
    famFile = '{0}.fam'.format(hwe_fileset)
    #
    if os.path.isfile(bedFile) and os.path.isfile(bimFile) and os.path.isfile(famFile):
        print '{0} OK'.format(bedFile)
        print '{0} OK'.format(bimFile)
        print '{0} OK'.format(famFile)
    else:
        cmd = '{0} --bfile {1} --hwe 0.05 --make-bed --allow-extra-chr --out {2}'.format(plink, maf_fileset, hwe_fileset)
        print cmd
        os.system(cmd)

    return bedFile, bimFile, famFile, hwe_fileset

def test_association(hwe_fileset, plink):
    out_fileset = '{0}.assoc'.format(hwe_fileset)
    #
    cmd = '{0} --bfile {1} --assoc --allow-no-sex --allow-extra-chr --out {2}'.format(plink, hwe_fileset, out_fileset)
    print cmd
    os.system(cmd)
    #
    cmd = '{0} --bfile {1} --assoc fisher --allow-no-sex --allow-extra-chr --out {2}'.format(plink, hwe_fileset, out_fileset)
    #cmd = '{0} --bfile {1} --assoc fisher-midp --allow-no-sex --out {2}'.format(plink, hwe_fileset, out_fileset)
    print cmd
    os.system(cmd)

def select_hetero(invcf, prefix):
    prefix = '{0}.hete'.format(prefix)
    outvcf = '{0}.vcf'.format(prefix)
    min_rate = 30.0
    max_rate = 50.0
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
            out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
        out.close()

    return outvcf, prefix

def main(args):
    print args
    avcf = args.vcf
    prefix = args.outprefix
    #
    avcf, prefix = select_snp(avcf, args.vcftools, prefix)
    avcf, prefix = select_depth(avcf, '10', args.vcftools, prefix)
    avcf, prefix = select_biallele(avcf, prefix)
    avcf, prefix = select_hetero(avcf, prefix)
    #
    pedFile, mapFile, plink_fileset = make_plinkInput(avcf, args.vcftools, prefix)
    bedFile, bimFile, famFile, binary_fileset = make_bed(plink_fileset, args.plink)
    # manual step and re-run
    famOriFile = '{0}.Ori.fam'.format(binary_fileset)
    if not os.path.isfile(famOriFile):
        print 'plz modify the fam file first.'
        print '    Step 1) cp SNP_filter_result.SNP.pb.fam SNP_filter_result.SNP.pb.Ori.fam'
        print '    Step 2) vim SNP_filter_result.SNP.pb.fam'
        import sys
        sys.exit()
    #
    bedFile, bimFile, famFile, geno_fileset = filter_geno(binary_fileset, args.plink)
    bedFile, bimFile, famFile, maf_fileset = filter_maf(geno_fileset, args.plink)
    bedFile, bimFile, famFile, hwe_fileset = filter_hwe(maf_fileset, args.plink)
    test_association(hwe_fileset, args.plink)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='specify input vcf file name',
            default='../multisample.snpeff.vcf')
    parser.add_argument('-vtp', '--vcftools', help='vcftools path',
            default='/BiO/BioTools/vcftools/current/bin/vcftools')
    parser.add_argument('-pp', '--plink', help='plink v1.9 path',
            #default='/BiO/BioTools/plink/v1.09/plink')
            default='/BiO/BioTools/plink/1.07/plink')
    parser.add_argument('-o', '--outprefix', help='output file name prefix',
            default='SNP_filter_result')
    args = parser.parse_args()
    main(args)
