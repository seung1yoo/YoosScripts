
import os
import sys

def fileChecker(ids, extDic, vcf):
    for id in ids:
        for cate, exts in extDic.iteritems():
            for ext in exts:
                afile = '{0}{1}'.format(id, ext)
                if not os.path.isfile(afile):
                    print id, cate, afile
                    return 1
                elif os.path.isfile(afile) and os.path.getsize(afile) < 10:
                    print id, cate, afile
                    return 1
    #
    if not os.path.isfile(vcf):
        print vcf
        return 1
    elif os.path.isfile(vcf) and os.path.getsize(vcf) < 10:
        print vcf
        return 1
    #
    return 0

def linker(ids, extDic, vcf):
    if not os.path.isdir('Upload'):
        os.makedirs('Upload')
    #
    for id in ids:
        for cate, exts in extDic.iteritems():
            targetDir = 'Upload/{0}/{1}'.format(id, cate)
            if not os.path.isdir(targetDir):
                os.makedirs(targetDir)
            for ext in exts:
                afile = os.path.abspath('{0}{1}'.format(id, ext))
                cmd = 'ln -s {0} {1}'.format(afile, targetDir)
                print cmd
                os.system(cmd)
    #
    cmd = 'ln -s {0} {1}'.format(os.path.abspath(vcf), 'Upload')
    print cmd
    os.system(cmd)


def main(args):
    extDic = {'01.CLEAN':['_1.clean.paired.fq.gz',
                          '_2.clean.paired.fq.gz'],
              '02.BAM':  ['.final.bam',
                          '.final.bai',
                          '.flagstat.txt'],
              '03.CMM':  ['.cmm.metrics.alignment_summary_metrics',
                          '.cmm.metrics.base_distribution_by_cycle.pdf',
                          '.cmm.metrics.base_distribution_by_cycle_metrics',
                          '.cmm.metrics.insert_size_histogram.pdf',
                          '.cmm.metrics.insert_size_metrics',
                          '.cmm.metrics.quality_by_cycle.pdf',
                          '.cmm.metrics.quality_by_cycle_metrics',
                          '.cmm.metrics.quality_distribution.pdf',
                          '.cmm.metrics.quality_distribution_metrics'],
              '04.CWM':  ['.cwm.metrics'],
              '05.UNMAP':['.unmapped_pair_R1.fq.gz',
                          '.unmapped_pair_R2.fq.gz',
                          '.unmapped_single_R1.fq.gz',
                          '.unmapped_single_R2.fq.gz']}

    vcf = 'multisamples.var.raw.vcf.gz'
    #
    fileCheck = fileChecker(args.ids, extDic, vcf)
    if fileCheck:
        print 'file checking'
        sys.exit()
    #
    linker(args.ids, extDic, vcf)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-is', '--ids', nargs='+',
        default=['Ig15255','Ig15256','Ig15257'])
    args = parser.parse_args()
    main(args)


