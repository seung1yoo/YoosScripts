
def subprocess_open(command):
    popen = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (stdoutdata, stderrdata) = popen.communicate()
    return stdoutdata, stderrdata

def main(args):
    print args
    cufflinks = '/BiO/BioTools/cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks'
    gff = args.gtf
    genome = args.genomeFasta
    for line in open(args.configuration_file):
        items = line.strip().split(',')
        n, sample, bam_file = items
        cmd = '{0} -o {1} -p 10 -g {2} --library-type fr-unstranded --multi-read-correct --frag-bias-correct {3} {4} > {1}.cufflinks.log'.format(cufflinks, sample, gff, genome, bam_file)
        print cmd
        os.system(cmd)

    cmd = 'find ./ -iname transcripts.gtf > gtf_list.txt'
    print cmd
    os.system(cmd)

    cuffmerge = '/BiO/BioTools/cufflinks/cufflinks-2.2.1.Linux_x86_64/cuffmerge' 
    cmd = '{0} -p 10 -s {1} -g {2} gtf_list.txt'.format(cuffmerge, genome, gff)
    print cmd
    os.system(cmd)

    merged_gtf = './merged_asm/merged.gtf'
    cuffquant = '/BiO/BioTools/cufflinks/cufflinks-2.2.1.Linux_x86_64/cuffquant' 
    for line in open(args.configuration_file):
        items = line.strip().split(',')
        n, sample, bam_file = items
        cmd = '{0} -p 10 -o {1} -b {2} -u --library-type fr-unstranded {3} {4}'.format(cuffquant, sample, genome, merged_gtf, bam_file)
        print cmd
        os.system(cmd)

    cmd = 'find ./ -iname abundances.cxb > cxb_list.txt'
    print cmd
    os.system(cmd)
    cxbs = []
    for line in open('cxb_list.txt'):
        cxbs.append(line.strip())
    print cxbs
    samples = []
    for cxb in cxbs:
        samples.append(cxb.split('/')[1])

    cuffnorm = '/BiO/BioTools/cufflinks/cufflinks-2.2.1.Linux_x86_64/cuffnorm' 
    cmd = '{0} -p 10 -o cuffnorm_normForm -library-norm-method classic-fpkm --library-type fr-unstranded -L {1} {2} {3}'.format(cuffnorm, ','.join(samples), merged_gtf, ' '.join(cxbs))
    print cmd
    os.system(cmd)
    cmd = '{0} -p 10 -o cuffnorm_diffForm --output-format cuffdiff --library-norm-method classic-fpkm --library-type fr-unstranded -L {1} {2} {3}'.format(cuffnorm, ','.join(samples), merged_gtf, ' '.join(cxbs))
    print cmd
    os.system(cmd)

    rawCountDic = dict()
    for line in open('./cuffnorm_diffForm/genes.read_group_tracking'):
        if line.startswith('tracking_id'):
            continue
        items = line.strip().split('\t')
        tracking, sample, rep, raw_frags = items[:4]
        rawCountDic.setdefault(tracking, {}).setdefault(sample, raw_frags)
    out = open('./cuffnorm_diffForm/genes.raw_fragments', 'w')
    out.write('tracking_id\t{0}\n'.format('\t'.join(samples)))
    for tracking, samDic in sorted(rawCountDic.iteritems()):
        units = [tracking]
        for sample in samples:
            units.append(samDic[sample])
        out.write('{0}\n'.format('\t'.join(units)))
    out.close()


if __name__=='__main__':
    import os
    import subprocess
    import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--configuration-file', help='tophat output sam location information',
            default='configuration.csv')
    parser.add_argument('-g', '--gtf')
    parser.add_argument('-f', '--genomeFasta')
    args = parser.parse_args()
    main(args)
