

def main():
    dirs = glob.glob('*')
    n = 0
    bams = []
    for dir in dirs:
        sam = '{0}/{0}.sam'.format(dir)
        if os.path.isdir(dir) and os.path.isfile(sam):
            n += 1
            cmd = 'samtools view -bT ../0.Database/TAIR10_chr_all.fas {0} > {0}.bam'.format(sam)
            bams.append('{0}.bam'.format(sam))
            #print n, cmd
            #os.system(cmd)

    covs = []
    for n, bam in enumerate(bams):
        ### **********!!!!!!!!!!!!!!!!!!!!!!!
        cmd = 'bedtools coverage -abam {0} -b ../0.Database/TAIR10.gtf.bed > {0}.TAIR10.bed.coverage'.format(bam)
        print n, cmd
        if os.path.isfile('{0}.TAIR10.bed.coverage'.format(bam)):
            covs.append('{0}.TAIR10.bed.coverage'.format(bam))
        else:
            os.system(cmd)
            


if __name__=='__main__':
    import os
    import glob
    main()
