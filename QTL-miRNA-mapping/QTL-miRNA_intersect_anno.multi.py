import os
import glob

def main():
    for afile in glob.glob('predicted-miRNAs/*.known.xls'):
        sample = afile.split('/')[-1].split('.')[0]

        #countFile = '0.Milk/{0}.intersect.count.xls'.format(sample)
        #tableFile = '0.Milk/{0}.intersect.table.xls'.format(sample)
        countFile = '1.Reproduction/{0}.intersect.count.xls'.format(sample)
        tableFile = '1.Reproduction/{0}.intersect.table.xls'.format(sample)
        annoKnown = 'predicted-miRNAs/{0}.mirdeep.result.known.xls'.format(sample)
        annoNovel = 'predicted-miRNAs/{0}.mirdeep.result.novel.xls'.format(sample)

        cmd = 'python QTL-miRNA_intersect_anno.py -cf {0} -tf {1} -ak {2} -an {3}'.format(countFile, tableFile, annoKnown, annoNovel)
        print cmd
        os.system(cmd)

if __name__=='__main__':
    main()
