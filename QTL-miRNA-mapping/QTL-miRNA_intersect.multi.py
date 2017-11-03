import os
import glob

def main():
    for abed in glob.glob('predicted-miRNAs/*.clean.bed'):
        sample = abed.split('/')[-1].split('.')[0]
        ## MILK
        cmd = 'python QTL-miRNA_intersect.py -rb cattle-QTLdb/cattle-QTLdb.Release33.Milk.All.sort.bed -tb {0} -s {1} > {1}.intersect.log'.format(abed, sample)
        print cmd
        os.system(cmd)
        ## Reproduct
        #cmd = 'python QTL-miRNA_intersect.py -rb cattle-QTLdb/cattle-QTLdb.Release33.Reproduction.All.sort.bed -tb {0} -s {1} > {1}.intersect.log'.format(abed, sample)
        #print cmd
        #os.system(cmd)

if __name__=='__main__':
    main()
