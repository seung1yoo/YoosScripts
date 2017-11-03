import os

class BEDtools():
    def __init__(self, bedtools):
        self.bedtools = bedtools

    def bed_stats(self, abed):
        lociDic = dict()
        for line in open(abed):
            items = line.strip().split('\t')
            loci = '-'.join(items[0:3])
            lociDic.setdefault(loci, '')
        print '# uniqe loci of {0} : {1}'.format(abed, len(lociDic))

    def intersect(self, ref_bed, target_bed, sample):
        self.intersect_wab = '{0}.intersect.wab.bed'.format(sample)
        self.intersect_overlap = '{0}.intersect.overlap.bed'.format(sample)
        self.intersect_count = '{0}.intersect.count.bed'.format(sample)
        #### STUDY
        ### Reference http://bedtools.readthedocs.io/en/latest/content/example-usage.html
        ## Report the base-pair overlap between -a and -b.
        #cmd = '{0} intersect -a {1} -b {2} > {3}.intersect.align.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        #cmd = '{0} intersect -a {1} -b {2} -u > {3}.intersect.u.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        #cmd = '{0} intersect -a {1} -b {2} -wa > {3}.intersect.wa.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        #cmd = '{0} intersect -a {1} -b {2} -wb > {3}.intersect.wb.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        #cmd = '{0} intersect -a {1} -b {2} -wa -wb -sorted > {3}.intersect.wab.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        ## How many base pairs of overlap were there?
        #cmd = '{0} intersect -a {1} -b {2} -wo > {3}.intersect.overlap.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        ## Counting the number of overlapping features.
        #cmd = '{0} intersect -a {1} -b {2} -c > {3}.intersect.count.bed'.format(self.bedtools, ref_bed, target_bed, sample)
        #print cmd
        #os.system(cmd)
        #### 
        if not os.path.isfile(self.intersect_wab):
            cmd = '{0} intersect -a {1} -b {2} -wa -wb -sorted > {3}'.format(self.bedtools, ref_bed, target_bed, self.intersect_wab)
            print cmd
            os.system(cmd)
        if not os.path.isfile(self.intersect_overlap):
            cmd = '{0} intersect -a {1} -b {2} -wo > {3}'.format(self.bedtools, ref_bed, target_bed, self.intersect_overlap)
            print cmd
            os.system(cmd)
        if not os.path.isfile(self.intersect_count):
            cmd = '{0} intersect -a {1} -b {2} -c > {3}'.format(self.bedtools, ref_bed, target_bed, self.intersect_count)
            print cmd
            os.system(cmd)

    def intersect_wab_parser(self, sample):
        if not os.path.isfile(self.intersect_wab):
            import sys
            print 'intersect_wab file is missing'
            sys.exit()
        #
        map_Dic = dict()
        for line in open(self.intersect_wab):
            items = line.strip().split('\t')
            #
            a_loci = '_'.join(items[0:3])
            a_desc = items[3].split('(')[0].strip()
            #
            b_loci = '_'.join(items[12:15])
            b_desc = items[15]
            #
            map_Dic.setdefault(a_desc, {}).setdefault(a_loci, {}).setdefault(b_desc, b_loci)
        #
        out = open('{0}.intersect.table.xls'.format(sample), 'w')
        for a_desc, a_lociDic in map_Dic.iteritems():
            for a_loci, b_descDic in a_lociDic.iteritems():
                for b_desc, b_loci in b_descDic.iteritems():
                    out.write('{0}\n'.format('\t'.join([a_desc, a_loci, b_desc, b_loci])))
        out.close()
        #
        out = open('{0}.intersect.count.xls'.format(sample), 'w')
        a_desc_Dic = dict()
        a_loci_Dic = dict()
        b_Dic = dict()
        for a_desc, a_lociDic in map_Dic.iteritems():
            a_desc_Dic.setdefault(a_desc, '')
            for a_loci, b_descDic in a_lociDic.iteritems():
                a_loci_Dic.setdefault(a_loci, '')
                b_members = []
                for b_desc, b_loci in b_descDic.iteritems():
                    b_Dic.setdefault(b_desc, b_loci)
                    b_members.append('{0}@{1}'.format(b_desc, b_loci))
                out.write('{0}\n'.format('\t'.join([a_desc, a_loci, str(len(b_members)), ','.join(b_members)])))
        print '# aligned QTL : {0}'.format(len(a_desc_Dic))
        print '# aligned QTL loci : {0}'.format(len(a_loci_Dic))
        print '# aligned {0} miRNA : {1}'.format(sample, len(b_Dic))
        out.close()


def main(args):
    bt = BEDtools(args.bedtools)
    bt.bed_stats(args.ref_bed)
    bt.bed_stats(args.target_bed)
    bt.intersect(args.ref_bed, args.target_bed, args.sample)
    bt.bed_stats(bt.intersect_wab)
    bt.bed_stats(bt.intersect_overlap)
    bt.bed_stats(bt.intersect_count)
    bt.intersect_wab_parser(args.sample)

if __name__=='__main__':
    import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-bt', '--bedtools', default='bedtools')
    parser.add_argument('-rb', '--ref-bed',
            default='../Databases/cattle-QTLdb/Release33/cattle-QTLdb.Release33.Milk.All.sort.bed')
    parser.add_argument('-tb', '--target-bed',
            default='./utils/example.mirdeep.result.clean.bed')
    parser.add_argument('-s', '--sample',
            default='example')
    args = parser.parse_args()
    main(args)

