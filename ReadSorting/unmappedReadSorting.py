
from Bio import SeqIO
import gzip

class UnmappedReadSorting:
    def __init__(self, args):
        self.read_dic = dict()
        self.r1_name = 'R1'
        self.r2_name = 'R2'
        self.add_target_read(args.unmap1, self.r1_name)
        self.add_target_read(args.unmap2, self.r2_name)
        self.write_target_read(args.raw1, args.sample, self.r1_name)
        self.write_target_read(args.raw2, args.sample, self.r2_name)

    def add_target_read(self, unmap, r_name):
        print(unmap, r_name)
        for record in SeqIO.parse(unmap, 'fastq'):
            self.read_dic.setdefault(record.id, {}).setdefault(r_name, None)
        return 1

    def write_target_read(self, raw, sample, r):
        outfh = open(f'{sample}_{r}.fastq', 'w')
        for record in SeqIO.parse(gzip.open(raw, 'rt'), 'fastq'):
            if record.id in self.read_dic:
                if self.r1_name in self.read_dic[record.id] and self.r2_name in self.read_dic[record.id]:
                    SeqIO.write(record, outfh, 'fastq')
        outfh.close()



def main(args):
    print(args)

    urs = UnmappedReadSorting(args)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-u1', '--unmap1', default='/data06/project/sjchoi/TBD210970_13156_Mouse_RNAref_20210930/Mouse/analysis/star/S_352/Unmapped.out.mate1')
    parser.add_argument('-u2', '--unmap2', default='/data06/project/sjchoi/TBD210970_13156_Mouse_RNAref_20210930/Mouse/analysis/star/S_352/Unmapped.out.mate2')
    parser.add_argument('-r1', '--raw1', default='/data06/project/sjchoi/TBD210970_13156_Mouse_RNAref_20210930/rawdata_merge/S_352_R1.fastq.gz')
    parser.add_argument('-r2', '--raw2', default='/data06/project/sjchoi/TBD210970_13156_Mouse_RNAref_20210930/rawdata_merge/S_352_R2.fastq.gz')
    parser.add_argument('-s', '--sample', default='S_352')
    args = parser.parse_args()
    main(args)

