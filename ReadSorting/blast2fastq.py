
from Bio import SeqIO
import gzip

class Blast2Fastq:
    def __init__(self, args):
        self.read_dic = dict()
        self.evalue = 1e-5
        self.readcov = 50.0
        self.alignsim = 70.0

        self.add_target_read(args.blast1, 'R1')
        self.add_target_read(args.blast2, 'R2')
        self.write_target_read(args.rawfq1, args.sample, 'R1')
        self.write_target_read(args.rawfq2, args.sample, 'R2')

    def add_target_read(self, blast, r_tag):
        for line in open(blast):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['queryNum']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            read_id = items[idx_dic['queryID']].split()[0]
            evalue = float(items[idx_dic['e-value']])
            readcov = float(items[idx_dic['queryCov']])
            alignsim = float(items[idx_dic['alignSim']])
            #
            if evalue > self.evalue:
                continue
            if readcov < self.readcov:
                continue
            if alignsim < self.alignsim:
                continue
            self.read_dic.setdefault(read_id, {}).setdefault(r_tag, None)
        return 1

    def write_target_read(self, raw, sample, r):
        outfh = open(f'{sample}_{r}.blast2.fastq', 'w')
        for record in SeqIO.parse(gzip.open(raw, 'rt'), 'fastq'):
            if record.id in self.read_dic:
                if 'R1' in self.read_dic[record.id] and 'R2' in self.read_dic[record.id]:
                    SeqIO.write(record, outfh, 'fastq')
        outfh.close()



def main(args):
    print(args)

    urs = Blast2Fastq(args)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1', '--blast1', default='/data07/project/TBD210878_13054_NIBR_Moss_Denovo_Assenbly_sjchoi_20210917/blast/ViviJP/CP_ViviJP_1.merge.1e-5.top1.tsv')
    parser.add_argument('-b2', '--blast2', default='/data07/project/TBD210878_13054_NIBR_Moss_Denovo_Assenbly_sjchoi_20210917/blast/ViviJP/CP_ViviJP_2.merge.1e-5.top1.tsv')
    parser.add_argument('-r1', '--rawfq1', default='/data07/project/TBD210878_13054_NIBR_Moss_Denovo_Assenbly_sjchoi_20210917/clean/ViviJP_1.fq.gz')
    parser.add_argument('-r2', '--rawfq2', default='/data07/project/TBD210878_13054_NIBR_Moss_Denovo_Assenbly_sjchoi_20210917/clean/ViviJP_2.fq.gz')
    parser.add_argument('-s', '--sample', default='ViviJP')
    args = parser.parse_args()
    main(args)

