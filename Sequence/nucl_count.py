

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def main(args):

    outfh = open(args.output, "w")
    for record in SeqIO.parse(open(args.input), 'fasta'):
        nucl_dic = dict()
        for nucl in str(record.seq):
            nucl_dic.setdefault(nucl, 0)
            nucl_dic[nucl] += 1
        print(record.id, sorted(nucl_dic))

        clean_seq = []
        for nucl in str(record.seq):
            if nucl in ['A','T','G','C']:
                clean_seq.append(nucl)
            else:
                clean_seq.append('N')
        nucl_dic = dict()
        for nucl in clean_seq:
            nucl_dic.setdefault(nucl, 0)
            nucl_dic[nucl] += 1
        print(record.id, sorted(nucl_dic))
        clean_record = SeqRecord(Seq(''.join(clean_seq)), id=record.id, description="")
        SeqIO.write(clean_record, outfh, 'fasta')
    outfh.close()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(args)
