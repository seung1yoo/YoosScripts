from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(args):
    print args
    out = open('{0}.{1}_{2}-{3}.fa'.format(args.prefix, args.seq_id, args.seq_start, args.seq_end), 'w')
    for record in SeqIO.parse(open(args.in_fa_fn), 'fasta'):
        if record.id in [args.seq_id]:
            new_Seq = record.seq[args.seq_start-1:args.seq_end]
            new_id = '{0}_{1}-{2}'.format(args.seq_id, args.seq_start, args.seq_end)
            new_record = SeqRecord(new_Seq, id=new_id, description="")
            SeqIO.write(new_record, out, 'fasta')
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_fa_fn', help='FASTA file name')
    parser.add_argument('--seq_id', help='Extraction Target ID')
    parser.add_argument('--seq_start', type=int, help='Extraction Target start')
    parser.add_argument('--seq_end', type=int, help='Extraction Target end')
    parser.add_argument('--prefix', help='Extraction OUTPUT FASTA file prefix')
    args = parser.parse_args()
    main(args)
