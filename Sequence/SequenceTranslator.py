
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(args):
    out = open(args.outPepFasta, 'w')
    for record in SeqIO.parse(open(args.inputNuclFasta), 'fasta'):
        new_record = SeqRecord(Seq(str(record.seq.translate())))
        new_record.id = record.id
        new_record.description = record.description
        SeqIO.write(new_record, out, 'fasta')
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-inf', '--inputNuclFasta', default='SequenceIDconverter.exam.in.fa')
    parser.add_argument('-opf', '--outPepFasta', default='SequenceIDconverter.exam.out.pep.fa')
    args = parser.parse_args()
    main(args)
