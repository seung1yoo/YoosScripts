from Bio import SeqIO

class IDconverter():
    def __init__(self, infa):
        self.seqn = self.seqcounter(infa)

    def seqcounter(self, infa):
        self.seqn = 0
        for line in open(infa):
            if line.startswith('>'):
                self.seqn += 1
        self.seqdigit = len(str(self.seqn))
        #

    def convert(self, infa, outfa, outmap, prefix):
        out_outfa = open(outfa, 'w')
        out_outmap = open(outmap, 'w')
        n = 0
        for record in SeqIO.parse(open(infa), 'fasta'):
            n += 1
            #
            new_id = '{0}{1:0{2}}'.format(prefix, n, self.seqdigit)
            out_outmap.write('{0} ==> {1}\n'.format(record.description, new_id))
            #
            record.id = new_id
            record.description = ''
            #
            SeqIO.write(record, out_outfa, 'fasta')
            #
        out_outfa.close()
        out_outmap.close()


def main(args):
    print(args)
    #
    idc = IDconverter(args.infa)
    idc.convert(args.infa, args.outfa, args.outmap, args.prefix)
    print('#DONE')


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--infa', default='SequenceIDconverter.exam.in.fa')
    parser.add_argument('-of', '--outfa', default='SequenceIDconverter.exam.out.fa')
    parser.add_argument('-om', '--outmap', default='SequenceIDconverter.exam.map.fa')
    parser.add_argument('-p', '--prefix', default='exam')
    args = parser.parse_args()
    main(args)
