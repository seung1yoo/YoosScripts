from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class SEQUENCE:
    def __init__(self):
        self.seqDic = dict()

    def get_list(self, listFn):
        for line in open(listFn):
            gene_id = line.rstrip('\n')
            if gene_id:
                self.seqDic.setdefault(gene_id, {})
            else:
                continue
            #
        #

    def get_seq(self, faFn, faKey):
        for record in SeqIO.parse(open(faFn), 'fasta'):
            gene_id = record.id.split('|')[1]
            if gene_id in self.seqDic:
                self.seqDic[gene_id].setdefault(faKey, record.seq)
            else:
                continue
            #
        #

    def make_fa(self, outPrefix, faKey):
        out = open('{0}.{1}.fa'.format(outPrefix, faKey),'w')
        for gene_id, faKeyDic in self.seqDic.iteritems():
            seq = faKeyDic[faKey]
            record = SeqRecord(seq, id=gene_id, description="")
            SeqIO.write(record, out, 'fasta')
        out.close()

    def make_table(self, outPrefix, *faKeys):
        out = open('{0}.table.xls'.format(outPrefix), 'w')
        out.write('gene_id\t{0}\n'.format('\t'.join(faKeys)))
        for gene_id, faKeyDic in self.seqDic.iteritems():
            items = [gene_id]
            items.extend([str(faKeyDic[x]) for x in faKeys])
            out.write('{0}\n'.format('\t'.join(items)))
        out.close()


def main():
    #
    listFn = 'list.txt'
    cdsfa = 'haliotis_0.1.CDS'
    pepfa = 'haliotis_0.1.protein'
    outPrefix = 'SelectedSeq'
    #
    seq = SEQUENCE()
    #
    seq.get_list(listFn)
    seq.get_seq(cdsfa, 'cds')
    seq.get_seq(pepfa, 'pep')
    #
    seq.make_fa(outPrefix, 'cds')
    seq.make_fa(outPrefix, 'pep')
    seq.make_table(outPrefix, 'cds', 'pep')



if __name__=='__main__':
    main()
