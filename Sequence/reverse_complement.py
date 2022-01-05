

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

outfh = open('Circularized_assembly_1_Dep_vir_01.RC.fasta', 'w')
for record in SeqIO.parse(open("Circularized_assembly_1_Dep_vir_01.fasta"), 'fasta'):
    if record.id in ['Contig1', 'Contig6', 'Contig4']:
        rc_seq = record.seq.reverse_complement()
        rc_record = SeqRecord(rc_seq, id="{0}_rc".format(record.id), description="")
        SeqIO.write(rc_record, outfh, 'fasta')
outfh.close()


