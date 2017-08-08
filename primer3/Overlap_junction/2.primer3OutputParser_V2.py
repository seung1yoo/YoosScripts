#!/usr/bin/python
#shin, younhee
#

import string, sys, os, re, time


class primer3Parser:
	
    def dicMaker (self, pos_file) :
	fd = open(pos_file, 'r')
	line = fd.readline()
	dicPos = {}
        while line != '' :
            unit = line.strip().split()
	    if len(unit) > 2 :
                dicPos[unit[0]] = unit[1]
            line = fd.readline()
        fd.close()
        return(dicPos)

    def iterResult(self, input_file):
        fr = open(input_file)
        lines = [fr.readline()]
        for line in fr:
            if not line.strip():
                continue
            elif line.startswith('PRIMER PICKING RESULTS'):
                yield lines
                lines = [line]
            else:
                lines.append(line)
        yield lines

    def parser (self, input_file, out_file, out_1_fasta, out_2_fasta, dicPos) :
        fw = open(out_file, 'w')
        fw2 = open(out_1_fasta, 'w')
        fw3 = open(out_2_fasta, 'w')
        for lines in self.iterResult(input_file):
            print lines
            seq_id = ''
            product_size = ''
            f_start = 0
            f_primer_inform = ''
            r_start = 0
            r_primer_inform = ''
            f_seq = ''
            r_seq = ''
            for line in lines:
                line = line.strip()
                if line.startswith("PRIMER PICKING RESULTS") :
                    seq_id = line.strip().split()[-1]
                elif line.startswith("LEFT PRIMER") :
                    f_primer_inform = '\t'.join(line.strip().split()[3:])
                    f_seq = line.strip().split()[-1]
                    f_start = int(dicPos[seq_id]) + int(line.strip().split()[2])
                elif line.startswith("RIGHT PRIMER") :
                    r_primer_inform = '\t'.join(line.strip().split()[3:])
                    r_seq = line.strip().split()[-1]
                    r_start = int(dicPos[seq_id]) + int(line.strip().split()[2])
                elif line.startswith("PRODUCT SIZE") :
                    product_size = line.strip().split()[2].split(',')[0]
                elif line.startswith("libprimer3 release 2.3.6") :
                    new_line = "%s\t%s\t%d\t%s\t%d\t%s\n"%(seq_id, product_size, f_start, f_primer_inform, r_start, r_primer_inform)
                    fw.write(new_line)
                    seq_line1 = '>%s/1\n%s\n'%(seq_id, f_seq)
                    seq_line2 = '>%s/2\n%s\n'%(seq_id, r_seq)
                    fw2.write(seq_line1)
                    fw3.write(seq_line2)
        fw.close()
        fw2.close()
        fw3.close()
        return(1)
            



#============ Parser ===================================================================
        
if __name__ == '__main__':
    executer = primer3Parser()
    input_file = sys.argv[1]
    pos_file = sys.argv[2]
    dicPos = executer.dicMaker(pos_file)
    out_file = sys.argv[1] + '.4excel'
    out_1_fasta = sys.argv[1] + '_1.fasta'
    out_2_fasta = sys.argv[1] + '_2.fasta'
    ans = executer.parser(input_file, out_file, out_1_fasta, out_2_fasta, dicPos)
    print "Your job was done"
