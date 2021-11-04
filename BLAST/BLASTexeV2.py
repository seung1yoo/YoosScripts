
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_num_fasta(fa_fn):
    count = 0
    for record in SeqIO.parse(open(fa_fn), 'fasta'):
        count += 1
    return count

def split_fasta(fa_fn, cpus, tmp_dir='./tmp'):
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    expCountPerSubFa = round(get_num_fasta(fa_fn)/float(cpus), 0)
    obsCountPerSubFa = 0
    numSubFa = 1

    splitted_fn_s = []
    splitted_fn = os.path.join(tmp_dir, 'splittedfasta_{0:0>3}.fa'.format(numSubFa))
    splitted_fh = open(splitted_fn, 'w')

    for record in SeqIO.parse(open(fa_fn), 'fasta'):
        obsCountPerSubFa += 1
        if obsCountPerSubFa <= expCountPerSubFa:
            SeqIO.write(record, splitted_fh, 'fasta')

        else:
            splitted_fn_s.append(splitted_fn)
            splitted_fh.close()

            obsCountPerSubFa = 1
            numSubFa += 1

            splitted_fn = os.path.join(tmp_dir, 'splittedfasta_{0:0>3}.fa'.format(numSubFa))
            splitted_fh = open(splitted_fn, 'w')

            SeqIO.write(record, splitted_fh, 'fasta')

    splitted_fn_s.append(splitted_fn)
    splitted_fh.close()

    return splitted_fn_s

def get_program_path(program):
    if program in ['blastx']:
        return 'blastx'
    elif program in ['blastn']:
        return 'blastn'
    elif program in ['blastp']:
        return 'blastp'
    elif program in ['megablast']:
        return 'megablast'

def makeCmd(splitted_fn_s, database, program, evalue, prefix):
    cmds = []
    for idx, splitted_fn in enumerate(splitted_fn_s):
        if program in ['blastx']:
            cmd = '{0} -query {1} -db {2} -out {4}.{5:0>3}.{3}.xml -evalue {3} -outfmt 5 -show_gis'.format(
                    get_program_path(program), splitted_fn, database, evalue, prefix, idx+1)
        elif program in ['blastn']:
            cmd = '{0} -query {1} -db {2} -out {4}.{5:0>3}.{3}.xml -evalue {3} -outfmt 5 -show_gis'.format(
                    get_program_path(program), splitted_fn, database, evalue, prefix, idx+1)
        elif program in ['blastp']:
            cmd = '{0} -query {1} -db {2} -out {4}.{5:0>3}.{3}.xml -evalue {3} -outfmt 5 -show_gis'.format(
                    get_program_path(program), splitted_fn, database, evalue, prefix, idx+1)
        elif program in ['megablast']:
            cmd = '{0} -i {1} -d {2} -o {4}.{5:0>3}.{3}.xml -e {3} -m 7'.format(
                    get_program_path(program), splitted_fn, database, evalue, prefix, idx+1)
        else:
            print('#ERROR : Check the program path.')
            sys.exit()
        cmds.append(cmd)
    return cmds

def execute(cmds, prefix):
    shell_fn = '{0}.pipe.sh'.format(prefix)
    shell_fh = open(shell_fn, 'w')
    shell_fh.write('#fork\n')
    for cmd in cmds:
        shell_fh.write('\t{0}\n'.format(cmd))
    shell_fh.write('#join\n')
    shell_fh.close()

def main(args):
    splitted_fn_s = split_fasta(args.query, args.cpus, os.path.join(args.prefix, 'tmp'))
    cmds = makeCmd(splitted_fn_s, args.database, args.program, args.evalue, args.prefix)
    execute(cmds, args.prefix)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser('BLAST executor')
    parser.add_argument('query', help='input file name as fasta format')
    parser.add_argument('database', help='blast database prefix')
    parser.add_argument('program', choices=('blastx','blastn','blastp','megablast'))
    parser.add_argument('evalue', help='1e-5')
    parser.add_argument('prefix', help='blast outfile prefix')
    parser.add_argument('cpus', help='minimum 2')
    args = parser.parse_args()
    main(args)
