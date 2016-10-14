
def fileSplit(file, cpus):
    count = 0
    for record in SeqIO.parse(open(file), 'fasta'):
        count += 1
    unit = round(count/float(cpus), 0)
    
    tag = 1
    count = 0
    outfile = '{0}.{1}.fna'.format(file, tag)
    outhandle = open(outfile, 'w')
    files = []
    for record in SeqIO.parse(open(file), 'fasta'):
        count += 1
        if count <= unit:
            SeqIO.write(record, outhandle, 'fasta')
            pass
        else:
            files.append(outfile)
            outhandle.close()
            tag += 1
            outfile = '{0}.{1}.fna'.format(file, tag)
            outhandle = open(outfile, 'w')
            count = 1
            SeqIO.write(record, outhandle, 'fasta')
    files.append(outfile)
    outhandle.close()
    return files

def makeCmd(files, database, program, evalue):
    cmds = []
    for file in files:
	## BLASTx command version
        if program in ['blastx']: ## not validate
            cmd = '{0} -query {1} -db {2} -out {1}_vs_{4}.{3}.xml -evalue {3} -outfmt 5 -show_gis &'.format(
                    program, file, database, evalue, database.split('/')[-1])
        elif program in ['blastn']:
            cmd = '{0} -query {1} -db {2} -out {1}_vs_{4}.{3}.xml -evalue {3} -outfmt 5 -show_gis &'.format(
                    program, file, database, evalue, database.split('/')[-1])
        elif program in ['blastp']:
            cmd = '{0} -query {1} -db {2} -out {1}_vs_{4}.{3}.xml -evalue {3} -outfmt 5 -show_gis &'.format(
                    program, file, database, evalue, database.split('/')[-1])
	## MEGABLAST command version
        elif program in ['megablast']:
            cmd = '{0} -i {1} -d {2} -o {1}_vs_{4}.{3}.xml -e {3} -m 7 &'.format(
                    program, file, database, evalue, database.split('/')[-1])
        else:
            print '# Check the program'
            sys.exit()
        cmds.append(cmd)
    return cmds

def execute(cmds):
    for cmd in cmds:
        print cmd
        os.system(cmd)

def main(args):
    print args
    files = fileSplit(args.query, args.cpus)
    cmds = makeCmd(files, args.database, args.program, args.evalue)
    execute(cmds)

if __name__=='__main__':
    import sys
    import os
    from Bio import SeqIO
    import argparse
    parser = argparse.ArgumentParser('BLAST executor')
    parser.add_argument('-q', '--query', default='Trinity_CD_Hit_500_min_Result.fasta.rename')
    parser.add_argument('-d', '--database', default='/tier1/Codes/database/nt/nt')
    parser.add_argument('-p', '--program', default='megablast')
    parser.add_argument('-e', '--evalue', default='1e-1')
    parser.add_argument('-c', '--cpus', default=8)
    args = parser.parse_args()
    main(args)
