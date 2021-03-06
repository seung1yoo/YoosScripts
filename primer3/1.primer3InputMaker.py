
def _seqDicMaker(reference, basketDic, outfile):
    trWithSsr = dict()
    for key, unitDic in basketDic.iteritems():
        trWithSsr.setdefault(unitDic['trID'], '')

    seqDic = dict()
    out = open('{0}.info'.format(outfile), 'w')
    for record in SeqIO.parse(open(reference), 'fasta'):
        if trWithSsr.has_key(record.id):
            seqDic.setdefault(record.id, str(record.seq))
            out.write('{0}\t0\t{1}\n'.format(record.id, len(str(record.seq))))
    out.close()
    return seqDic

def seqDicMaker(reference, basketDic, outfile):
    trWithSsr = dict()
    for key, unitDic in basketDic.iteritems():
        trWithSsr.setdefault(unitDic['trID'], '')

    seqDic = dict()
    for record in SeqIO.parse(open(reference), 'fasta'):
        if trWithSsr.has_key(record.id):
            seqDic.setdefault(record.id, str(record.seq))
    return seqDic

def makeBasket(infile):
    basketDic = dict()
    for n, line in enumerate(open(infile)):
        items = line.rstrip('\n').split('\t')
        key = '{0}'.format('-'.join([items[0],items[1],items[2]]))
        basketDic.setdefault(key, {}).setdefault('trID', items[0])
        basketDic.setdefault(key, {}).setdefault('start', items[1])
        basketDic.setdefault(key, {}).setdefault('end', items[2])
        basketDic.setdefault(key, {}).setdefault('infos', items[3:])
    print 'make basket total : {0}'.format(n)
    return basketDic

def makePRIMER3_input(seqDic, basketDic, amplicon_size, outfile, config):
    out = open('{0}.{1}'.format(outfile,amplicon_size), 'w')
    outInfo = open('{0}.info'.format(outfile), 'w')
    for key, unitDic in basketDic.iteritems():
        length = int(unitDic['end']) - int(unitDic['start']) + 1
        info = '''\
SEQUENCE_ID={0}
SEQUENCE_TEMPLATE={1}
SEQUENCE_TARGET={2},{3}
PRIMER_TASK=generic
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=55.0
PRIMER_MAX_TM=65.0
PRIMER_OPT_GC_PERCENT=50.0
PRIMER_MIN_GC=40.0
PRIMER_MAX_GC=60.0
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PRODUCT_SIZE_RANGE={4}
P3_FILE_FLAG=1
PRIMER_NUM_RETURN=7
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH={5}
=
'''.format(
        '-'.join([unitDic['trID'], unitDic['start'], unitDic['end'], '-'.join(unitDic['infos'])]),
        seqDic[unitDic['trID']][int(unitDic['start'])-500:int(unitDic['end'])+500],
        str(499),
        length,
        amplicon_size,
        config)
        out.write(info)
        outInfo.write('{0}\t{1}\t{2}\n'.format(
            '-'.join([unitDic['trID'], unitDic['start'], unitDic['end'], '-'.join(unitDic['infos'])]),
            str(int(unitDic['start'])-499),
            str(int(unitDic['end'])+499)))
    out.close()
    outInfo.close()

def executePrimer3(infile, outfile, amplicon_size):
    cmd = 'primer3_core -format_output -output={1}.{2} {0}.{2}'.format(infile, outfile, amplicon_size)
    print cmd
    os.system(cmd)

def main(args):
    basketDic = makeBasket(args.table)
    seqDic = seqDicMaker(args.reference, basketDic, args.outfile)
    makePRIMER3_input(seqDic, basketDic, args.amplicon_size, args.infile, args.config)
    executePrimer3(args.infile, args.outfile, args.amplicon_size)

if __name__=='__main__':
    from Bio import SeqIO
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--table',
            help = 'chr\tstart\tend\tsample\trefSample\tref\talt\t[whatever...]',
            default='snpEffxls_primerSet.xls')
    parser.add_argument('-r', '--reference',
            default='genome.scf.fasta')
    parser.add_argument('-s', '--amplicon-size',
            default='100-300')
    parser.add_argument('-i', '--infile',
            default='primer3.input')
    parser.add_argument('-o', '--outfile',
            default='primer3.output')
    parser.add_argument('-c', '--config')
    args = parser.parse_args()
    main(args)

#PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/data/Bioinformatics/Tools/src/primer3-2.3.4/src/primer3_config/
