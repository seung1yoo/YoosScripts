from Bio import SeqIO

def seqSpliter(fasta, count):
    seqDic = dict()
    n = 0
    num_part = 1
    for record in SeqIO.parse(open(fasta), 'fasta'):
        if n <= count:
            n += 1
        else:
            num_part += 1
            n = 0
        seqDic.setdefault(num_part, {}).setdefault(n, record)
    return seqDic


def main(args):
    print('{0} split by {1}...'.format(args.fasta, args.count))
    seqDic = seqSpliter(args.fasta, args.count)

    for num_part, nDic in seqDic.items():
        print('writing part {0}'.format(str(num_part)))
        out = open('{0}.part{1}.fa'.format(args.outfasta_prefix, str(num_part)), 'w')
        for n, seq_record in nDic.items():
            SeqIO.write(seq_record, out, 'fasta')
        out.close()

    print('DONE!')


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='fasta file name')
    parser.add_argument('-c', '--count', type=int, help='number of seq count to split')
    parser.add_argument('-o', '--outfasta-prefix', help='prefix of out fasta file')
    args = parser.parse_args()
    main(args)


