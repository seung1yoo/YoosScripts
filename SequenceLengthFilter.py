def main(input, output, cutoff):
    for record in SeqIO.parse(input, 'fasta'):
        if len(record.seq) >= int(cutoff):
            SeqIO.write(record, output, 'fasta')

if __name__=='__main__':
    from Bio import SeqIO
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), help='input file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output file')
    parser.add_argument('-c', '--cutoff', type=str, help='cutoff value')
    args = parser.parse_args()
    main(args.input, args.output, args.cutoff)
