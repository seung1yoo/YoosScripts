from Bio import SeqIO

def listDicMaker(idListFile):
    idDic = dict()
    for line in open(idListFile):
        items = line.rstrip('\n').split('\t')
        idDic.setdefault(items[0])
    return idDic

def main(args):
    print args
    idDic = listDicMaker(args.idList)

    out = open(args.output)
    for record in SeqIO.parse(open(args.input), args.inputType):
        if record.id in idDic:
            SeqIO.write(out, record, args.inputType)
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-l', '--idList')
    parser.add_argument('-t', '--inputType', choices=('fastq','fasta'))
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args)
