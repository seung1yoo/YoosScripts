import csv

def main(args):

    out = open(args.xls, 'w')
    for items in csv.reader(open(args.csv)):
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv')
    parser.add_argument('--xls')
    args = parser.parse_args()
    main(args)
