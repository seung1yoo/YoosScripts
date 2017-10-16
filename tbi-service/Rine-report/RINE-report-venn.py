

def main(args):
    print args

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='genes.de.xls',
        default='/BiO/BioProjects/TBD170758-gntech-Pig-RNAref-20171011/Rine_Quant/report/DEG/genes.de.xls')
    parser.add_argument('-o', '--output', help='name of output file',
        default='venndiagram_table.xls')
    args = parser.parse_args()
    main(args)

