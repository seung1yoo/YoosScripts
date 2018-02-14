
class GCS():
    def __init__(self):
        pass

    def headerList_parser(self, hlfn):
        headers = list()
        for line in open(hlfn):
            item = line.rstrip('\n')
            headers.append(item)
        return headers

    def select_col(self, igfn, headerID, headers):
        cols = list()
        for line in open(igfn):
            items = line.rstrip('\n').split('\t')
            if items[0] in [headerID]:
                for idx, item in enumerate(items):
                    if item in headers:
                        cols.append(idx)
                    else:
                        pass
            else:
                pass
        return cols

    def writeXls(self, igfn, cols, ogfn):
        out = open(ogfn, 'w')
        for line in open(igfn):
            items = line.rstrip('\n').split('\t')
            new_items = []
            for col in cols:
                new_items.append(str(items[col]))
            out.write('{0}\n'.format('\t'.join(new_items)))
        out.close()


def main(args):
    tor = GCS()
    headers = tor.headerList_parser(args.headerList)
    cols = tor.select_col(args.inputGenes, args.headerId, headers)
    tor.writeXls(args.inputGenes, cols, args.outfilename)
    print 'DONE'

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ig', '--inputGenes', help='genes.xls',
            default='genes.xls')
    parser.add_argument('-hl', '--headerList',
            default='genes.header')
    parser.add_argument('-hi', '--headerId',
            default='Order')
    parser.add_argument('-ofn', '--outfilename',
            default='genes.outfilename.xls')
    args = parser.parse_args()
    main(args)
