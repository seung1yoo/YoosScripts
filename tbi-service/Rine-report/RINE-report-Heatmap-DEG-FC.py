def idxFinder(items):
    idxDic = dict()
    for idx, item in enumerate(items):
        if item in ['GeneId']:
            idxDic.setdefault('id', idx)
        elif item.startswith('DEG:') and item.endswith(':SELECT'):
            idxDic.setdefault('sel', []).append(idx)
        elif item.startswith('DEG:') and item.endswith(':Log2FC'):
            idxDic.setdefault('fc', []).append(idx)
    return idxDic

def main(args):
    print(args)

    out = open(args.outxls, 'w')
    for line in open(args.xls):
        items = line.rstrip('\n').split('\t')
        new_items = list()
        if line.startswith('Order'):
            idxDic = idxFinder(items)
            print(idxDic)
            new_items.append(items[idxDic['id']])
            new_items.extend([items[idx] for idx in idxDic['fc']])
            #print(new_items)
            out.write('{0}\n'.format('\t'.join(new_items)))
            continue
        sels = [items[idx] for idx in idxDic['sel']]
        if not 'Y' in sels:
            print('FILTER-OUT : {0} : {1}'.format(items[idxDic['id']], sels))
            continue
        fcs = [items[idx] for idx in idxDic['fc']]
        fcs = ['0.0' if value in ['-'] else value for value in fcs]
        new_items.append(items[idxDic['id']])
        new_items.extend(fcs)
        #print(new_items)
        out.write('{0}\n'.format('\t'.join(new_items)))
    out.close()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xls', help='RINE report genes.xls file')
    parser.add_argument('-o', '--outxls', help='outfile name')
    args = parser.parse_args()
    main(args)
