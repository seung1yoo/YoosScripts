

def selector(criteria, threshold, items, idxDic):
    targets = []
    for col, item in idxDic[criteria].items():
        targets.append(float(items[col]))
    filtereds = []
    for target in targets:
        if criteria in ['Log2FC']:
            if abs(target) >= threshold:
                filtereds.append(target)
        elif criteria in ['PV', 'QV']:
            if target <= threshold:
                filtereds.append(target)
    return targets, filtereds

def main(args):
    print(args)

    out = open(args.out_xls, 'w')
    for line in open(args.xls):
        select_tags = []
        items = line.rstrip('\n').split('\t')
        if items[0] in ['Order']:
            out.write(line)
            idxDic = dict()
            for idx, item in enumerate(items):
                if item.endswith(':Log2FC'):
                    idxDic.setdefault('Log2FC', {}).setdefault(idx, item)
                if item.endswith(':PV'):
                    idxDic.setdefault('PV', {}).setdefault(idx, item)
                if item.endswith(':QV'):
                    idxDic.setdefault('QV', {}).setdefault(idx, item)
            #CHECK-idxDic
            for criteria, colDic in idxDic.items():
                for col, item in colDic.items():
                    print('#CHECK-idxDic =', criteria, col, item.split(':')[-1])
                    #print('#CHECK-idxDic =', criteria, col, item)
                    pass
            continue
        ## VER 1. Selected DEGs must be in all DEG sets
        ##selector by "Log2FC"
        #targets, filtereds = selector('Log2FC', args.foldChange, items, idxDic)
        #if len(targets) in [len(filtereds)]:
        #    select_tags.append('Log2FC')
        ##selector by "PV"
        #targets, filtereds = selector('PV', args.p_value, items, idxDic)
        #if len(targets) in [len(filtereds)]:
        #    select_tags.append('PV')
        ##selector by "QV"
        #targets, filtereds = selector('QV', args.q_value, items, idxDic)
        #if len(targets) in [len(filtereds)]:
        #    select_tags.append('QV')
        #if select_tags in [['Log2FC', 'PV', 'QV']]:
        #    out.write(line)
        #    print(items[0:2], select_tags)

        ## VER 2. Selected DEGs must be in at least one DEG set
        degDic = dict()
        for criteria, colDic in idxDic.items():
            for col, item in colDic.items():
                deg_tag = item.split(':')[1]
                degDic.setdefault(deg_tag, {}).setdefault(criteria, float(items[col]))

        selected_degs = []
        for deg_tag, criteriaDic in degDic.items():
            #print(deg_tag, criteriaDic['Log2FC'])
            #print(deg_tag, criteriaDic['PV'])
            #print(deg_tag, criteriaDic['QV'])
            if abs(criteriaDic['Log2FC']) >= args.foldChange and criteriaDic['PV'] <= args.p_value and criteriaDic['QV'] <= args.q_value:
                selected_degs.append(deg_tag)

        if selected_degs:
            out.write(line)
    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xls', help='xls file as Rine report format')
    parser.add_argument('-fc', '--foldChange', type=float, default=1.0, help='log2foldChange value for select.')
    parser.add_argument('-p', '--p-value', type=float, default=1.0)
    parser.add_argument('-q', '--q-value', type=float, default=1.0)
    parser.add_argument('-o', '--out-xls', help='out xls file name')
    args = parser.parse_args()
    main(args)
