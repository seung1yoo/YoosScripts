import os


def annoFileParser(annoFile):
    annoDic = dict()
    for line in open(annoFile):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['mature miRBase miRNAs detected by miRDeep2', 'novel miRNAs predicted by miRDeep2']:
            print '# annoFileParser : {0}'.format(items[0])
            continue
        if items[0] in ['tag id', 'provisional id']:
            idxDic = dict()
            for idx, item in enumerate(items):
                idxDic.setdefault(idx, item)
                print '# index col.names : {0} {1}'.format(idx, item)
            annoColnames = items
            continue
        _id = items[0]
        for idx, colname in idxDic.iteritems():
            annoDic.setdefault(_id, {}).setdefault(colname, items[idx])
    return annoDic, annoColnames

def annotator_countFile(countFile, knownDic, novelDic):
    out_fn = '{0}.anno{1}'.format(os.path.splitext(countFile)[0], os.path.splitext(countFile)[1])
    out = open(out_fn, 'w')
    new_cols = ['#QTL', 'Locus', 'Count', 'Members']
    out.write('{0}\n'.format('\t'.join(new_cols)))
    for line in open(countFile):
        items = line.rstrip('\n').split('\t')
        qtl = items[0]
        locus = items[1]
        count = items[2]
        new_items = [qtl, locus, count]
        #
        members = items[3].split(',')
        new_members = []
        for member in members:
            _mirna = member.split('@')[0]
            mi_locus = member.split('@')[1]
            tag = _mirna.split(':')[0]
            mirna = _mirna.split(':')[1]
            #
            if tag in ['known']:
                new_members.append('{0}:{1}@{2}'.format(tag,
                    knownDic[mirna]['mature miRBase miRNA'], mi_locus))
            elif tag in ['novel']:
                new_members.append(member)
        new_items.append(','.join(new_members))
        out.write('{0}\n'.format('\t'.join(new_items)))
    out.close()

def annotator_tableFile(tableFile, knownDic, novelDic, knownCols, novelCols):
    out_fn = '{0}.anno{1}'.format(os.path.splitext(tableFile)[0], os.path.splitext(tableFile)[1])
    out = open(out_fn, 'w')
    new_cols = ['#QTL', 'QTL_locus', 'predicted_miRNA', 'predicted_miRNA_locus']
    new_cols.extend(novelCols[1:])
    out.write('{0}\n'.format('\t'.join(new_cols)))
    for line in open(tableFile):
        items = line.rstrip('\n').split('\t')
        qtl = items[0]
        qtl_locus = items[1]
        mirna_tag = items[2].split(':')[0]
        mirna_id = items[2].split(':')[1]
        mirna_locus = items[3]
        #
        if mirna_tag in ['known']:
            for col in knownCols[1:]:
                items.append(knownDic[mirna_id][col])
        elif mirna_tag in ['novel']:
            for col in novelCols[1:]:
                items.append(novelDic[mirna_id][col])
        #
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()


def main(args):
    print args
    knownDic, knownCols = annoFileParser(args.annoKnown)
    novelDic, novelCols = annoFileParser(args.annoNovel)
    print '# no. known miRNA : {0}'.format(len(knownDic))
    print '# no. novel miRNA : {0}'.format(len(novelDic))
    #
    annotator_countFile(args.countFile, knownDic, novelDic)
    annotator_tableFile(args.tableFile, knownDic, novelDic, knownCols, novelCols)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-cf', '--countFile', help='count.xls file of QTL-miRNA_intersect.py results',
            default='example.intersect.count.xls')
    parser.add_argument('-tf', '--tableFile', help='table.xls file of QTL-miRNA_intersect.py results',
            default='example.intersect.table.xls')
    parser.add_argument('-ak', '--annoKnown', help='known.xls file of predicted-miRNA-csv_cleaner.py results',
            default='utils/example.mirdeep.result.known.xls')
    parser.add_argument('-an', '--annoNovel', help='novel.xls file of predicted-miRNA-csv_cleaner.py results',
            default='utils/example.mirdeep.result.novel.xls')
    args = parser.parse_args()
    main(args)
