



class DAVID:
    def __init__(self, david_fn, anno_cate):
        self.david_fn = david_fn
        self.anno_cate = anno_cate

    def parse_david_anno_table(self):
        #
        self.gene_dic = dict()
        self.anno_dic = dict()
        #
        for line in open(self.david_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['ID']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            if len(items) <= idx_dic[self.anno_cate]:
                continue
            _id_s = [x.strip() for x in items[idx_dic['ID']].split(',')]
            _anno_s = self.parse_anno_string(items[idx_dic[self.anno_cate]])
            for _id in _id_s:
                for _anno in _anno_s:
                    if not _anno:
                        continue
                    anno_id, anno_desc = _anno.split(':')
                    self.anno_dic.setdefault(anno_id, anno_desc)
                    self.gene_dic.setdefault(_id, {}).setdefault(anno_id, anno_desc)
        #
        print('annotated gene ==> {0}'.format(len(self.gene_dic)))
        print('annotation ==> {0}'.format(len(self.anno_dic)))

    def parse_anno_string(self, string):
        _anno_s = list()
        pre_sub = ''
        for sub in string.split(','):
            if not pre_sub and ':' in sub:
                pre_sub = sub
            elif pre_sub and ':' in sub:
                _anno_s.append(pre_sub)
                pre_sub = sub
            elif not sub:
                _anno_s.append(pre_sub)
                pre_sub = sub
            elif pre_sub and not ':' in sub:
                pre_sub = '{0},{1}'.format(pre_sub, sub)
            else:
                print('ERROR', pre_sub, sub)
                import sys
                sys.exit()
        return _anno_s

    def write_david_anno_table(self):
        out_fh = open('{0}.table'.format(self.david_fn), 'w')
        for anno_id, anno_desc in self.anno_dic.items():
            out_fh.write('{0}\t{1}\t{2}\n'.format(anno_id, self.anno_cate, anno_desc))
        out_fh.close()

    def write_anno_genes(self, genes_fn, out_fn):
        out_fh = open(out_fn, 'w')
        for line in open(genes_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['GeneId']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                #
                out_fh.write('{0}\t{1}\n'.format('\t'.join(items), self.anno_cate))
                continue
            _id = items[idx_dic['GeneAcc']]
            if _id in self.gene_dic:
                items.append(','.join(sorted(self.gene_dic[_id].keys())))
            else:
                items.append('-')
            out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()








def main(args):
    david = DAVID(args.david_fn, args.anno_cate)
    david.parse_david_anno_table()
    david.write_david_anno_table()
    david.write_anno_genes(args.genes_fn, args.out_fn)





if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--david-fn',
            default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180681/GSEA_20190409/tr_C0A76205A3D21554775920385.txt')
    parser.add_argument('--anno-cate',
            default='KEGG_PATHWAY')
    parser.add_argument('--genes-fn',
            default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180681/GSEA_20190409/gene.de.xls')
    parser.add_argument('--out-fn',
            default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180681/GSEA_20190409/gene.de.kegg.xls')
    args = parser.parse_args()
    main(args)
