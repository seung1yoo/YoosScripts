

class GOSEQ_INPUT_MAKER:
    def __init__(self, args):
        self.cuffdiff = args.cuffdiff
        self.selectde = args.selectde
        self.genetogo = args.genetogo
        self.degset = args.degset
        self.outdir = args.outdir
        #
    def make_len_f(self):
        out = open('{0}/{1}.goseq.len.xls'.format(self.outdir, self.degset),'w')
        for line in open(self.cuffdiff):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['test_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            g_id = items[idx_dic['test_id']]
            g_locus = sorted(items[idx_dic['locus']].split(':')[1].split('-'))
            g_len = abs(int(g_locus[-1])-int(g_locus[0]))+1
            out.write('{0}\t{1}\n'.format(g_id, g_len))
        out.close()

    def make_deg_f(self):
        col_test = 3
        test_idx_dic = dict()
        #
        out = open('{0}/{1}.goseq.deg.xls'.format(self.outdir, self.degset),'w')
        for line in open(self.selectde):
            if line.startswith('# TEST'):
                col_test += 1
                test = line.rstrip('\n').split()[-1]
                test_idx_dic.setdefault(test, col_test)
                continue
            elif line.startswith('#'):
                continue
            #
            items = line.rstrip('\n').split('\t')
            g_id = items[0]
            infos = items[test_idx_dic[self.degset]].split(';')
            if infos[-1] in ['False']:
                out.write('{0}\t0\n'.format(g_id))
            elif infos[-1] in ['True']:
                out.write('{0}\t1\n'.format(g_id))
        out.close()

    def make_map_f(self):
        out = open('{0}/{1}.goseq.mapping.xls'.format(self.outdir, self.degset),'w')
        for line in open(self.genetogo):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['# Order']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            g_id = items[idx_dic['ID']]
            g_go = items[idx_dic['GO_ID_S']]
            out.write('{0}\t{1}\n'.format(g_id, g_go))
        out.close()


def main(args):
    gim = GOSEQ_INPUT_MAKER(args)
    gim.make_len_f()
    gim.make_deg_f()
    gim.make_map_f()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--cuffdiff', default="/BiO/BioProjects/TBD180874-PSU-Tomato-RNAref-20190114/Rine_Quant/analysis/diff-S0001-S0002.cuffdiff/gene_exp.diff")
    parser.add_argument('--selectde', default="/BiO/BioProjects/TBD180874-PSU-Tomato-RNAref-20190114/Rine_Quant/analysis/select_de.gene.xls")
    parser.add_argument('--genetogo', default='/BiO/BioProjects/TBD180874-PSU-Tomato-RNAref-20190114/Rine_Quant/output/genes_to_go_with_ancestor.xls')
    parser.add_argument('--degset', default='S0001-S0002')
    parser.add_argument('--outdir', default='/BiO/BioProjects/TBD180874-PSU-Tomato-RNAref-20190114/Rine_Quant/output/GO')
    args = parser.parse_args()
    main(args)
