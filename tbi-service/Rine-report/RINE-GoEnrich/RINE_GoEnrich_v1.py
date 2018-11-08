#!/usr/bin/python
# @seung1.yoo


class RINE_SAMPLE_IN_PARSER:
    def mk_sam_label_dic(self, sam_in_fn):
        s_l_dic = dict()
        for line in open(sam_in_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split()
            if items[0] in ['INPUT']:
                sample_id = 'S{0:0>4}'.format(items[1])
                label = items[-1]
                s_l_dic.setdefault(sample_id, label)
        return s_l_dic

    def mk_label_sam_dic(self, sam_in_fn):
        l_s_dic = dict()
        for line in open(sam_in_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split()
            if items[0] in ['INPUT']:
                sample_id = 'S{0:0>4}'.format(items[1])
                label = items[-1]
                l_s_dic.setdefault(label, sample_id)
        return l_s_dic


class RINE_GENES_PARSER:
    headerkey = 'Order'
    tbiu_colname = 'GeneId'
    go_colname = 'GeneOntology'
    cutoff = 'PV=0.05'

    def mk_deg_dic(self, genes_fn):
        deg_dic = dict()
        for line in open(genes_fn):
            items = line.rstrip('\n').split('\t')

            if items[0] in [RINE_GENES_PARSER.headerkey]:
                deg_design_dic = self.ext_deg_design(items)
                #
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue

            tbiu = items[idx_dic[RINE_GENES_PARSER.tbiu_colname]]
            go = items[idx_dic[RINE_GENES_PARSER.go_colname]]

            for design_pair in deg_design_dic.keys():
                deg_dic.setdefault(design_pair,{})
                #
                is_deg = items[idx_dic['DEG:{0}:{1}:SELECT'.format(design_pair, RINE_GENES_PARSER.cutoff)]]
                if is_deg in ['Y']:
                    deg_dic[design_pair].setdefault(tbiu, go)

        return deg_dic

    def ext_deg_design(self, items):
        deg_design_dic = dict()
        for item in items:
            if item.startswith('DEG:'):
                design_pair = item.split(':')[1]
                pair_s = design_pair.split(' vs ')
                ctrl_s = pair_s[0].split(',')
                case_s = pair_s[1].split(',')
                deg_design_dic.setdefault(design_pair,{}).setdefault('ctrl',ctrl_s)
                deg_design_dic.setdefault(design_pair,{}).setdefault('case',case_s)
            else:
                continue
        return deg_design_dic


class RINE_GO2GENES_PARSER:
    def mk_go_d_dic(self, fn):
        go_d_dic = dict()
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['GO_ID']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            go_d_dic.setdefault(items[idx_dic['GO_ID']], {})
            go_d_dic[items[idx_dic['GO_ID']]].setdefault('c',items[idx_dic['Category']])
            go_d_dic[items[idx_dic['GO_ID']]].setdefault('n',items[idx_dic['Name']])
            go_d_dic[items[idx_dic['GO_ID']]].setdefault('d',items[idx_dic['Description']])
        return go_d_dic


class RINE_GOENRICH_PARSER:
    def mk_enrich_dic(self, go_enrich_fn):
        enrich_dic = dict()
        for lines in self.fileiter(go_enrich_fn):
            if not lines:
                continue
            design = self.ext_design(lines)
            go_p_dic = self.ext_go_pvalue(lines)
            enrich_dic.setdefault(design, go_p_dic)
        return enrich_dic

    def ext_design(self, lines):
        for line in lines:
            items = line.rstrip('\n').split('\t')
            if items[0] in ['NAME']:
                return items[1]
        return False

    def ext_go_pvalue(self, lines):
        go_p_dic = dict()
        for line in lines:
            items = line.rstrip('\n').split('\t')
            if items[0] in ['GO']:
                go_id = items[1]
                pvalue = items[2]
                go_p_dic.setdefault(go_id, pvalue)
        return go_p_dic

    def fileiter(self, fn):
        lines = list()
        for line in open(fn):
            if line.rstrip('\n') in ['//']:
                yield lines
                lines = []
                continue
            lines.append(line)
        yield lines


class GOENRICH_MAKER(RINE_GENES_PARSER, RINE_SAMPLE_IN_PARSER, RINE_GOENRICH_PARSER, RINE_GO2GENES_PARSER):
    def __init__(self, args):
        self.sam_in_fn = args.sample_in
        self.genes_fn = args.genes
        self.go2gene_fn = args.go_to_genes
        self.go_enrich_fn = args.go_enrich
        self.outprefix = args.outprefix
        #
        self.sam_label_dic = self.mk_sam_label_dic(self.sam_in_fn)
        self.label_sam_dic = self.mk_label_sam_dic(self.sam_in_fn)
        #
        self.deg_dic = self.mk_deg_dic(self.genes_fn)
        #for design, gene_dic in self.deg_dic.iteritems():
        #    print design, len(gene_dic)
        self.dego_dic = self.mk_dego_dic()
        #for design, go_dic in self.dego_dic.iteritems():
        #    print design, len(go_dic)
        #
        self.go_d_dic = self.mk_go_d_dic(self.go2gene_fn)
        #
        self.enrich_dic = self.mk_enrich_dic(self.go_enrich_fn)
        #for _design, go_p_dic in self.enrich_dic.iteritems():
        #    design = self.cvrt_design_id(_design)
        #    print _design, design, len(go_p_dic)
        #
        #
        # writing!!

    def writing(self):
        for _design, go_p_dic in self.enrich_dic.iteritems():
            go_dic = self.dego_dic[self.cvrt_design_id(_design)]
            #
            out = open('{0}_{1}.xls'.format(self.outprefix, _design),'w')
            header = ['GO_ID', 'Category', 'Name', 'Description',
                      'GO(P-val):{0}'.format(self.cvrt_design_id(_design)),
                      'GO(P-val<0.0010):{0}'.format(self.cvrt_design_id(_design))]
            out.write('{0}\n'.format('\t'.join(header)))
            for go_id in go_dic:
                if go_id in ['-']:
                    continue
                category = self.go_d_dic[go_id]['c']
                name = self.go_d_dic[go_id]['n']
                desc = self.go_d_dic[go_id]['d']
                pvalue = go_p_dic[go_id]
                is_p = self.is_p_yn(pvalue)
                #
                out.write('{0}\n'.format('\t'.join([go_id,category,name,desc,pvalue,is_p])))
            out.close()

    def is_p_yn(self, pvalue):
        if float(pvalue) < 0.001:
            return 'Y'
        else:
            return 'N'

    def cvrt_design_id(self, _id):
        _pair_s = _id.split('-')
        _ctrl_s = _pair_s[0].split('_')
        _case_s = _pair_s[1].split('_')
        #
        ctrl_s = [self.sam_label_dic[x] for x in _ctrl_s]
        case_s = [self.sam_label_dic[x] for x in _case_s]
        #
        return ' vs '.join([','.join(ctrl_s), ','.join(case_s)])

    def mk_dego_dic(self):
        dego_dic = dict()
        for design, gene_dic in self.deg_dic.iteritems():
            dego_dic.setdefault(design, {})
            for gene, go in gene_dic.iteritems():
                gos = go.split(',')
                for go in gos:
                    dego_dic[design].setdefault(go, [])
                    if not gene in dego_dic[design][go]:
                        dego_dic[design][go].append(gene)
                #
            #
        return dego_dic



def main(args):

    obj_samin = RINE_SAMPLE_IN_PARSER()
    obj_genes = RINE_GENES_PARSER()
    obj_gotogenes = RINE_GO2GENES_PARSER()
    obj_goenrich = RINE_GOENRICH_PARSER()

    obj = GOENRICH_MAKER(args)
    obj.writing()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='To make go_enrich.xls')
    parser.add_argument('--sample-in', help="./sample.denovo.in",
            default='./sample.denovo.in')
    parser.add_argument('--genes', help="./Rine_Denovo/report/Files/genes.xls",
            default='./Rine_Denovo/report/Files/genes.xls')
    parser.add_argument('--go-to-genes', help="./Rine_Denovo/report/Files/go_to_genes_with_ancestor.xls",
            default='./Rine_Denovo/report/Files/go_to_genes_with_ancestor.xls')
    parser.add_argument('--go-enrich', help="./Rine_Denovo/output/GO/go_enrich.txt",
            default='./Rine_Denovo/output/GO/go_enrich.txt')
    parser.add_argument('--outprefix', help="outprefix of output files",
            default='go_enrich')
    args = parser.parse_args()
    main(args)


