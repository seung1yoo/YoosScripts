
import sys
from Bio import SeqIO

class GOSEQ:
    def __init__(self, cds_fa, geneanno_fn, geneanno_iprscan_fn, selected_list_fn, outprefix):
        self.cds_fa = cds_fa
        self.geneanno_fn = geneanno_fn
        self.geneanno_iprscan_fn = geneanno_iprscan_fn
        self.selected_list_fn = selected_list_fn
        self.outprefix = outprefix

        self.g2t_dic = dict()
        self.t2g_dic = dict()
        self.anno_dic = dict()
        self.get_map()
        self.get_length()
        self.get_go()
        self.get_deg()

        #print(self.anno_dic['GENE00126.1'])

        self.write_deg()
        self.write_mapping()
        self.write_len()

    def write_len(self):
        out_fn = '{0}.len.tsv'.format(self.outprefix)
        out_fh = open(out_fn, 'w')
        for acc, info_dic in sorted(self.anno_dic.items()):
            items = [acc]
            items.append(info_dic['len'])
            out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()


    def write_deg(self):
        out_fn = '{0}.deg.tsv'.format(self.outprefix)
        out_fh = open(out_fn, 'w')
        for acc, info_dic in sorted(self.anno_dic.items()):
            items = [acc]
            items.append(info_dic['deg'])
            out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()

    def write_mapping(self):
        out_fn = '{0}.mapping.tsv'.format(self.outprefix)
        out_fh = open(out_fn, 'w')
        for acc, info_dic in sorted(self.anno_dic.items()):
            items = [acc]
            if len(info_dic['go']) in [0]:
                items.append('NA')
            else:
                items.append(','.join(info_dic['go']))
            out_fh.write('{0}\n'.format('\t'.join(items)))
        out_fh.close()

    def get_deg(self):
        deg_s = list()
        for line in open(self.selected_list_fn):
            deg_s.append(line.rstrip('\n'))
        for acc, info_dic in self.anno_dic.items():
            if acc in deg_s:
                self.anno_dic[acc]['deg'] = '1'
            else:
                self.anno_dic[acc]['deg'] = '0'

    def get_go(self):
        for line in open(self.geneanno_iprscan_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#Query']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            t_acc = items[idx_dic['#Query']]
            if len(items) in [10,11]:
                continue
            elif len(items) in [12,13]:
                t_go_s = items[idx_dic['GO_Annotation']].split('|')
                for go in t_go_s:
                    if not go.strip():
                        continue
                    if not go in self.anno_dic[t_acc]['go']:
                        self.anno_dic[t_acc]['go'].append(go)

    def _get_go(self):
        go_dic = dict()
        for line in open(self.geneanno_iprscan_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#Query']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            t_acc = items[idx_dic['#Query']]
            g_acc = self.t2g_dic[t_acc]
            if len(items) in [10,11]:
                continue
            elif len(items) in [12,13]:
                t_go_s = items[idx_dic['GO_Annotation']].split('|')
                go_dic.setdefault(g_acc, {}).setdefault(t_acc, t_go_s)
        for g_acc, t_dic in go_dic.items():
            for t_acc, t_go_s in t_dic.items():
                for go in t_go_s:
                    if not go.strip():
                        continue
                    self.anno_dic[g_acc].setdefault('go', [])
                    if not go in self.anno_dic[g_acc]['go']:
                        self.anno_dic[g_acc]['go'].append(go)

    def get_length(self):
        for record in SeqIO.parse(open(self.cds_fa), 'fasta'):
            t_acc = record.id
            t_len = len(record.seq)
            self.anno_dic[t_acc]['len'] = str(t_len)

    def _get_length(self):
        len_dic = dict()
        for line in open(self.geneanno_iprscan_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#Query']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            t_acc = items[idx_dic['#Query']]
            g_acc = self.t2g_dic[t_acc]
            t_len = items[idx_dic['SeqLength']]
            len_dic.setdefault(g_acc, {}).setdefault(t_acc, t_len)
        for g_acc, t_dic in len_dic.items():
            len_s = list()
            for t_acc, t_len in t_dic.items():
                len_s.append(int(t_len))
            if g_acc in self.anno_dic:
                self.anno_dic[g_acc].setdefault('len',max(len_s))
            else:
                print('ERROR : {0}'.format(g_acc))
                sys.exit()

    def get_map(self):
        for line in open(self.geneanno_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['GeneAcc']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            g_acc = items[idx_dic['GeneAcc']]
            t_acc = items[idx_dic['TranAcc']]
            self.g2t_dic.setdefault(g_acc, t_acc)
            self.t2g_dic.setdefault(t_acc, g_acc)
            #self.anno_dic.setdefault(g_acc, {})
            self.anno_dic.setdefault(t_acc, {}).setdefault('len', '0')
            self.anno_dic.setdefault(t_acc, {}).setdefault('go', [])
            self.anno_dic.setdefault(t_acc, {}).setdefault('deg', '0')




def main(args):
    print(args)

    if sys.argv[1] in ['goseq']:
        goseq = GOSEQ(args.cds_fa, args.geneannotation, args.geneannotation_iprscan, args.selected_list, args.outprefix)




if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser = subparsers.add_parser('goseq')
    subparser.add_argument('--cds-fa',
        default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD200405/20200602_gsea/geneset.cds.fa')
    subparser.add_argument('--geneannotation',
        default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD200405/20200602_gsea/GeneAnnotation.xls')
    subparser.add_argument('--geneannotation-iprscan',
        default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD200405/20200602_gsea/GeneAnnotation.iprscan.xls')
    subparser.add_argument('--selected-list',
        default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD200405/20200602_gsea/selected_list.txt')
    subparser.add_argument('--outprefix',
        default='/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD200405/20200602_gsea/goseq_input')

    args = parser.parse_args()
    main(args)
