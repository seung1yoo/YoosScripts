
import os
from RINE_UpDownGoSeq import ParserSelectDE


def make_go_dic(go_enrich):
    go_dic = dict()
    for line in open(go_enrich):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['GO_ID']:
            idx_dic = dict()
            for idx, item in enumerate(items):
                idx_dic.setdefault(item, idx)
            continue
        go_id = items[idx_dic['GO_ID']]
        category = items[idx_dic['Category']]
        name = items[idx_dic['Name']]
        desc = items[idx_dic['Description']]

        go_dic.setdefault(go_id, {}).setdefault('category', category)
        go_dic.setdefault(go_id, {}).setdefault('name', name)
        go_dic.setdefault(go_id, {}).setdefault('desc', desc)

    return go_dic

def make_goseq_dic(goseq_results):
    go_seq_dic = dict()
    for fn_result in goseq_results:
        deg_name = '.'.join(fn_result.split('.')[:-3])
        deg_pair, updown = deg_name.split('/')[-1].split('.')
        go_seq_dic.setdefault(deg_pair, {}).setdefault(updown, {})
        for line in open(fn_result):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['category']:
                continue
            go_id = items[1]
            pvalue = items[2]
            go_seq_dic[deg_pair][updown].setdefault(go_id, pvalue)
    return go_seq_dic

def make_report(go_dic, goseq_dic, p_cut, name_dic):
    out = open('go_enrich.UpDown.xls', 'w')
    titles = ['GO_ID', 'Category', 'Name', 'Description']
    orders = []
    for deg_pair, updown_dic in sorted(goseq_dic.iteritems()):
        for updown, result_dic in sorted(updown_dic.iteritems()):
            titles.append('{0}:{1}:pvalue'.format(deg_pair, updown).replace(deg_pair, name_dic[deg_pair]['cus']))
            titles.append('{0}:{1}:pvalue<{2}'.format(deg_pair, updown, p_cut).replace(deg_pair, name_dic[deg_pair]['cus']))
            orders.append('{0}:{1}'.format(deg_pair, updown))
    out.write('{0}\n'.format('\t'.join(titles)))


    for go_id, info_dic in go_dic.iteritems():
        items = []
        items.append(go_id)
        items.append(info_dic['category'])
        items.append(info_dic['name'])
        items.append(info_dic['desc'])
        for order in orders:
            deg_pair, updown = order.split(':')
            result_dic = goseq_dic[deg_pair][updown]
            if go_id in result_dic:
                items.append(result_dic[go_id])
                if float(result_dic[go_id]) < float(p_cut):
                    items.append('Y')
                else:
                    items.append('N')
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()

    return 'go_enrich.UpDown.xls'

def uploader(name_dic, fn_report, goseq_results):
    if not os.path.exists('./Upload'):
        os.mkdir('./Upload')
    if not os.path.exists('./Upload/{0}'.format(fn_report)):
        os.symlink(os.path.abspath(fn_report), './Upload/{0}'.format(fn_report))
    for goseq_result in goseq_results:
        if not os.path.exists('./Upload/{0}'.format(goseq_result)):
            os.symlink(os.path.abspath(goseq_result), './Upload/{0}'.format(goseq_result))
    #
    out = open('./Upload/Sample_ID_Match.txt', 'w')
    for deg_pair, info_dic in name_dic.iteritems():
        items = list()
        items.append(deg_pair)
        items.append(info_dic['tbi'])
        items.append(info_dic['cus'])
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()

def main(args):
    go_dic = make_go_dic(args.go_enrich)
    goseq_dic = make_goseq_dic(args.goseq_results)
    #
    psde = ParserSelectDE(args)
    #
    fn_report = make_report(go_dic, goseq_dic, args.p_cut, psde.name_dic)
    #
    uploader(psde.name_dic, fn_report, args.goseq_results)


    

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('select_deg',
            help='Rine_Quant/analysis/select_de.gene.xls')
    parser.add_argument('go_enrich',
            help='/BiO/BioProjects/TBD180681-AJOU-Mouse-RNAref-20181115/Rine_Quant/report/GO/go_enrich.xls')
    parser.add_argument('goseq_results', nargs='+',
            help='*.goseq.wall.xls')
    parser.add_argument('p_cut', type=str)
    args = parser.parse_args()
    main(args)
