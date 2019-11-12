


from bs4 import BeautifulSoup
import glob
import os


class SumupGSEA:
    def __init__(self):
        pass
        
    def load_index_html(self, index_html):
        self.index_html = index_html

    def parse_index_html(self):
        self.gs_stats_dic = dict()
        bs = BeautifulSoup(open(self.index_html), 'html.parser')
        for unit_div in bs.find_all("div"):
            for unit_h4 in unit_div.find_all("h4"):
                if "Enrichment in phenotype" in unit_h4.contents[0]:
                    phenotype = unit_h4.b.contents[0]
                    #print(phenotype)
                    self.gs_stats_dic.setdefault(phenotype, {})
            for unit_ul in unit_div.find_all("ul"):
                for unit_li in unit_ul.find_all("li"):
                    if "gene sets are upregulated in phenotype" in unit_li.contents[0]:
                        gs_enriched = unit_li.contents[0].split()[0]
                        gs_target  = unit_li.contents[0].split()[2]
                        #print(gs_enriched, gs_target)
                        self.gs_stats_dic[phenotype].setdefault('enriched', gs_enriched)
                        self.gs_stats_dic[phenotype].setdefault('target', gs_target)
                    if "at FDR < 25%" in unit_li.contents[0]:
                        gs_fdr25 = unit_li.contents[0].split()[0]
                        #print(gs_fdr25)
                        self.gs_stats_dic[phenotype].setdefault('FDR25', gs_fdr25)
                    if "at nominal pvalue < 1%" in unit_li.contents[0]:
                        gs_p1 = unit_li.contents[0].split()[0]
                        #print(gs_p1)
                        self.gs_stats_dic[phenotype].setdefault('P001', gs_p1)
                    if "at nominal pvalue < 5%" in unit_li.contents[0]:
                        gs_p5 = unit_li.contents[0].split()[0]
                        #print(gs_p5)
                        self.gs_stats_dic[phenotype].setdefault('P005', gs_p5)

    def write_gs_stats(self):
        outfh = open('{0}_{1}.stats.xls'.format(args.prefix, args.nametag),'w')
        headers = ['Nametag']
        headers.append('Phenotype')
        headers.append('No.target gene-set')
        headers.append('No.enriched gene-set')
        headers.append('FDR<25%')
        headers.append('NOM p<1%')
        headers.append('NOM p<5%')
        outfh.write('{0}\n'.format('\t'.join(headers)))
        for phynotype, info_dic in sorted(self.gs_stats_dic.items()):
            a = info_dic
            items = [args.nametag]
            items.append(phynotype)
            if a:
                items.append(a['target'])
                items.append(a['enriched'])
                items.append(a['FDR25'])
                items.append(a['P001'])
                items.append(a['P005'])
            else:
                items.append('0')
                items.append('0')
                items.append('0')
                items.append('0')
                items.append('0')
            outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()

    def load_xls(self, cont_fn, case_fn):
        self.cont_xls = cont_fn
        self.cont_tag = cont_fn.split('/')[-1].split('_')[3]
        self.case_xls = case_fn
        self.case_tag = case_fn.split('/')[-1].split('_')[3]

    def parse_xls(self, filter_name, filter_value):
        self.gs_table_dic = dict()
        self.read_xls(self.cont_xls, self.cont_tag, filter_name, filter_value)
        self.read_xls(self.case_xls, self.case_tag, filter_name, filter_value)
    
    def read_xls(self, xls_fn, tag, filter_name, filter_value):
        n = 0
        for line in open(xls_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['NAME']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            if not items[idx_dic[filter_name]]:
                continue
            if float(items[idx_dic[filter_name]]) < float(filter_value):
                fv = float(items[idx_dic[filter_name]])
                n += 1 
                self.gs_table_dic.setdefault(tag, {})
                self.gs_table_dic[tag].setdefault(fv, {})
                self.gs_table_dic[tag][fv].setdefault(n, {})
                self.gs_table_dic[tag][fv][n].setdefault('name', items[idx_dic['NAME']])
                self.gs_table_dic[tag][fv][n].setdefault('size', items[idx_dic['SIZE']])
                self.gs_table_dic[tag][fv][n].setdefault('es', items[idx_dic['ES']])
                self.gs_table_dic[tag][fv][n].setdefault('nom_p', items[idx_dic['NOM p-val']])
                self.gs_table_dic[tag][fv][n].setdefault('fdr_q', items[idx_dic['FDR q-val']])
            else:
                continue

    def write_gs_table(self, filter_tag):
        outfh = open('{0}_{1}.{2}.table.xls'.format(args.prefix, args.nametag, filter_tag),'w')
        outfh_top5 = open('{0}_{1}.{2}.TOP5.table.xls'.format(args.prefix, args.nametag, filter_tag),'w')
        headers = ['Nametag']
        headers.append('Phenotype')
        headers.append('Rank')
        headers.append('gene-set name')
        headers.append('gene-set size')
        headers.append('Enrichment score')
        headers.append('NOM p')
        headers.append('FDR q')
        outfh.write('{0}\n'.format('\t'.join(headers)))
        outfh_top5.write('{0}\n'.format('\t'.join(headers)))
        for phynotype, fv_dic in self.gs_table_dic.items():
            rank = 0
            for fv, n_dic in sorted(fv_dic.items()):
                for n, info_dic in sorted(n_dic.items()):
                    rank += 1
                    a = info_dic
                    items = [args.nametag]
                    items.append(phynotype)
                    items.append(str(rank))
                    items.append(a['name'])
                    items.append(a['size'])
                    items.append(a['es'])
                    items.append(a['nom_p'])
                    items.append(a['fdr_q'])
                    outfh.write('{0}\n'.format('\t'.join(items)))
                    if rank <= 5:
                        outfh_top5.write('{0}\n'.format('\t'.join(items)))
        outfh.close()
        outfh_top5.close()


def main(args):
    sumup = SumupGSEA()
    if args.dirname and args.cont_tag and args.case_tag:
        args.index_html = os.path.join(args.dirname, 'index.html')
        args.cont_xls = glob.glob('{0}/gsea_report_for_{1}_*.xls'.format(
                                   args.dirname, args.cont_tag))[0]
        args.case_xls = glob.glob('{0}/gsea_report_for_{1}_*.xls'.format(
                                   args.dirname, args.case_tag))[0]
    if args.index_html:
        sumup.load_index_html(args.index_html)
        sumup.parse_index_html()
        sumup.write_gs_stats()
    if args.cont_xls and args.case_xls:
        sumup.load_xls(args.cont_xls, args.case_xls)
        sumup.parse_xls('NOM p-val', 0.05)
        sumup.write_gs_table('P005')
    else:
        print('ERROR! PLEASE CHECK THE args')




if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('nametag', default='DEG001_GO_03C_vs_03DP')
    parser.add_argument('prefix', default='SumupGSEA')

    parser.add_argument('--dirname')
    parser.add_argument('--cont-tag')
    parser.add_argument('--case-tag')

    parser.add_argument('--index-html')

    parser.add_argument('--cont-xls')
    parser.add_argument('--case-xls')

    args = parser.parse_args()
    main(args)
