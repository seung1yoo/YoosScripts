import os
import sys

class ParserSelectDE:
    def __init__(self, args):
        self.in_fn = args.select_deg
        self.name_dic = self.make_name_dic()
        self.deg_dic = self.make_deg_dic()

    def make_deg_dic(self):
        deg_dic = dict()
        for line in open(self.in_fn):
            if line.startswith('# METHOD'):
                continue
            elif line.startswith('# NAME'):
                continue
            elif line.startswith('# TEST'):
                continue
            items = line.rstrip('\n').split('\t')
            gid_tbi = items[0]
            gid_acc = items[1]
            gid_sym = items[2]
            num_deg = items[3]
            info_degs = items[4:]
            for deg_pair, infoDic in self.name_dic.iteritems():
                _idx = int(infoDic['idx']) - 1
                info_deg = info_degs[_idx]
                info_dic = self.grep_deg_info(info_deg)
                deg_dic.setdefault(deg_pair, {})
                deg_dic[deg_pair].setdefault(gid_tbi, info_dic)
        return deg_dic

    def grep_deg_info(self, info_deg):
        info_dic = dict()
        infos = info_deg.split(';')
        info_dic.setdefault('exp_s1', infos[0])
        info_dic.setdefault('exp_s2', infos[1])
        info_dic.setdefault('log2fc', infos[2])
        info_dic.setdefault('p', infos[3])
        info_dic.setdefault('q', infos[4])
        info_dic.setdefault('yn', infos[5])
        return info_dic

    def make_name_dic(self):
        name_dic = dict()
        for line in open(self.in_fn):
            if line.startswith('# METHOD'):
                continue
            elif line.startswith('# NAME'):
                items = line.rstrip('\n').split('\t')
                deg_idx = items[1]
                deg_samId_pair = items[2]
                deg_tbiId_pair = items[3]
                deg_cusId_pair = items[4]

                name_dic.setdefault(deg_samId_pair, {}).\
                         setdefault('idx', deg_idx)
                name_dic.setdefault(deg_samId_pair, {}).\
                         setdefault('tbi', deg_tbiId_pair)
                name_dic.setdefault(deg_samId_pair, {}).\
                         setdefault('cus', deg_cusId_pair)
            else:
                continue
        return name_dic

    def make_goseq_deg(self, deg_pair, goseq_deg):
        fn_up = '{0}.goseq.deg.up.xls'.format(deg_pair)
        fn_dw = '{0}.goseq.deg.dw.xls'.format(deg_pair)

        out_up = open(fn_up, 'w')
        out_dw = open(fn_dw, 'w')
        #for gid_tbi, info_dic in sorted(self.deg_dic[deg_pair].iteritems()):
        for line in open(goseq_deg):
            items = line.rstrip('\n').split('\t')
            gid_tbi = items[0]
            if gid_tbi in self.deg_dic[deg_pair]:
                info_dic = self.deg_dic[deg_pair][gid_tbi]
                if info_dic['yn'] in ['True'] and info_dic['log2fc'].startswith('-'):
                    out_up.write('{0}\t0\n'.format(gid_tbi))
                    out_dw.write('{0}\t1\n'.format(gid_tbi))
                elif info_dic['yn'] in ['True'] and not info_dic['log2fc'].startswith('-'):
                    out_up.write('{0}\t1\n'.format(gid_tbi))
                    out_dw.write('{0}\t0\n'.format(gid_tbi))
                elif info_dic['yn'] in ['False']:
                    out_up.write('{0}\t0\n'.format(gid_tbi))
                    out_dw.write('{0}\t0\n'.format(gid_tbi))
                else:
                    print(info_dic)
                    print('ERROR')
                    sys.exit()
            else:
                out_up.write('{0}\t0\n'.format(gid_tbi))
                out_dw.write('{0}\t0\n'.format(gid_tbi))

        out_up.close()
        out_dw.close()

        return fn_up, fn_dw

class GoseqRunner:
    def __init__(self, args, deg_fn, updown):
        self.deg_pair_name = '{0}.{1}'.format(args.deg_pair, updown)

        self.deg_fn = deg_fn
        self.len_fn = args.goseq_len
        self.map_fn = args.goseq_mapping

        self.script_temp = '''\
library(goseq)

Input_genes<-read.table("deg_fn", head=F, sep="\t")
Input_go<-read.table("map_fn", head=F, sep="\t")
Input_genes_length<-read.table("len_fn", head=F, sep="\t")

genes<-as.integer(Input_genes[,2])
names(genes)<-Input_genes[,1]
go<-strsplit(as.character(Input_go[,2]), ",")
names(go)<-Input_go[,1]
length<-as.numeric(Input_genes_length[,2])

pwf=nullp(genes, bias.data=length, plot.fit=F)

GO.wall<-goseq(pwf, gene2cat=go, method="Wallenius")
write.table(GO.wall, "deg_pair_name.goseq.wall.xls", quote=F, sep="\t")

GO.nobias<-goseq(pwf, gene2cat=go, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Hypergeometric")
write.table(GO.nobias, "deg_pair_name.goseq.obias.xls", quote=F, sep="\t")

GO.samp<-goseq(pwf, gene2cat=go, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Sampling", repcnt=2000)
write.table(GO.samp, "deg_pair_name.goseq.sampling.xls", quote=F, sep="\t")
'''
    def write_script(self):
        out = open('{0}.goseq.R'.format(self.deg_pair_name), 'w')
        script = self.script_temp.replace('deg_fn', self.deg_fn).\
                                  replace('map_fn', self.map_fn).\
                                  replace('len_fn', self.len_fn).\
                                  replace('deg_pair_name', self.deg_pair_name)
        out.write(script)
        out.close()
        return '{0}.goseq.R'.format(self.deg_pair_name)

    def run(self, fn_rscript):
        # make sure R PATH
        print('/BiO/BioTools/R/R-3.1.1/bin/R CMD BATCH {0}'.format(fn_rscript))
        os.system('/BiO/BioTools/R/R-3.1.1/bin/R CMD BATCH {0}'.format(fn_rscript))

def main(args):
    psde = ParserSelectDE(args)
    fn_up, fn_dw = psde.make_goseq_deg(args.deg_pair, args.goseq_deg)

    gr = GoseqRunner(args, fn_up, 'up')
    fn_rscript = gr.write_script()
    gr.run(fn_rscript)

    gr = GoseqRunner(args, fn_dw, 'dw')
    fn_rscript = gr.write_script()
    gr.run(fn_rscript)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('select_deg',
            help='Rine_Quant/analysis/select_de.gene.xls')
    parser.add_argument('goseq_deg',
            help='Rine_Quant/output/GO/*.goseq.deg.xls')
    parser.add_argument('goseq_len',
            help='Rine_Quant/output/GO/*.goseq.len.xls')
    parser.add_argument('goseq_mapping',
            help='Rine_Quant/output/GO/*.goseq.mapping.xls')
    parser.add_argument('deg_pair',
            help='S0001_S0002-S0003')
    args = parser.parse_args()
    main(args)
