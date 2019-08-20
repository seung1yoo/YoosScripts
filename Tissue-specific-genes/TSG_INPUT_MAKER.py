





class TSG_INPUT_MAKER:
    def __init__(self):
        self.load_mappedreadtable()
        self.load_genesxls()
        self.load_outprefix()
        self.load_samples()
        #
        self.gene_dic = dict()
        self.add_len()
        self.add_fpkm()
        self.add_mappedread()
        #print(self.gene_dic['TBIG000001'])
        #
        self.make_input()
        self.make_anno()

    def load_mappedreadtable(self): self.mappedreadtable = args.mappedreadtable
    def load_genesxls(self): self.genesxls = args.genesxls
    def load_outprefix(self): self.outprefix = args.outprefix
    def load_samples(self): self.samples = args.samples.split(':')
    def add_len(self):
        for line in open(self.genesxls):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['Order']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            g_id = items[idx_dic['GeneId']]
            start = items[idx_dic['Start']]
            end   = items[idx_dic['End']]
            pos_s = sorted([int(start), int(end)])
            g_len = pos_s[-1] - pos_s[0] + 1
            self.gene_dic.setdefault(g_id, {}).setdefault('len', str(g_len))
        #
    def add_fpkm(self):
        for line in open(self.genesxls):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['Order']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            g_id = items[idx_dic['GeneId']]
            for sample in self.samples:
                self.gene_dic.setdefault(g_id, {}).setdefault('{0}_fpkm'.format(sample), items[idx_dic['EXP:{0}:FPKM'.format(sample)]])
            #
        #
    def add_mappedread(self):
        for line in open(self.mappedreadtable):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['ID']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            g_id = items[idx_dic['ID']]
            for sample in self.samples:
                self.gene_dic.setdefault(g_id, {}).setdefault('{0}_count'.format(sample), items[idx_dic['{0}'.format(sample)]])
            #
        #
    def make_input(self):
        out = open('{0}.input'.format(self.outprefix),'w')
        headers = ['TrID','len']
        for sample in self.samples:
            headers.append('{0}_count'.format(sample))
        for sample in self.samples:
            headers.append('{0}_fpkm'.format(sample))
        out.write('{0}\n'.format('\t'.join(headers)))
        for g_id, info_dic in sorted(self.gene_dic.items()):
            items = [g_id]
            items.append(info_dic['len'])
            for sample in self.samples:
                items.append(info_dic['{0}_count'.format(sample)])
            for sample in self.samples:
                items.append(info_dic['{0}_fpkm'.format(sample)])
            out.write('{0}\n'.format('\t'.join(items)))
        out.close()

    def make_anno(self):
        out = open('{0}.anno'.format(self.outprefix),'w')
        for line in open(self.genesxls):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['Order']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            out.write('{0}\n'.format('\t'.join(items[1:])))
        out.close()



def main(args):
    tim = TSG_INPUT_MAKER()




if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('mappedreadtable')
    parser.add_argument('genesxls')
    parser.add_argument('outprefix')
    parser.add_argument('samples', help='sep=":"')
    args = parser.parse_args()
    main(args)
