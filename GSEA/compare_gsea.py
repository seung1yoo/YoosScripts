#!/usr/bin/python

import math

class COMPARE_GSEA:
    def __init__(self, conf, p_cut, prefix):
        self.conf = conf
        self.p_cut = float(p_cut)
        self.prefix = prefix
        #
        self.gseaDic = self.merge_report()
        result_fn = self.compare_gsea_writer(self.gseaDic, 'origin')
        #
        self.filterP(result_fn, 'P{0}'.format(str(self.p_cut)).replace('.',''))

    def filterP(self, fn, subprefix):
        out = open('{0}.{1}.xls'.format(self.prefix, subprefix), 'w')
        out_ml = open('{0}.{1}.minuslog.xls'.format(self.prefix, subprefix), 'w')

        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['GENESET_TERM']:
                out.write(line)
                out_ml.write(line)
                continue

            ps = [float(x) for x in items[1:]]

            if min(ps) < self.p_cut:
                out.write(line)
                #
                new_items = [items[0]]
                new_items.extend([self.p_to_minuslog(p) for p in ps])
                out_ml.write('{0}\n'.format('\t'.join([str(item) for item in new_items])))
            else:
                pass

        out.close()
        out_ml.close()

    def p_to_minuslog(self, p, _max=10.0):
        if p in [0.0]:
            p_minuslog = _max
        elif p in [1.0]:
            p_minuslog = 0.00001
        else:
            p_minuslog = round(-math.log(p, 10), 5)
        if p_minuslog > _max:
            p_minuslog = _max
        return p_minuslog

    def compare_gsea_writer(self, aDic, subprefix):
        out = open('{0}.{1}.xls'.format(self.prefix, subprefix), 'w')
        titles = ['GENESET_TERM']
        for gs_term, gidDic in aDic.iteritems():
            for gsea_id, epDic in sorted(gidDic.iteritems()):
                for enrich_pheno, p in sorted(epDic.iteritems()):
                    titles.append('{0}|{1}'.format(gsea_id, enrich_pheno))
            break # Do not delete
        out.write('{0}\n'.format('\t'.join(titles)))
        #
        for gs_term, gidDic in aDic.iteritems():
            new_items = [gs_term]
            for title in titles:
                if '|' not in title:
                    continue
                gsea_id, enrich_pheno = title.split('|')
                new_items.append(gidDic[gsea_id][enrich_pheno])
            out.write('{0}\n'.format('\t'.join([str(x) for x in new_items])))
        out.close()
        print 'FILE GENERATED : {0}.{1}.xls'.format(self.prefix, subprefix)
        return '{0}.{1}.xls'.format(self.prefix, subprefix)

    def merge_report(self):
        gseaDic = dict()
        gseaSubDic = dict()
        for line in open(self.conf):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split()
            gsea_id = items[0]
            enrich_pheno = items[1]
            gseaSubDic.setdefault(gsea_id, \
                    {}).setdefault(enrich_pheno, None)
            #
            fn = items[2]
            gsDic = self.gsea_report_parser(fn)
            for gs_term, p in gsDic.iteritems():
                gseaDic.setdefault(gs_term, {})
            #
        #
        for gs_term, subDic in gseaDic.iteritems():
            for gsea_id, epDic in gseaSubDic.iteritems():
                for enrich_pheno, tmp in epDic.iteritems():
                    gseaDic[gs_term].setdefault(gsea_id,\
                            {}).setdefault(enrich_pheno, 1)
                #
            #
        #
        for line in open(self.conf):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split()
            gsea_id = items[0]
            enrich_pheno = items[1]
            fn = items[2]
            gsDic = self.gsea_report_parser(fn)
            for gs_term, p in gsDic.iteritems():
                if gsea_id in gseaDic[gs_term]:
                    gseaDic[gs_term][gsea_id][enrich_pheno] = p
                #
            #
        return gseaDic

    def gsea_report_parser(self, fn):
        gsDic = dict()
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['NAME']:
                #NAME
                #GS<br> follow link to MSigDB
                #GS DETAILS
                #SIZE
                #ES
                #NES
                #NOM p-val
                #FDR q-val
                #FWER p-val
                #RANK AT MAX
                #LEADING EDGE
                idxDic = dict()
                for idx, item in enumerate(items):
                    idxDic.setdefault(item, idx)
                continue
            gs_term = items[idxDic['NAME']]
            p = float(items[idxDic['NOM p-val']])
            #p = float(items[idxDic['NES']])
            gsDic.setdefault(gs_term, p)


        return gsDic

def main(args):
    print args
    compair_gsea = COMPARE_GSEA(args.conf, args.p_cut, args.prefix)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('conf')
    parser.add_argument('p_cut')
    parser.add_argument('prefix')
    args = parser.parse_args()
    main(args)


