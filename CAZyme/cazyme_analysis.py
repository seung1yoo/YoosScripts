import re

##############################################################################################

def cazyme_filter(inputGenes):
    outfile = 'genes.only.cazyme.xls'
    out = open(outfile, 'w')
    for line in open(inputGenes):
        items = line.rstrip('\n').split('\t')
        if items[0].startswith('Order'):
            out.write(line)
            for idx, item in enumerate(items):
                if item in ['CAZyme']:
                    cazyme_idx = idx
            continue
        #
        cazyme = items[cazyme_idx]
        if cazyme in ['-']:
            continue
        #
        out.write(line)
    out.close()
    return outfile

##############################################################################################

def classes_profiling_writer(aDic, outfile):
    out = open(outfile, 'w')
    out.write('CAZyme_class\tCAZyme_family\tCAZyme_subfamily\tCount.genes\tMember.genes\n')
    for cazyme_class, familyDic in aDic.iteritems():
        for family, geneDic in familyDic.iteritems():
            subfamilys = []
            genes = []
            for gene_id, cazyme_ids in geneDic.iteritems():
                genes.append('{0}({1})'.format(gene_id, ';'.join(cazyme_ids)))
                for cazyme_id in cazyme_ids:
                    if cazyme_id not in subfamilys:
                        subfamilys.append(cazyme_id)
                #
            out.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(cazyme_class, family, ','.join(subfamilys),\
                str(len(geneDic)), ','.join(genes)))
    out.close()

def cazyme_parser(cazyme_ids):
    cazymeDic = dict()
    for cazyme_id in cazyme_ids:
        cazyme_class = re.sub('[0-9\_]', '', cazyme_id)
        cazyme_family = cazyme_id.split('_')[0]
        cazymeDic.setdefault(cazyme_class, {}).setdefault(cazyme_family, []).append(cazyme_id)
    return cazymeDic

def classes_profiling_aDicUpdater(aDic, items, gene_idx, cazyme_idx):
    gene_id = items[gene_idx]
    cazyme_ids = items[cazyme_idx].split(',')
    cazymeDic = cazyme_parser(cazyme_ids)
    for cazyme_class, cazyme_family_dic in cazymeDic.iteritems():
        for cazyme_family, cazyme_ids in cazyme_family_dic.iteritems():
            aDic.setdefault(cazyme_class, {}).setdefault(\
                    cazyme_family, {}).setdefault(gene_id, cazyme_ids)
        #
    return aDic

def classes_profiling(cazyme_file, outfile):
    aDic = dict()
    for line in open(cazyme_file):
        items = line.rstrip('\n').split('\t')
        if items[0].startswith('Order'):
            for idx, item in enumerate(items):
                if item in ['CAZyme']:
                    cazyme_idx = idx
                elif item in ['GeneAcc']:
                    gene_idx = idx
            continue
        #
        aDic = classes_profiling_aDicUpdater(aDic, items, gene_idx, cazyme_idx)
    #
    classes_profiling_writer(aDic, outfile)

def classes_profiling_deg(cazyme_file, outfile, select_idxs):
    aDic = dict()
    for line in open(cazyme_file):
        items = line.rstrip('\n').split('\t')
        if items[0].startswith('Order'):
            for idx, item in enumerate(items):
                if item in ['CAZyme']:
                    cazyme_idx = idx
                elif item in ['GeneAcc']:
                    gene_idx = idx
            continue
        #
        selections = []
        for select_idx in select_idxs:
            selections.append(items[select_idx])
        if 'Y' in selections:
            aDic = classes_profiling_aDicUpdater(aDic, items, gene_idx, cazyme_idx)
        else:
            continue
        #
    #
    classes_profiling_writer(aDic, outfile)

##############################################################################################

def exp_parser(items, gene_id, exp_idx_dic):
    expDic = dict()
    for exp_idx, _sample in exp_idx_dic.iteritems():
        sample = _sample.split(':')[1]
        exp = items[exp_idx]
        expDic.setdefault(gene_id, {}).setdefault(sample, exp)
    return expDic

def cazyme_selector(cazyme_file, select_idxs):
    cazymes = []
    cazymeDic = dict()
    for line in open(cazyme_file):
        items = line.rstrip('\n').split('\t')
        if items[0].startswith('Order'):
            for idx, item in enumerate(items):
                if item in ['GeneAcc']:
                    gene_idx = idx
                if item in ['CAZyme']:
                    cazyme_idx = idx
            continue
        #
        selections = []
        for select_idx in select_idxs:
            selections.append(items[select_idx])
        if 'Y' in selections:
            gene_id = items[gene_idx]
            cazyme_ids = items[cazyme_idx].split(',')
            _cazymeDic = cazyme_parser(cazyme_ids)
            for cazyme_class, _familyDic in _cazymeDic.iteritems():
                for family, cazyme_ids in _familyDic.iteritems():
                    cazymeDic.setdefault(cazyme_class,
                        {}).setdefault(family,
                        {}).setdefault(gene_id, 0)
        else:
            continue
    class_countDic = dict()
    family_countDic = dict()
    for cazyme_class, familyDic in cazymeDic.iteritems():
        class_countDic.setdefault(cazyme_class, 0)
        class_countDic[cazyme_class] += len(familyDic)
        for family, gene_idDic in familyDic.iteritems():
            family_countDic.setdefault(cazyme_class, {}).setdefault(family, 0)
            family_countDic[cazyme_class][family] += len(gene_idDic)
        #
    for cazyme_class, c_count in sorted(class_countDic.iteritems(), key=lambda (k,v):(v,k), reverse=True):
        cazymes.append(cazyme_class) ## class level
        for family, f_count in sorted(family_countDic[cazyme_class].iteritems(), key=lambda (k,v):(v,k), reverse=True):
            print cazyme_class, c_count, family, f_count
            #cazymes.append(family) ## family level
    return cazymes

def expression_profiling_aDicUpdater(aDic, items, gene_idx, cazyme_idx, exp_idx_dic):
    gene_id = items[gene_idx]
    #
    expDic = exp_parser(items, gene_id, exp_idx_dic)
    for sample, exp in expDic[gene_id].iteritems():
        aDic.setdefault(gene_id, {}).setdefault(
                'exp', {}).setdefault(sample, 0)
        aDic[gene_id]['exp'][sample] += float(exp)
    #
    cazyme_ids = items[cazyme_idx].split(',')
    cazymeDic = cazyme_parser(cazyme_ids)
    for cazyme_class, familyDic in cazymeDic.iteritems():
        for family, cazyme_ids in familyDic.iteritems():
            for cazyme_id in cazyme_ids:
                aDic.setdefault(gene_id, {}).setdefault(\
                    'cazyme', {}).setdefault(\
                    cazyme_class, {}).setdefault(\
                    family, []).append(cazyme_id)
    return aDic

def expression_profiling_deg(cazyme_file, outfile, select_idxs):
    aDic = dict()
    cazymes = cazyme_selector(cazyme_file, select_idxs)
    print cazymes
    for line in open(cazyme_file):
        items = line.rstrip('\n').split('\t')
        if items[0].startswith('Order'):
            exp_idx_dic = dict()
            samples = []
            for idx, item in enumerate(items):
                if item in ['CAZyme']:
                    cazyme_idx = idx
                elif item in ['GeneAcc']:
                    gene_idx = idx
                elif item.startswith('EXP'):
                    exp_idx_dic.setdefault(idx, item)
                    samples.append(item.split(':')[1])
            continue
        #
        gene_id = items[gene_idx]
        if gene_id in ['-']:
            continue
        #
        selections = []
        for select_idx in select_idxs:
            selections.append(items[select_idx])
        if 'Y' in selections:
            aDic = expression_profiling_aDicUpdater(aDic, items, gene_idx, cazyme_idx, exp_idx_dic)
        else:
            continue
        #
    out = open(outfile, 'w')
    out.write('{0}\t{1}\t{2}\n'.format(\
        '\t'.join(['gene_id']),\
        '\t'.join(samples), '\t'.join(cazymes)))
    for gene_id, infoDic in aDic.iteritems():
        items = [gene_id.upper()]
        for sample in samples:
            exp = infoDic['exp'][sample]
            items.append(exp)
        for cazyme in cazymes:
            if cazyme in infoDic['cazyme']:
                items.append('1')
            else:
                items.append('0')
        out.write('{0}\n'.format('\t'.join([str(x) for x in items])))
    out.close()

##############################################################################################

def main(args):
    cazyme_file = cazyme_filter(args.inputGenes)
    #
    classes_profiling(cazyme_file, 'cazyme_calsses_profiling_total.xls')
    classes_profiling_deg(cazyme_file, 'cazyme_calsses_profiling_W-Bp.xls', [21])
    classes_profiling_deg(cazyme_file, 'cazyme_calsses_profiling_W-B.xls', [25])
    classes_profiling_deg(cazyme_file, 'cazyme_calsses_profiling_Bp-B.xls', [29])
    classes_profiling_deg(cazyme_file, 'cazyme_calsses_profiling_DEG.xls', [21,25,29])
    #
    #expression_profiling(cazyme_file, 'cazyme_expression_profiling_total.xls')
    expression_profiling_deg(cazyme_file, 'cazyme_expression_profiling_DEG.xls', [21,25,29])

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ig', '--inputGenes',
            help='RINE genes.xls was annotated CAZyme',
            default='../Assembly/genes.anno.clean.CAZyme.xls')
    args = parser.parse_args()
    main(args)
