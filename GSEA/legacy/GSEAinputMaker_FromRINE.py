
def goDescDicMaker(go_table):
    goDescDic = dict()
    for line in open(go_table):
        items = line.strip().split('\t')
        goDescDic.setdefault(items[0], items[2])
    return goDescDic

def annoDicMaker(xml):
    annoDic = dict()
    for line in open(xml):
        items = line.strip().split('\t')
        if items[0] in ['Order']: ## genes.xls version
        #if items[0] in ['GeneId']: ## genes.de.xls version
            #print (dir(items))
            go_idx = items.index('GeneOntology')
            probe_idx = items.index('GeneId')
            #gene_idx = items.index('GeneName') ## GeneName version
            gene_idx = items.index('GeneAcc') ## GeneAcc version
            desc_idx = items.index('Desc')
            continue
        if items[go_idx] in ['-']:
            continue
        for go in items[go_idx].split(','):
            go = go.strip("\"")
            annoDic.setdefault(go, {}).setdefault(items[probe_idx], (items[gene_idx], items[desc_idx]))
    return annoDic

def gmtMaker(species, goDescDic, annoDic):
    out = open('{0}.gmt'.format(species), 'w')
    for go, probeDic in annoDic.items():
        gs_genesDic = dict()
        for probe, infos in probeDic.items():
            gs_genesDic.setdefault(infos[0], 0)
        if go in goDescDic:
            go_desc = goDescDic[go].upper()
        else:
            print ('{0} has no description !!!!!'.format(go))
            import sys
            sys.exit()
        import re
        go_desc = go_desc.replace("/","").replace("'","").replace("#","").replace("?","")
        out.write('{0}\t{1}\t{2}\n'.format(go_desc, go, '\t'.join(gs_genesDic.keys())))
    out.close()
    return '{0}.gmt'.format(species)

def chipMaker(species, annoDic):
    chipDic = dict()
    for go, probeDic in annoDic.items():
        for probe, infos in probeDic.items():
            chipDic.setdefault(probe, infos)

    print ('probe# in chip: {0}'.format(len(chipDic)))
    out = open('{0}.chip'.format(species), 'w')
    out.write('Probe Set ID\tGene Symbol\tGene Title\n')
    for probe, infos in chipDic.items():
        out.write('{0}\t{1}\n'.format(probe, '\t'.join(infos)))
    out.close()
    return '{0}.chip'.format(species)

def lister(idxs, items):
    some_list = []
    for idx in idxs:
        some_list.append(items[idx])
    try: some_list = [float(x) for x in some_list]
    except ValueError: some_list = [float(0) for x in some_list] 
    return some_list

def expFilter(xml, log2fc, statistics, value):
    expDic = dict()
    expOutDic = dict()
    for line in open(xml):
        filtered = 0
        items = line.strip().split('\t')
        if items[0] in ['Order']: ## genes.xls version
        #if items[0] in ['GeneId']: ## genes.de.xls version
            import re
            probe_idx = items.index('GeneId')
            #gene_idx = items.index('GeneName') ## GeneName version
            gene_idx = items.index('GeneAcc') ## GeneAcc version
            desc_idx = items.index('Desc')
            exp_idxs = [i for i, item in enumerate(items) if re.search('^EXP', item)]
            samples = [items[x] for x in exp_idxs]
            log2fc_idxs = [i for i, item in enumerate(items) if re.search('Log2FC$', item.strip("\""))]
            log2fc_cols = [items[x] for x in log2fc_idxs]
            if statistics in ['p']:
                v_idxs = [i for i, item in enumerate(items) if re.search('PV$', item.strip("\""))]
                v_cols = [items[x] for x in v_idxs]
            elif statistics in ['q']:
                v_idxs = [i for i, item in enumerate(items) if re.search('QV$', item.strip("\""))]
                v_cols = [items[x] for x in v_idxs]
            print (probe_idx, gene_idx, desc_idx, exp_idxs, log2fc_idxs, v_idxs)
            continue

        exp_list = lister(exp_idxs, items)
        if max(exp_list) < 0.3:
            filtered = 1
        log2fc_list = lister(log2fc_idxs, items)
        if float('inf') in log2fc_list or float('-inf') in log2fc_list:
            exp_list_temps = [x+0.1 for x in exp_list]
            #log2fc_temps = []
            log2fc_list = []
            for v1 in exp_list_temps:
                for v2 in exp_list_temps:
                    import math
                    #log2fc_temp = math.log(v1,2)-math.log(v2,2)
                    #log2fc_temps.append(log2fc_temp)
                    log2fc_temp = math.log(v1,2)-math.log(v2,2)
                    log2fc_list.append(log2fc_temp)
            log2fc_list = [abs(x) for x in log2fc_list]
            if max(log2fc_list) < float(log2fc):
                filtered = 2
        else:
            log2fc_list = [abs(x) for x in log2fc_list]
            if max(log2fc_list) < float(log2fc):
                filtered = 2
        v_list = lister(v_idxs, items)
        if min(v_list) > float(value):
            filtered = 3

        if filtered in [0]:
            expDic.setdefault(items[gene_idx], {}).setdefault('desc', items[desc_idx])
            ##GeneAcc version.
            expDic.setdefault(items[gene_idx], {}).setdefault('exp', exp_list)
            '''
            ##GeneName version.
            ##there is duplicate names.
            ##because, GeneName column is annotated gene name.
            expDic.setdefault(items[gene_idx], {}).setdefault('exp', [0.0]*len(exp_list))
            expDic[items[gene_idx]]['exp'] = [exp+exp_list[i] \
                    for i, exp in enumerate(expDic[items[gene_idx]]['exp'])]
            '''
            expDic.setdefault(items[gene_idx], {}).setdefault('log2fc', log2fc_list)
            expDic.setdefault(items[gene_idx], {}).setdefault('v', v_list)
        else:
            expOutDic.setdefault(items[gene_idx], {}).setdefault('desc', items[desc_idx])
            ##GeneAcc version.
            expOutDic.setdefault(items[gene_idx], {}).setdefault('exp', exp_list)
            '''
            ##GeneName version.
            ##there is duplicate names.
            ##because, GeneName column is annotated gene name.
            expOutDic.setdefault(items[gene_idx], {}).setdefault('exp', [0.0]*len(exp_list))
            expOutDic[items[gene_idx]]['exp'] = [exp+exp_list[i] \
                    for i, exp in enumerate(expDic[items[gene_idx]]['exp'])]
            '''
            expOutDic.setdefault(items[gene_idx], {}).setdefault('log2fc', log2fc_list)
            expOutDic.setdefault(items[gene_idx], {}).setdefault('v', v_list)

    return expDic, expOutDic, samples, log2fc_cols, v_cols

def expMaker(species, expDic, samples):
    out = open('{0}.exp.txt'.format(species), 'w')
    out.write('NAME\tDESCRIPTION\t{0}\n'.format('\t'.join([x.split(':')[1] for x in samples])))
    for gene, infoDic in expDic.items():
        out.write('{0}\t{1}\t{2}\n'.format(gene, infoDic['desc'], 
            '\t'.join([str(x) for x in infoDic['exp']])))
    out.close()
    return '{0}.exp.txt'.format(species)

def heatmap_draw(exp_file):
    # case 1
    #cmd = 'python2.7 /Users/Yoo/Documents/YooScripts/hierarchical_clustering.py --i {0}'.format(exp_file)
    #os.system(cmd)
    outDir = '/'.join(exp_file.split('/')[:-1])
    cmd = "R CMD BATCH --no-save --no-restore '--args {0} {1}' /Users/Yoo/Documents/YooScripts/Heat-map/HeatmapMulti_forBATCH.R {1}/Heatmap.Rout".format(exp_file, outDir)
    print (cmd)
    os.system(cmd)


def split_go(gmt_file, expDic, expOutDic, goInfo_dir, samples, log2fc_cols, v_cols):
    for line in open(gmt_file):
        items = (line.strip().split('\t'))
        go_desc = items[0]
        go_dir = '{0}/{1}'.format(goInfo_dir, go_desc)
        go_dir = go_dir.replace("(", "").replace(")", "").replace(" ", "_")
        if os.path.isdir('{0}'.format(go_dir)): pass
        else: os.system('''mkdir "{0}"'''.format(go_dir))
        go_id = items[1]
        go_units = items[2:]
        #print (go_desc)
        #print (go_id)
        #print (go_units)

        out = open('{0}/genes.xls'.format(go_dir), 'w')
        out.write('{0}\t{1}\tDEG?\t{2}\t{3}\n'.format('\t'.join(['go_id', 'go_desc', 'gene']),
            '\t'.join(samples), '\t'.join(log2fc_cols), '\t'.join(v_cols)))
        out2 = open('{0}/exp.txt'.format(go_dir), 'w')
        out2.write('Genes\t{0}\n'.format('\t'.join([x.split(':')[1] for x in samples])))
        for gene in go_units:
            if gene in expDic:
                deg_tag = 'Yes'
                infoDic = expDic[gene]
                #print (gene, expDic[gene])
            elif gene in expOutDic:
                deg_tag = 'No'
                infoDic = expOutDic[gene]
                #print (gene, expOutDic[gene])
            line_units = []
            line_units.append(go_id)
            line_units.append(go_desc)
            line_units.append(gene)
            line_units.extend(infoDic['exp'])
            line_units.append(deg_tag)
            line_units.extend(infoDic['log2fc'])
            line_units.extend(infoDic['v'])
            out.write('{0}\n'.format('\t'.join([str(x) for x in line_units])))
            out2.write('{0}\t{1}\n'.format(gene, '\t'.join([str(x) for x in infoDic['exp']])))
        out.close()
        out2.close()
        if not len(go_units) in [1]:
            heatmap_draw('{0}/exp.txt'.format(go_dir))

def main(args):
    print (args)
    #### gmt & chip 
    goDescDic = goDescDicMaker(args.go_table)
    print ('go# in goDescDic: {0}'.format(len(goDescDic)))
    annoDic = annoDicMaker(args.xml)
    print ('go# in annotated {0}: {1}'.format(args.species, len(annoDic)))
    gmt_file = gmtMaker(args.species, goDescDic, annoDic)
    chip_file = chipMaker(args.species, annoDic)

    #### exp
    expDic, expOutDic, samples, log2fc_cols, v_cols = expFilter(
            args.xml, args.log2fc, args.statistics, args.value)
    print ('expDic: {0}'.format(len(expDic)))
    exp_file = expMaker(args.species, expDic, samples)

    #### exp, deg information add to gmt
    #goInfo_dir = 'goInfo'
    #if os.path.isdir(goInfo_dir): pass
    #else: os.system('mkdir {0}'.format(goInfo_dir))
    #split_go(gmt_file, expDic, expOutDic, goInfo_dir, samples, log2fc_cols, v_cols)

if __name__=='__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--go-table', help='My database',
            default='../Ontology_Archives/gene_ontology_edit.obo.2016-06-01.table')
    parser.add_argument('-x', '--xml', help='genes.xml',
            default='/Volumes/TBI_siyoo/TBI_Research/06_DKU/00_Prof_HanGD/Tapes/RNAseq/DKU_Han-Venerupis-2016-03_V1_Ref_P005/Files/genes.xls')
    parser.add_argument('-s', '--species', help='species name or prefix name',
            default='Tapes')
    parser.add_argument('-f', '--log2fc', help='exp table filter',
            default=1)
    parser.add_argument('-m', '--statistics', choices=('p','q'),
            default='p')
    parser.add_argument('-v', '--value', help='filter value p-value or q-value',
            default=0.05)
    args = parser.parse_args()
    main(args)
