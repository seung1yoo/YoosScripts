import os
import math

def goDescDicMaker(go_table):
    goDescDic = dict()
    for line in open(go_table):
        items = line.strip().split('\t')
        goDescDic.setdefault(items[0], items[2])
    #
    return goDescDic

def annoDicMaker(xls, title_key, go_key, probe_key, gene_key, desc_key, go_delimiter):
    annoDic = dict()
    for line in open(xls):
        items = line.strip().split('\t')
        if items[0] in [title_key]:
            go_idx = items.index(go_key)
            probe_idx = items.index(probe_key)
            gene_idx = items.index(gene_key)
            desc_idx = items.index(desc_key)
            continue
        if items[go_idx] in ['-']:
            continue
        for go in items[go_idx].split(go_delimiter):
            go = go.strip("\"")
            annoDic.setdefault(go, {}).setdefault(\
		items[probe_idx], \
		(items[gene_idx], items[desc_idx]))
    #
    print ('## No. Gene Ontology in {0} : {1}'.format(xls, len(annoDic)))
    return annoDic

def gmtMaker(outprefix, goDescDic, annoDic):
    out = open('{0}.gmt'.format(outprefix), 'w')
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
        go_desc = go_desc.replace("/","").replace("'",\
			"").replace("#","").replace("?","")
        out.write('{0}\t{1}\t{2}\n'.format(\
		go_desc, go, '\t'.join(gs_genesDic.keys())))
    out.close()
    #
    return '{0}.gmt'.format(outprefix)

def chipMaker(outprefix, annoDic):
    chipDic = dict()
    for go, probeDic in annoDic.items():
        for probe, infos in probeDic.items():
            chipDic.setdefault(probe, infos)
    #
    out = open('{0}.chip'.format(outprefix), 'w')
    out.write('Probe Set ID\tGene Symbol\tGene Title\n')
    for probe, infos in chipDic.items():
        out.write('{0}\t{1}\n'.format(probe, '\t'.join(infos)))
    out.close()
    #
    return '{0}.chip'.format(outprefix)

def lister(idxs, items):
    some_list = []
    for idx in idxs:
        some_list.append(items[idx])
    try: some_list = [float(x) for x in some_list]
    except ValueError: some_list = [float(0) for x in some_list] 
    return some_list

def expFilter(xls, log2fc, statistics, value, title_key, probe_key, gene_key, desc_key):
    expDic = dict()
    expOutDic = dict()
    for line in open(xls):
        filtered = 0
        items = line.strip().split('\t')
        if items[0] in [title_key]:
            import re
            probe_idx = items.index(probe_key)
            gene_idx = items.index(gene_key)
            desc_idx = items.index(desc_key)
            exp_idxs = [i for i, item in enumerate(items) if re.search('^EXP', item)]
            samples = [items[x].split(':')[1] for x in exp_idxs]
            #log2fc_idxs = [i for i, item in enumerate(items) if re.search('Log2FC$', item.strip("\""))]
            log2fc_idxs = [i for i, item in enumerate(items) if re.search('log2fc$', item.strip("\"").lower())]
            log2fc_cols = [items[x] for x in log2fc_idxs]
            if statistics in ['p']:
                v_idxs = [i for i, item in enumerate(items) if re.search('PV$', item.strip("\""))]
                v_cols = [items[x] for x in v_idxs]
            elif statistics in ['q']:
                v_idxs = [i for i, item in enumerate(items) if re.search('QV$', item.strip("\""))]
                v_cols = [items[x] for x in v_idxs]
            continue
        #
        exp_list = lister(exp_idxs, items)
        if max(exp_list) < 0.3:
            filtered = 1
        #
        log2fc_list = lister(log2fc_idxs, items)
        if float('inf') in log2fc_list or float('-inf') in log2fc_list:
            exp_list_temps = [x+0.1 for x in exp_list]
            log2fc_list = []
            for v1 in exp_list_temps:
                for v2 in exp_list_temps:
                    log2fc_temp = math.log(v1,2)-math.log(v2,2)
                    log2fc_list.append(log2fc_temp)
            log2fc_list = [abs(x) for x in log2fc_list]
            if max(log2fc_list) < float(log2fc):
                filtered = 2
        else:
            log2fc_list = [abs(x) for x in log2fc_list]
            if max(log2fc_list) < float(log2fc):
                filtered = 2
        #
        v_list = lister(v_idxs, items)
        if min(v_list) > float(value):
            filtered = 3
        #
        if filtered in [0]:
            expDic.setdefault(items[gene_idx], {}).setdefault('desc', items[desc_idx])
            expDic.setdefault(items[gene_idx], {}).setdefault('exp', exp_list)
            expDic.setdefault(items[gene_idx], {}).setdefault('log2fc', log2fc_list)
            expDic.setdefault(items[gene_idx], {}).setdefault('v', v_list)
        else:
            expOutDic.setdefault(items[gene_idx], {}).setdefault('desc', items[desc_idx])
            expOutDic.setdefault(items[gene_idx], {}).setdefault('exp', exp_list)
            expOutDic.setdefault(items[gene_idx], {}).setdefault('log2fc', log2fc_list)
            expOutDic.setdefault(items[gene_idx], {}).setdefault('v', v_list)
        #
    return expDic, expOutDic, samples, log2fc_cols, v_cols

def expMaker(outprefix, expDic, samples, orderedSamples):
    out = open('{0}.exp.txt'.format(outprefix), 'w')
    out.write('NAME\tDESCRIPTION\t{0}\n'.format('\t'.join(orderedSamples)))
    #
    ordered_idxs = []
    for asample in orderedSamples:
        ordered_idxs.append(samples.index(asample))
    #
    for gene, infoDic in expDic.items():
        exps = [str(x) for x in infoDic['exp']]
        ordered_exps = []
        for idx in ordered_idxs:
            ordered_exps.append(exps[idx])
        out.write('{0}\t{1}\t{2}\n'.format(gene, infoDic['desc'], 
            '\t'.join(ordered_exps)))
    out.close()
    #
    return '{0}.exp.txt'.format(outprefix)

def gseaSampleParser(gsea_sample):
    gseaSampleDic = dict()
    for sample_group in gsea_sample.split(':'):
        sample, group = sample_group.split(',')
        gseaSampleDic.setdefault(group, []).append(sample)
    #
    orderedSamples = []
    for group, samples in sorted(gseaSampleDic.items()):
        for sample in sorted(samples):
            orderedSamples.append(sample)
    #
    return gseaSampleDic, orderedSamples

def clsMaker(outprefix, gseaSampleDic, orderedSamples):
    out = open('{0}.cls'.format(outprefix), 'w')
    out.write('{0} {1} 1\n'.format(len(orderedSamples), len(gseaSampleDic)))
    out.write('# {0}\n'.format(' '.join(sorted([x.replace('-','') for x in gseaSampleDic.keys()]))))
    #
    orderedGroups = []
    for asample in orderedSamples:
        for group, samples in sorted(gseaSampleDic.items()):
            if asample in samples:
                orderedGroups.append(group)
    #
    out.write('{0}\n'.format(' '.join([x.replace('-','') for x in orderedGroups])))
    out.close()
    #
    return '{0}.cls'.format(outprefix)

def main(args):
    title_key = args.primarykeys[0]
    go_key = args.primarykeys[1]
    probe_key = args.primarykeys[2]
    gene_key = args.primarykeys[3]
    desc_key = args.primarykeys[4]
    #### take GO infos.
    goDescDic = goDescDicMaker(args.go_table)
    annoDic = annoDicMaker(args.xls, title_key, go_key, probe_key, gene_key, desc_key, args.go_delimiter)
    #### make the output prefix
    if not os.path.isdir('{0}'.format(args.outdir)):
        os.system('mkdir -p {0}'.format(args.outdir))
    outprefix = '{0}/{1}'.format(args.outdir, args.species)
    #### make gmt & chip
    gmt_file = gmtMaker(outprefix, goDescDic, annoDic)
    chip_file = chipMaker(outprefix, annoDic)
    #### make cls
    gseaSampleDic, orderedSamples = gseaSampleParser(args.gsea_sample)
    cls_file = clsMaker(outprefix, gseaSampleDic, orderedSamples)
    #### make exp
    expDic, expOutDic, samples, log2fc_cols, v_cols = expFilter(\
            args.xls, args.log2fc, args.statistics, args.value, \
            title_key, probe_key, gene_key, desc_key)
    exp_file = expMaker(outprefix, expDic, samples, orderedSamples)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--go-table', help='My go database',
            default='~/YoosScripts/Databases/Ontology_Archives/gene_ontology_edit.obo.2016-06-01.table')
    parser.add_argument('-x', '--xls', help='genes.xls',
            default='gene.xls')
    parser.add_argument('-o', '--outdir', help='outdir name',
            default='GSEA')
    parser.add_argument('-s', '--species', help='species name or prefix name',
            default='test_species')
    parser.add_argument('-f', '--log2fc', help='exp table filter',
            default=0)
    parser.add_argument('-m', '--statistics', choices=('p','q'),
            default='p')
    parser.add_argument('-v', '--value', help='filter value p-value or q-value',
            default=1)
    parser.add_argument('-gd', '--go-delimiter',
            default=',')
    parser.add_argument('-gs', '--gsea-sample',
            default='S0001,G01:S0002,G01:S0003,G02:S0004,G02')
    parser.add_argument('-pks', '--primarykeys', nargs='+', help='Title_key GO_key Probe_key Gene_key Desc_key',
            default=['Order', 'GeneOntology', 'GeneAcc', 'GeneId', 'Desc'])
    args = parser.parse_args()
    main(args)
