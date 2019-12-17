import os
import math

def gsDicMaker(anno_table, xls):
    gsDic = dict()
    for line in open(anno_table):
        items = line.strip('\n').split('\t')
        gsDic.setdefault(items[0], items[2])
    return gsDic

def annoDicMaker(xls, title_key, gs_key, probe_key, gene_key, desc_keys, delimiter):
    annoDic = dict()
    for line in open(xls):
        items = line.strip().split('\t')
        if items[0] in [title_key]:
            gs_idx = items.index(gs_key)
            probe_idx = items.index(probe_key)
            gene_idx = items.index(gene_key)
            desc_idxs = []
            for desc_key in desc_keys:
                desc_idxs.append(items.index(desc_key))
            continue
        if items[gs_idx] in ['-', '0']:
            continue
        descs = []
        for desc_idx in desc_idxs:
            descs.append(items[desc_idx])
        desc = '@'.join(descs)
        for gs in items[gs_idx].split(delimiter):
            gs = gs.strip("\"")
            annoDic.setdefault(gs, {}).setdefault(\
		items[probe_idx], \
		(items[gene_idx], desc))
    return annoDic

def gmtMaker(outprefix, gsDic, annoDic):
    out = open('{0}.gmt'.format(outprefix), 'w')
    for gs, probeDic in annoDic.items():
        gs_genesDic = dict()
        for probe, infos in probeDic.items():
            gs_genesDic.setdefault(infos[0], 0)
        if gs in gsDic:
            gs_desc = gsDic[gs].lower()
        else:
            print ('{0} has no description !!!!!'.format(gs))
            import sys
            sys.exit()
            continue
        import re
        gs_desc = gs_desc.replace("/","").replace("'",\
			"").replace("#","").replace("?","")
        if len(gs_genesDic.keys()) < 5:
            print ('{0} has {1} members !!!!!'.format(gs, str(len(gs_genesDic.keys()))))
            continue
        out.write('{0}\t{1}\t{2}\n'.format(\
		gs_desc, gs, '\t'.join(gs_genesDic.keys())))
    out.close()
    #
    return '{0}.gmt'.format(outprefix)

def chipMaker(outprefix, annoDic):
    chipDic = dict()
    for gs, probeDic in annoDic.items():
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
    if len(some_list) in [0]: # This is for other format of RINE-report(genes.xls)
        some_list = [0.0] # this format has not log2fc$, PV$, QV$
    return some_list

def expFilter(xls, log2fc, statistics, value, title_key, probe_key, gene_key, desc_keys):
    expDic = dict()
    expOutDic = dict()
    for line in open(xls):
        filtered = 0
        items = line.strip().split('\t')
        if items[0] in [title_key]:
            import re
            probe_idx = items.index(probe_key)
            gene_idx = items.index(gene_key)
            desc_idxs = []
            for desc_key in desc_keys:
                desc_idxs.append(items.index(desc_key))
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
        descs = []
        for desc_idx in desc_idxs:
            descs.append(items[desc_idx])
        desc = '@'.join(descs)
        if filtered in [0]:
            expDic.setdefault(items[gene_idx], {}).setdefault('desc', desc)
            expDic.setdefault(items[gene_idx], {}).setdefault('exp', exp_list)
            expDic.setdefault(items[gene_idx], {}).setdefault('log2fc', log2fc_list)
            expDic.setdefault(items[gene_idx], {}).setdefault('v', v_list)
        else:
            expOutDic.setdefault(items[gene_idx], {}).setdefault('desc', desc)
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

def gseaSampleParser(sample_info):
    gseaSampleDic = dict()
    for sample_group in sample_info.split(','):
        sample, group = sample_group.split(':')
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
    gs_key = args.primarykeys[1]
    probe_key = args.primarykeys[2]
    gene_key = args.primarykeys[3]
    desc_keys = args.primarykeys[4].split('@')

    #### take Gene-Set infos.
    gsDic = gsDicMaker(args.anno_table, args.xls)
    annoDic = annoDicMaker(args.xls, title_key, gs_key, probe_key, gene_key, desc_keys, args.delimiter)

    #### make the output prefix
    if not os.path.isdir('{0}'.format(args.outdir)):
        os.system('mkdir -p {0}'.format(args.outdir))
    outprefix = '{0}/{1}'.format(args.outdir, args.prefix)

    #### make gmt & chip
    gmt_file = gmtMaker(outprefix, gsDic, annoDic)
    chip_file = chipMaker(outprefix, annoDic)

    #### make cls
    gseaSampleDic, orderedSamples = gseaSampleParser(args.sample_info)
    cls_file = clsMaker(outprefix, gseaSampleDic, orderedSamples)

    #### make exp
    expDic, expOutDic, samples, log2fc_cols, v_cols = expFilter(\
            args.xls, args.log2fc, args.statistics, args.value, \
            title_key, probe_key, gene_key, desc_keys)
    exp_file = expMaker(outprefix, expDic, samples, orderedSamples)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--anno-table', default='./DO_HumanDO.obo.totDOLst.txt')
    parser.add_argument('--xls', default='./01.gene.de.mod.DO.xls')
    parser.add_argument('--outdir', default='./gsea')
    parser.add_argument('--prefix', default='mouse')
    parser.add_argument('--log2fc', default=0)
    parser.add_argument('--statistics', choices=('p','q'), default='p')
    parser.add_argument('--value', default=1)
    parser.add_argument('--delimiter', default=',')
    parser.add_argument('--sample-info', default='Test_1:Con,Test_2:Con,Test_3:Case,Test_4:Case')
    parser.add_argument('--primarykeys', nargs='+', help='Title_key GS_key Probe_key Gene_key Desc_keys(Name@Desc)',
            default=['GeneId', 'DO clust', 'GeneId', 'GeneAcc', 'GeneName@Desc'])
    args = parser.parse_args()
    main(args)
