
def kaas_q1_parser(infile):
    keggAnnoDic = dict()
    keggTableDic = dict()
    for line in open(infile):
        if not line[0] in ['A','B','C','D']:
            continue
        #
        if line.startswith('A'):
            hierarchy_a = line.rstrip('\n').split('>')[1].split('<')[0]
        elif line.startswith('B') and '>' in line:
            hierarchy_b = line.rstrip('\n').split('>')[1].split('<')[0]
        elif line.startswith('C'):
            items = line[1:].rstrip('\n').split()
            keggBrite_id = items[0] # ko id
            keggBrite_desc = ' '.join(items[1:]) # ko desc
        elif line.startswith('D'):
            items = line[1:].rstrip('\n').split(';')
            items = [x.strip() for x in items]
            myGene = items[0]
            keggOrthology_id = items[1].split()[0] # K(kegg Orthology) id
            keggEnzyme_alias = ' '.join(items[1].split()[1:]) # EC alias
            keggEnzyme_desc = items[2] # EC desc
            keggEnzyme_id = items[2].split('[')[-1].split(']')[0] # EC id
            #
            keggAnnoDic.setdefault(myGene, {}).setdefault('keggBrite', {}).setdefault(keggBrite_id, keggBrite_desc)
            keggAnnoDic.setdefault(myGene, {}).setdefault('keggOrthology', {}).setdefault(keggOrthology_id, '')
            keggAnnoDic.setdefault(myGene, {}).setdefault('keggEnzyme', {}).setdefault(keggEnzyme_id, '')
            #
            keggTableDic.setdefault(hierarchy_a, {}).setdefault(hierarchy_b, {}).setdefault(keggBrite_id, keggBrite_desc)
        #
    return keggAnnoDic, keggTableDic

def keggUpdate(keggAnnoDic, genes_in, genes_out, title_key, myGene_key, delimiter):
    out = open(genes_out, 'w')
    for line in open(genes_in):
        items = line.rstrip('\n').split('\t')
        if items[0] in [title_key]:
            items.extend(['keggBrite_id', 'keggBrite_desc', 'keggOrthology', 'keggEnzyme'])
            out.write('{0}\n'.format('\t'.join(items)))
            myGene_idx = items.index(myGene_key)
            continue
        #
        myGene = items[myGene_idx]
        if not myGene in keggAnnoDic:
            items.extend(['-','-','-'])
            out.write('{0}\n'.format('\t'.join(items)))
        else:
            keggBrite_ids = []
            keggBrite_descs = []
            for keggBrite_id, keggBrite_desc in keggAnnoDic[myGene]['keggBrite'].items():
                keggBrite_ids.append(keggBrite_id)
                keggBrite_descs.append(keggBrite_desc)
            items.append(delimiter.join(keggBrite_ids))
            items.append(delimiter.join(keggBrite_descs))
            #
            keggOrthology_ids = []
            for keggOrthology_id, temp in keggAnnoDic[myGene]['keggOrthology'].items():
                keggOrthology_ids.append(keggOrthology_id)
            items.append(delimiter.join(keggOrthology_ids))
            #
            keggEnzyme_ids = []
            for keggEnzyme_id, temp in keggAnnoDic[myGene]['keggEnzyme'].items():
                keggEnzyme_ids.append(keggEnzyme_id)
            items.append(delimiter.join(keggEnzyme_ids))
            #
            out.write('{0}\n'.format('\t'.join(items)))
    out.close()

def keggTableMaker(keggTableDic, db_table_out):
    out = open(db_table_out, 'w')
    for hierarchy_a, bDic in keggTableDic.items():
        for hierarchy_b, keggBriteDic in bDic.items():
            for keggBrite_id, keggBrite_desc in keggBriteDic.items():
                out.write('{0}\t{1}@{2}\t{3}\n'.format(keggBrite_id, hierarchy_a, hierarchy_b, keggBrite_desc))
    out.close()

def main(args):
    print(args)
    #
    keggAnnoDic, keggTableDic = kaas_q1_parser(args.kaas_q1)
    keggUpdate(keggAnnoDic, args.genes_in, args.genes_out, 'Order', 'GeneId', ',')
    keggTableMaker(keggTableDic, args.db_table_out)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-kq', '--kaas-q1', default='q00001.keg')
    parser.add_argument('-gi', '--genes-in', default='genes.xls')
    parser.add_argument('-go', '--genes-out', default='genes.kegg.xls')
    parser.add_argument('-dt', '--db-table-out', default='kegg.db.table')
    args = parser.parse_args()
    main(args)
