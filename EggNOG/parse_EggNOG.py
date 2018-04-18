class EggNOG():
    def __init__(self, ifns, ofx):
        self.ifns = sorted(ifns)
        self.ofx = ofx
        #
        self.ome_fn = self.make_ome()
        self.stats_fn = self.make_stats()

    def make_ome(self):
        titles = ['query_name', 'seed_eggNOG_ortholog', 'seed_ortholog_evalue',
                  'seed_ortholog_score', 'predicted_gene_name', 'GO_terms',
                  'KEGG_KO', 'BiGG_Reactions', 'Annotation_tax_scope', 'Matching_OGs',
                  'best_OG|evalue|score', 'COG functional categories', 'eggNOG_HMM_model_annotation']
        ofn = '{0}.ome.xls'.format(self.ofx)
        fh = open(ofn, 'w')
        fh.write('{0}\n'.format('\t'.join(titles)))
        for ifn in self.ifns:
            for line in open(ifn):
                fh.write(line)
        fh.close()
        return ofn

    def make_stats(self):
        #http://clovr.org/docs/clusters-of-orthologous-groups-cogs/
        classDic = {'D':'CELLULAR PROCESSES AND SIGNALING',
                    'M':'CELLULAR PROCESSES AND SIGNALING',
                    'N':'CELLULAR PROCESSES AND SIGNALING',
                    'O':'CELLULAR PROCESSES AND SIGNALING',
                    'T':'CELLULAR PROCESSES AND SIGNALING',
                    'U':'CELLULAR PROCESSES AND SIGNALING',
                    'V':'CELLULAR PROCESSES AND SIGNALING',
                    'W':'CELLULAR PROCESSES AND SIGNALING',
                    'Y':'CELLULAR PROCESSES AND SIGNALING',
                    'Z':'CELLULAR PROCESSES AND SIGNALING',
                    'A':'INFORMATION STORAGE AND PROCESSING',
                    'B':'INFORMATION STORAGE AND PROCESSING',
                    'J':'INFORMATION STORAGE AND PROCESSING',
                    'K':'INFORMATION STORAGE AND PROCESSING',
                    'L':'INFORMATION STORAGE AND PROCESSING',
                    'C':'METABOLISM',
                    'E':'METABOLISM',
                    'F':'METABOLISM',
                    'G':'METABOLISM',
                    'H':'METABOLISM',
                    'I':'METABOLISM',
                    'P':'METABOLISM',
                    'Q':'METABOLISM',
                    'R':'POORLY CHARACTERIZED',
                    'S':'POORLY CHARACTERIZED',
                    'NoHit':'NoHit'}
        clusterDic = {'D':'Cell cycle control, cell division, chromosome partitioning',
                      'M':'Cell wall/membrane/envelope biogenesis',
                      'N':'Cell motility',
                      'O':'Post-translational modification, protein turnover, and chaperones',
                      'T':'Signal transduction mechanisms',
                      'U':'Intracellular trafficking, secretion, and vesicular transport',
                      'V':'Defense mechanisms',
                      'W':'Extracellular structures',
                      'Y':'Nuclear structure',
                      'Z':'Cytoskeleton',
                      'A':'RNA processing and modification',
                      'B':'Chromatin structure and dynamics',
                      'J':'Translation, ribosomal structure and biogenesis',
                      'K':'Transcription',
                      'L':'Replication, recombination and repair',
                      'C':'Energy production and conversion',
                      'E':'Amino acid transport and metabolism',
                      'F':'Nucleotide transport and metabolism',
                      'G':'Carbohydrate transport and metabolism',
                      'H':'Coenzyme transport and metabolism',
                      'I':'Lipid transport and metabolism',
                      'P':'Inorganic ion transport and metabolism',
                      'Q':'Secondary metabolites biosynthesis, transport, and catabolism',
                      'R':'General function prediction only',
                      'S':'Function unknown',
                      'NoHit':'NoHit'}
        cogDic = dict()
        for line in open(self.ome_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['query_name']:
                idxDic = dict()
                for idx, item in enumerate(items):
                    idxDic.setdefault(item, idx)
                continue
            query = items[idxDic['query_name']]
            cog = items[idxDic['COG functional categories']]
            for acog in cog.split(','):
                acog = acog.strip()
                if acog in ['']:
                    acog = 'NoHit'
                cogDic.setdefault(acog, {}).setdefault(query, None)
        #
        ofn = '{0}.stats.xls'.format(self.ofx)
        fh = open(ofn, 'w')
        fh.write('{0}\n'.format('\t'.join(['COG','CLASS','CLUSTER','COUNT'])))
        for acog, queryDic in cogDic.iteritems():
            fh.write('{0}\n'.format('\t'.join([acog, classDic[acog], clusterDic[acog], str(len(queryDic))])))
        fh.close()
        return ofn


def main(args):
    egg = EggNOG(args.ifns, args.ofx)

if __name__=='__main__':
    import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ifns', nargs='+',
        default=glob.glob('*.emapper.annotations'))
    parser.add_argument('-ofx',
        default='EggNOG')
    args = parser.parse_args()
    main(args)

