
import os

class TargetGO():
    def __init__(self):
        self.smrnaDic = dict()
        self.isdegDic = dict()
        self.tgDic = dict()
        self.tgAnnoDic = dict()

    def targetPredictionDEG_parser(self, input):
        for line in open(input): 
            items = line.rstrip('\n').split('\t')
            if line.startswith('miRNA'):
                idxDic = dict()
                idxDic.setdefault('miRNA', items.index('miRNA'))
                idxDic.setdefault('estimatedDEG', items.index('estimatedDEG'))
                continue
            miRNA = items[idxDic['miRNA']]
            isdeg = items[idxDic['estimatedDEG']]
            self.isdegDic.setdefault(miRNA, isdeg)
            self.smrnaDic.setdefault(miRNA, {})

    def addTarget(self, database_files):
        for dbfile in database_files:
            for line in open(dbfile):
                items = line.rstrip('\n').split('\t')
                if line.startswith('ID'):
                    idxDic = dict()
                    idxDic.setdefault('smrna', items.index('SMRNAID'))
                    idxDic.setdefault('gene', items.index('GENESYMBOL'))
                    continue
                smrna = items[idxDic['smrna']]
                gene = items[idxDic['gene']]
                if smrna in self.smrnaDic:
                    self.tgDic.setdefault(gene, {})
                    self.smrnaDic[smrna].setdefault(gene, {})
                    #break # for test

    def addGeneInfo(self, genesInfo):
        for line in open(genesInfo):
            items = line.rstrip('\n').split('\t')
            if line.startswith('#StableID'):
                idxDic = dict()
                idxDic.setdefault('GeneName', items.index('GeneName'))
                idxDic.setdefault('GO', items.index('GO'))
                idxDic.setdefault('GeneStart', items.index('GeneStart'))
                idxDic.setdefault('GeneEnd', items.index('GeneEnd'))
                continue
            geneName = items[idxDic['GeneName']]
            if geneName in self.tgDic:
                go = items[idxDic['GO']]
                start = int(items[idxDic['GeneStart']])
                end = int(items[idxDic['GeneEnd']])
                posis = sorted([start, end], reverse=True)
                geneLen = posis[0] - posis[1] + 1
                self.tgDic[geneName].setdefault('go', go)
                self.tgDic[geneName].setdefault('geneLen', geneLen)
                self.tgAnnoDic.setdefault(geneName, 0)
        #
        for smrna, geneDic in self.smrnaDic.iteritems():
            for gene, aDic in geneDic.iteritems():
                self.smrnaDic[smrna][gene].update(self.tgDic[gene])
          
    def makeTable(self, outprefix):
        self.tableFile = '{0}.targetGoTable'.format(outprefix)
        out = open(self.tableFile, 'w')
        titles = ['#smrna','isDEG','targetGeneSymbol','targetGeneLen','targetGeneGO']
        out.write('{0}\n'.format('\t'.join([str(x) for x in titles])))
        for smrna, geneDic in self.smrnaDic.iteritems():
            if not geneDic:
                items = [smrna, self.isdegDic[smrna], 'NotFound', 'NotFound', 'NotFound']
                out.write('{0}\n'.format('\t'.join([str(x) for x in items])))
            else:
                for gene, infoDic in geneDic.iteritems():
                    if 'geneLen' in infoDic:
                        items = [smrna, self.isdegDic[smrna], gene, infoDic['geneLen'], infoDic['go']]
                        out.write('{0}\n'.format('\t'.join([str(x) for x in items])))
                    else:
                        items = [smrna, self.isdegDic[smrna], gene, 'NotFound', 'NotFound']
                        out.write('{0}\n'.format('\t'.join([str(x) for x in items])))
        out.close()

    def goseqInputMaker(self, outprefix):
        goseqDic = dict()
        for line in open(self.tableFile):
            items = line.rstrip('\n').split('\t')
            if line.startswith('#smrna'):
                idxDic = dict()
                idxDic.setdefault('gene', items.index('targetGeneSymbol'))
                idxDic.setdefault('isdeg', items.index('isDEG'))
                idxDic.setdefault('len', items.index('targetGeneLen'))
                idxDic.setdefault('go', items.index('targetGeneGO'))
                continue
            gene = items[idxDic['gene']]
            isdeg = items[idxDic['isdeg']]
            geneLen = items[idxDic['len']]
            go = items[idxDic['go']]
            if go in ['NotFound']:
                continue
            else:
                goseqDic.setdefault(gene, {}).setdefault('isdeg', isdeg)
                goseqDic.setdefault(gene, {}).setdefault('len', geneLen)
                goseqDic.setdefault(gene, {}).setdefault('go', go)
        #
        self.goseq_genes = '{0}.goseq.deg.xls'.format(outprefix)
        self.goseq_go = '{0}.goseq.mapping.xls'.format(outprefix)
        self.goseq_genes_length = '{0}.goseq.len.xls'.format(outprefix)
        #
        out_genes = open(self.goseq_genes, 'w')
        out_go = open(self.goseq_go, 'w')
        out_genes_length = open(self.goseq_genes_length, 'w')
        for gene, infoDic in goseqDic.iteritems():
            out_genes.write('{0}\t{1}\n'.format(gene, infoDic['isdeg']))
            out_go.write('{0}\t{1}\n'.format(gene, infoDic['go']))
            out_genes_length.write('{0}\t{1}\n'.format(gene, infoDic['len']))
        out_genes.close()
        out_go.close()
        out_genes_length.close()
        
    def goseq_executor(self, outprefix):
        goseq_code = """
library(goseq)

Input_genes<-read.table("{0}", head=F, sep="\t")
Input_go<-read.table("{1}", head=F, sep="\t")
Input_genes_length<-read.table("{2}", head=F, sep="\t")

genes<-as.integer(Input_genes[,2])
names(genes)<-Input_genes[,1]
go<-strsplit(as.character(Input_go[,2]), ",")
names(go)<-Input_go[,1]
length<-as.numeric(Input_genes_length[,2])

pwf=nullp(genes, bias.data=length, plot.fit=F)

GO.wall<-goseq(pwf, gene2cat=go, method="Wallenius")
write.table(GO.wall, "{3}.goseq.wall.xls", quote=F, sep="\t")

GO.nobias<-goseq(pwf, gene2cat=go, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Hypergeometric")
write.table(GO.nobias, "{3}.goseq.obias.xls", quote=F, sep="\t")

GO.samp<-goseq(pwf, gene2cat=go, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Sampling", repcnt=2000)
write.table(GO.samp, "{3}.goseq.sampling.xls", quote=F, sep="\t")
""".format(self.goseq_genes, self.goseq_go, self.goseq_genes_length, outprefix)
        #
        goseq_script = '{0}.goseq.R'.format(outprefix)
        out = open(goseq_script, 'w')
        out.write(goseq_code)
        out.close()

        cmd = "/BiO/BioPeople/smrna/analysis_pack/Tools/R-3.2.4/bin/R CMD BATCH --no-save --no-resotre {0} {0}out".format(goseq_script)
        print cmd
        os.system(cmd)

def main(args):
    targetgo = TargetGO()
    targetgo.targetPredictionDEG_parser(args.input)
    print '#No.smrna : {0}'.format(len(targetgo.smrnaDic))
    targetgo.addTarget(args.databases)
    print '#No.Target gene : {0}'.format(len(targetgo.tgDic))
    targetgo.addGeneInfo(args.genesInfo)
    print '#No.Annotated Target gene : {0}'.format(len(targetgo.tgAnnoDic))
    targetgo.makeTable(args.outprefix)
    print '#Target GENE,GENE_LEN,GO Table : {0}'.format(targetgo.tableFile)
    print '#if not found target gene and information(len, go) -> NotFound tagging in Table'
    #
    targetgo.goseqInputMaker(args.outprefix)
    targetgo.goseq_executor(args.outprefix)
    

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='small RNA TargetPrediction_DEG *.out file',
        default='TargetPrediction_DEG_001/DE.DEG_001.targetscan.out')
    parser.add_argument('-db', '--databases', nargs='+', help='small RNA Target Gene Database files',
        default=['/BiO/BioPeople/smrna/analysis_pack/targetgo/database/db.mirdb.cnv.refined.txt',
                 '/BiO/BioPeople/smrna/analysis_pack/targetgo/database/db.targetscanhuman.cnv.refined.txt'])
    parser.add_argument('-g', '--genesInfo', help='Gene Information files',
        default='/BiO/BioPeople/smrna/analysis_pack/targetgo/input/genes.ENS72.info')
    parser.add_argument('-o', '--outprefix', help='prefix of result file',
        default='TargetPrediction_DEG_001/DE.DEG_001.targetscan.out')
    args = parser.parse_args()
    main(args)




