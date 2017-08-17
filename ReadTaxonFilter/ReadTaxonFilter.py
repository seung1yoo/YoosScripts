

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import glob
import os

class ReadTaxonFilter():
    def __init__(self, sample, reads, annoFile, taxons, outDir):
        self.sample = sample
        self.taxons = [x.upper() for x in taxons]
        print '# list of filter out taxons : {0}'.format(\
                ','.join(self.taxons))
        #
        self.annoFile = self.annoFile_checker(annoFile)
        self.cleanReadDic = self.annoFile_parser()
        #
        self.readDic = self.reads_checker(reads)
        self.reads_filter(outDir)
        #

    def reads_filter(self, outDir):
        if not os.path.isdir(outDir):
            os.mkdir(outDir)
        #
        for idx, aDic in self.readDic.iteritems():
            file_name, file_ext = os.path.splitext(os.path.split(aDic['file'])[-1])
            outFile = '{0}/{1}.clean{2}'.format(outDir, file_name, file_ext)
            print '# reads_filter : {0} make OK'.format(outFile)
            #
            with gzip.open(outFile, 'wb') as out_handle:
                total_n = 0
                filter_n = 0
                clean_n = 0
                for title, seq, qual in FastqGeneralIterator(aDic['handle']):
                    total_n += 1
                    readId = title.split()[0]
                    if readId in self.cleanReadDic:
                        out_handle.write("@{0}\n{1}\n+\n{2}\n".format(title, seq, qual))
                        clean_n += 1
                    else:
                        filter_n += 1
                print '# reads_filter : {3} : clean / total * 100 = {0}/{1}*100 = {2}%'.format(\
                        clean_n, total_n, clean_n/float(total_n)*100, self.sample)
            #

    def reads_checker(self, reads):
        readDic = dict()
        for idx, read in enumerate(reads):
            if os.path.isfile(read):
                print '# reads_checker : {0} OK'.format(read)
            else:
                print '# reads_checker : {0} FAIL'.format(read)
                import sys
                sys.exit()
            #
            if read.endswith('.gz'):
                print '# reads_checker : {0} (format = gz)'.format(read)
                readDic.setdefault(idx, {}).setdefault('file', read)
                readDic.setdefault(idx, {}).setdefault('handle', gzip.open(read, 'r'))
                readDic.setdefault(idx, {}).setdefault('type', 'gz')
                #
            elif read.endswith('.fastq') or read.endswith('.fq'):
                print '# reads_checker : {0} (format = fastq)'.format(read)
                readDic.setdefault(idx, {}).setdefault('file', read)
                readDic.setdefault(idx, {}).setdefault('handle', open(read, 'r'))
                readDic.setdefault(idx, {}).setdefault('type', 'fastq')
                #
            elif read.endswith('.fasta') or read.endswith('.fa'):
                print '# reads_checker : {0} (format = fasta)'.format(read)
                readDic.setdefault(idx, {}).setdefault('file', read)
                readDic.setdefault(idx, {}).setdefault('handle', open(read, 'r'))
                readDic.setdefault(idx, {}).setdefault('type', 'fasta')
                #
            else:
                print '# reads_checker : {0} (format = FAIL)'.format(read)
                import sys
                sys.exit()
        #
        return readDic

    def annoFile_checker(self, annoFile):
        if os.path.isfile(annoFile):
            print '# annoFile_checker : {0} OK'.format(annoFile)
        else:
            print '# annoFile_checker : {0} FAIL'.format(annoFile)
            import sys
            sys.exit()
        #
        with open(annoFile, 'r') as annoFile_handle:
            file_format_check = 1 # negative control
            file_title_check = 0 # positive control
            for line in annoFile_handle:
                if line.startswith('#queryID'):
                    file_title_check = 1
                #
                items = line.rstrip('\n').split('\t')
                if not len(items) in [12]:
                    print items
                    file_format_check = 0
            #
            if file_format_check in [0] or file_title_check in [0]:
                print '# annoFile_checker : {0} FORMAT_FAIL'.format(annoFile)
                import sys
                sys.exit()
        #
        return annoFile

    def annoFile_parser(self):
        cleanReadDic = dict()
        with open(self.annoFile, 'r') as annoHandle:
            for line in annoHandle:
                items = line.rstrip('\n').split('\t')
                if line.startswith('#queryID'):
                    idxDic = dict()
                    for idx, item in enumerate(items):
                        idxDic.setdefault(item, idx)
                    continue
                #
                readId = items[idxDic['#queryID']]
                filterOut = 0
                for taxon in items[idxDic['TaxName']].split(';'):
                    if taxon.upper() in self.taxons:
                        filterOut = 1
                #
                if filterOut in [0]:
                    cleanReadDic.setdefault(readId, None)
                #
        print '# annoFile_parser : No. clean read {1}'.format(self.annoFile, len(cleanReadDic))
        return cleanReadDic



def main(args):
    rtf = ReadTaxonFilter(args.sample, args.reads, args.annoFile, args.taxons, args.outDir)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample', default='PO-0-1_1')
    parser.add_argument('-o', '--outDir', default='ReadTaxonFilter')
    parser.add_argument('-r', '--reads', nargs='+',
            default=['/BiO/BioProjects/TBD160494-SSU-Metagenome-20170609/00.RawData/PO-0-1_1.fastq.gz',
                     '/BiO/BioProjects/TBD160494-SSU-Metagenome-20170609/00.RawData/PO-0-1_2.fastq.gz'])
    parser.add_argument('-a', '--annoFile', help='Taxonomy profiling result',
            default='/BiO/BioProjects/TBD160494-SSU-Metagenome-20170609/01.remove_plantea/PO-0-1_1/Scaffold_taxonomy_Annotation.xls')
    parser.add_argument('-t', '--taxons', nargs='+',
            default=['Viridiplantae'])
    args = parser.parse_args()
    main(args)
