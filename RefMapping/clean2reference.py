class CLEAN2REF():
    def __init__(self, configure, program, gff, ref):
        massage = '''
        This Source code was made for Reference mapping.
        CONDITIONs : only paired-end data, non-strand-specific data
        A. you have to prepare configure file. (-c, --configure)
        there are several colums (comma-separated)
         a. No.Sample
         b. Sample Name(id)
         c. Paired Tag (1 or 2)
         d. Path of clean data (full)

        B. Tool (-p, --program) - I will update...
         a. Tophat2

        C. Indexed Reference - I will update with "B.Tool"
         a. bowtie2-build
        '''
        print massage
        self.configure = self.file_checker(configure)
        self.cleanDic = self.configure_parser(configure)
        self.program = self.program_selector(program)
        self.gff = self.file_checker(gff)
        self.ref = self.file_checker('{0}.1.bt2'.format(ref))
        if self.ref:
            self.ref = ref

    def file_checker(self, source):
        if os.path.isfile(source):
            print '#FILE CHECK : OKAY : {0}'.format(source)
        else:
            print '#FILE CHECK : FAIL : {0}'.format(source)
            print '#STOP process.'
            sys.exit()
        return source

    def configure_parser(self, configure):
        cleanDic = dict()
        for line in open(configure):
            if line.startswith('#'):
                continue
            items = line.strip().split(',')
            cleanDic.setdefault(items[1], {}).setdefault(items[2], items[3])
        for sample, tagDic in cleanDic.iteritems():
            for tag, clean_file in tagDic.iteritems():
                print '#CONFIGURE : {0} : {1} : {2}'.format(sample, tag, clean_file)
                self.file_checker(clean_file)
        return cleanDic

    def program_selector(self, program):
        program_path_dic = {
                'Tophat2':'/BiO/BioTools/TopHat/tophat-2.0.13.Linux_x86_64/tophat'
                }
        print '#PROGRAM : {0}'.format(program_path_dic[program])
        return program_path_dic[program]

    def bowtie2_exe(self, sample, tagDic):
        cmd = '{0} -o {1} -p 30 -r 300 --mate-std-dev 200 --library-type fr-unstranded -G {2} '\
                '{3} {4},{5} {6},{7}'.format(self.program, sample, self.gff,
                 self.ref, tagDic['1P'], tagDic['1U'], tagDic['2P'], tagDic['2U'])
        print cmd
        os.system(cmd)

def main(args):
    clean2ref = CLEAN2REF(args.configure, args.program, args.gff, args.ref)
    for sample, tagDic in clean2ref.cleanDic.iteritems():
        print sample, tagDic
        if args.program in ['Tophat2']:
            clean2ref.bowtie2_exe(sample, tagDic)

if __name__=='__main__':
    import os
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--configure', default='configure_example.csv')
    parser.add_argument('-p', '--program', choices=['Tophat2'], default='Tophat2')
    parser.add_argument('-f', '--gff', default='/BiO/BioPeople/siyoo/01.FNP-Strawberry-Transcriptome-201604/00.Ref/ncbi/ref_FraVesHawaii_1.0_top_level.keep.gff3')
    parser.add_argument('-r', '--ref', default='/BiO/BioPeople/siyoo/01.FNP-Strawberry-Transcriptome-201604/00.Ref/ncbi/101020_ref_FraVesHawaii_1.0.all.rename')
    args = parser.parse_args()
    main(args)

