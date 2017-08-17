
class BLAST_CHECKER():
    def __init__(self, inblast, infasta, tempFile):
        self.inblast = inblast
        self.infasta = infasta
        self.tempFile = tempFile
        self.last_title = self.last_title_finder()
        self.process_rate = self.process_checker()

    def last_title_finder(self):
        cmd = 'tail -n900000 {0} | grep title | tail -n1 > {1}'.format(self.inblast, self.tempFile)
        #print cmd
        os.system(cmd)
        #
        for line in open(self.tempFile):
            if line.strip().startswith('title'):
                last_title = line.strip().split('"')[-1].split()[0]
        #
        if not last_title:
            print("There is no title.")
            import sys
            sys.exit()
        #
        return last_title

    def process_checker(self):
        process_rate = 0.0
        process_n = 0
        process_total = 0
        #
        for record in SeqIO.parse(open(self.infasta), 'fasta'):
            process_total += 1
            #
            if self.last_title in record.id:
                process_n = process_total
        #
        if process_n in [0]:
            print("Title {0} in not found".format(self.last_title))
            import sys
            sys.exit()
        #
        process_rate = round(process_n / float(process_total) * 100, 2)
        #
        return process_rate

def main(args):
    bc = BLAST_CHECKER(args.inblast, args.infasta, args.tempFile)
    print("Last Title of blast : {0}".format(bc.last_title))
    print("Process Rate of {1} : {0}".format(bc.process_rate, bc.inblast))

if __name__=='__main__':
    from Bio import SeqIO
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ib', '--inblast', help='11 format of blast in RINE process')
    parser.add_argument('-if', '--infasta', help='fasta of blastx')
    parser.add_argument('-t', '--tempFile', default='blast_checker.temp')
    args = parser.parse_args()
    main(args)
