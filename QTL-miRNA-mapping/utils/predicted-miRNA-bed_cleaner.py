import os
import glob

def main():
    beds = glob.glob('./*/*.bed')
    for bed in beds:
        print bed
        sample = bed.split('/')[-2]
        print sample
        out = open('{0}.mirdeep.result.temp.bed'.format(sample), 'w')
        for line in open(bed):
            out.write('Chr.{0}'.format(line))
        out.close()
        #
        cmd = 'bedtools sort -i {0}.mirdeep.result.temp.bed > {0}.mirdeep.result.clean.bed'.format(sample)
        print cmd
        os.system(cmd)
        #




if __name__=='__main__':
    main()
