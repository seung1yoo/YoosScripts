import os
import glob

def main(args):
    #beds = glob.glob('./*/*.bed')
    for idx, bed in enumerate(args.beds):
        print bed
        #sample = bed.split('/')[-2]
        sample = args.samples[idx]
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
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-bs', '--beds', nargs='+', help='miRDeep2 results (bed format)')
    parser.add_argument('-ss', '--samples', nargs='+', help='sample names')
    args = parser.parse_args()
    main(args)
