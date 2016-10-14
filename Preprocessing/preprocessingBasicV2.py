


def main(args):
    cmd = 'java -jar {0} PE '\
            '-threads 10 '\
            '{2} {3} '\
            '{4}.Clean.1P.fq.gz '\
            '{4}.Clean.1U.fq.gz '\
            '{4}.Clean.2P.fq.gz '\
            '{4}.Clean.2U.fq.gz '\
            'ILLUMINACLIP:{1}:2:30:10 '\
            'SLIDINGWINDOW:4:20 '\
            'LEADING:3 '\
            'TRAILING:3 '\
            'MINLEN:36 '\
            ''.format(args.program, args.adapter, args.forward, args.reverse, args.sample)
    print cmd
    os.system(cmd)

if __name__=='__main__':
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--program',
            default = '/BiO/BioPeople/siyoo/00.Tools/Trimmomatic-0.33/trimmomatic-0.33.jar')
    parser.add_argument('-a', '--adapter',
            default = '/BiO/BioPeople/siyoo/00.Tools/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa')
    parser.add_argument('-f', '--forward')
    parser.add_argument('-r', '--reverse')
    parser.add_argument('-s', '--sample', default='mySampleFiltered.fq.gz')
    args = parser.parse_args()
    main(args)
