def adapterTrim(files_pair, outhandle, adapters):
    adapterTrimSeq = '{0}_nonadapter.fq'.format(outhandle)
    if not os.path.isfile(adapterTrimSeq):
        cmd = 'clc_adapter_trim -r -i {0} -g {1} '\
                               '-a {2} > {1}.err'.format(' '.join(files_pair), 
                                       adapterTrimSeq, ' -a '.join(adapters))
        print cmd
        os.system(cmd)
    else:
        print '{0} file exist'.format(adapterTrimSeq)
    return adapterTrimSeq

def qualityTrim(adapterTrimSeq):
    reportfile = '{0}.q20.report'.format(adapterTrimSeq.split('.fq')[0])
    outfile_pair = '{0}.q20.pair.fq'.format(adapterTrimSeq.split('.fq')[0])
    outfile_single = '{0}.q20.single.fq'.format(adapterTrimSeq.split('.fq')[0])
    if not os.path.isfile(reportfile):
        #for nextSeq 150x0.6 = 90bp
        #for HiSeq 101x0.9 = 90bp
        cmd = 'clc_quality_trim -r {0} -f 33 -c 20 -m 90 '\
                               '-l 0.9 -b 0.1 -p {1} '\
                               '-o {2} > {3}'.format(adapterTrimSeq, outfile_pair, 
                                                     outfile_single, reportfile)
        print cmd
        os.system(cmd)
    else:
        print '{0} file exist'.format(reportfile)
    return outfile_pair, outfile_single

def main(outdir, sample, adapters, files_pair):
    if os.path.isdir(outdir):
        pass
    else:
        os.system('mkdir {0}'.format(outdir))
    outhandle = '{0}/{1}'.format(outdir, sample)
    adapterTrimSeq = adapterTrim(files_pair, outhandle, adapters)
    qualityTrimPair, qualityTrimSingle = qualityTrim(adapterTrimSeq)

if __name__ == "__main__":
    import glob
    import os
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Paired-end preprocessing')
    parser.add_argument('-o', '--out-dir', type=str)
    parser.add_argument('-s', '--sample-name', type=str)
    parser.add_argument('-a', '--adapters', nargs='+', default=[
                   'GATCGGAAGAGC']) # TruSeq nano
    parser.add_argument('-l', '--libraries-pair', nargs='+', help='-l '\
            '/tier2/sample1_R1.fq /tier2/sample1_R2.fq')
    args = parser.parse_args()
    main(args.out_dir, args.sample_name, args.adapters, args.libraries_pair)
