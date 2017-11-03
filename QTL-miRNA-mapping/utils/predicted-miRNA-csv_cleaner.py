import glob

def sample_name_parser(csv):
    sample = csv.split('/')[-2]
    return sample

def file_preprocessing(sample, csv):
    out_novel = open('{0}.mirdeep.result.novel.xls'.format(sample), 'w')
    out_known = open('{0}.mirdeep.result.known.xls'.format(sample), 'w')
    novel = 0
    known = 0
    for line in open(csv):
        if not line.rstrip('\n'):
            novel = 0
            known = 0
            continue
        #
        if line.startswith('novel miRNAs predicted by miRDeep2'):
            novel = 1
            known = 0
        elif line.startswith('mature miRBase miRNAs detected by miRDeep2'):
            novel = 0
            known = 1
        elif line.startswith('miRDeep2 score'):
            novel = 0
            known = 0
        elif line.startswith('#miRBase miRNAs detected by miRDeep2'):
            novel = 0
            known = 0
        #
        if novel in [1] and known in [0]:
            out_novel.write(line)
        elif novel in [0] and known in [1]:
            out_known.write(line)
        else:
            continue

    out_novel.close()
    out_known.close()
    return '{0}.mirdeep.result.novel.xls'.format(sample), '{0}.mirdeep.result.known.xls'.format(sample)


def main():
    csvs = glob.glob('./*/*.csv')
    for csv in csvs:
        sample = sample_name_parser(csv)
        print sample
        #
        novel_fn, known_fn = file_preprocessing(sample, csv)
        print novel_fn
        print known_fn




if __name__=='__main__':
    main()
