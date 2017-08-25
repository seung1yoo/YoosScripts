
def sample_in_parser(samplein):
    samDic = dict()
    for line in open(samplein):
        items = line.rstrip('\n').split()
        if line.startswith('INPUT'):
            samNum = items[1]
            samNum = "S{0:0>4}".format(samNum)
            samName = '_'.join(items[4].split('/')[-1].split('_')[:-1])
            print samName
            samDic.setdefault(samName, samNum)
    return samDic

def main(args):
    samDic = sample_in_parser(args.samplein)
    #
    aDic = dict()
    for alog in args.inputs:
        samName = alog.split('/')[-1].split('.log')[0]
        for line in open(alog):
            items = line.rstrip('\n').split('\t')
            if line.startswith('N_RAW_READ'):
                aDic.setdefault(samDic[samName], {}).setdefault('RAW.N_READ', items[1])
            elif line.startswith('N_SELECT'):
                aDic.setdefault(samDic[samName], {}).setdefault('RAW.N_SELECT', items[1])
            elif line.startswith('N_UNSELECT'):
                aDic.setdefault(samDic[samName], {}).setdefault('RAW.N_UNSELECT', items[1])
            elif line.startswith('N_BASEPAIR'):
                aDic.setdefault(samDic[samName], {}).setdefault('RAW.N_BASE', items[2])
            else:
                continue
    #
    for sample, infoDic in sorted(aDic.iteritems()):
        for info, value in infoDic.iteritems():
            print "INSERT INTO samples_info VALUES('{0}','{1}','{2}');".format(sample, info, value)

if __name__=='__main__':
    import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', nargs='+', default=glob.glob('./*.log'))
    parser.add_argument('-s', '--samplein', default='../sample.denovo.in')
    args = parser.parse_args()
    main(args)
