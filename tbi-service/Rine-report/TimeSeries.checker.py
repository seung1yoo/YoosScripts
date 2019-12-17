import glob
import os


def ts_output_checker(path, ts_id):
    pattern_dic = dict()
    for fn_path in glob.glob(os.path.join(path, '*')):
        fn = fn_path.split('/')[-1]
        items = fn.split('.')
        pattern = items[1]
        if pattern in ['design', 'xls']:
            continue
        if items[0].startswith('TS'):
            pattern_dic.setdefault(pattern, {}).setdefault('GO.R',0)
            pattern_dic.setdefault(pattern, {}).setdefault('GO.wall',0)
    #
    for pattern, info_dic in pattern_dic.iteritems():
        if os.path.exists(os.path.join(path, '{0}.{1}.GO.R'.format(ts_id, pattern))):
            pattern_dic[pattern]['GO.R'] = 1
        if os.path.exists(os.path.join(path, '{0}.{1}.GO.wall.xls'.format(ts_id, pattern))):
            pattern_dic[pattern]['GO.wall'] = 1
    #
    return pattern_dic


def main(args):
    pattern_dic = ts_output_checker(args.path, args.ts_id)
    for pattern, info_dic in pattern_dic.iteritems():
        if info_dic['GO.R'] and info_dic['GO.wall']:
            print '{0} results are existed'.format(pattern)
        if not info_dic['GO.R'] and not info_dic['GO.wall']:
            print 'There are not genes in {0}'.format(pattern)
        if info_dic['GO.R'] and not info_dic['GO.wall']:
            print '{0} have to re-run'.format(pattern)
            cmd = '{0} CMD BATCH {1}'.format(args.r_path,
                    os.path.join(args.path, '{0}.{1}.GO.R'.format(args.ts_id, pattern)))
            print cmd
            os.system(cmd)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', default='./')
    parser.add_argument('--ts-id', default='TS001')
    parser.add_argument('--r-path', default='/export/home/siyoo/miniconda2/envs/goseq/bin/R')
    args = parser.parse_args()
    main(args)
