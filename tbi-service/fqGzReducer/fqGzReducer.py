
import os

class REDUCER:
    def __init__(self, conf, outdir):
        self.conf_fn = conf
        self.outdir = outdir
        #
        self.conf_dic = dict()
        self.conf2dic()
        #
        self.f_dic = dict()

    def conf2dic(self):
        for line in open(self.conf_fn):
            items = line.rstrip('\n').split()
            sample_name = items[0]
            target_read_count = items[1]
            self.conf_dic.setdefault(sample_name, target_read_count)

    def f_finder(self, tag_1='_1', tag_2='_2', f_ext='.fq.gz'):
        for sample_name, target_read_count in self.conf_dic.iteritems():
            for tag in [tag_1, tag_2]:
                f_name = '{0}{1}{2}'.format(sample_name, tag, f_ext)
                if os.path.exists(f_name):
                    self.f_dic.setdefault(sample_name, []).append(f_name)
                else:
                    print('cannot found {0}'.format(f_name))

    def make_cmd(self):
        cmds = []
        for sample_name, fn_s in self.f_dic.iteritems():
            for fn in fn_s:
                cmd = ['zcat']
                cmd.append(fn)
                cmd.append('|')
                cmd.append('head -n')
                cmd.append(str(int(self.conf_dic[sample_name])*4))
                cmd.append('/dev/stdin')
                cmd.append('|')
                cmd.append('gzip -c')
                cmd.append('>')
                cmd.append('{0}/{1}'.format(self.outdir.rstrip('/'),fn))
                cmds.append(' '.join(cmd))
        return cmds


def main(args):
    instance = REDUCER(args.conf, args.outdir)
    instance.f_finder()
    cmds = instance.make_cmd()
    for cmd in cmds:
        print(cmd)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf', help="[SAMPLE_NAME] [NUM_READ/each for&rev]",
        default='fqGzReducer.conf')
    parser.add_argument('--outdir', help="outdir",
        default='../')
    args = parser.parse_args()
    main(args)
