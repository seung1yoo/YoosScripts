
class DEFIND_SEGMENT:
    def __init__(self, args):
        self.in_fn = args.in_fn
        self.out_fn = args.out_fn
        self.idxs = [int(x) for x in args.idxs]
        #
        out = open(self.out_fn, 'w')
        #out.write('{0}\n'.format('\t'.join(['#CHROM','START','END'
        for segment_s in self.iter_segment():
            if not segment_s:
                continue
            chrom = segment_s[0][0]
            start = segment_s[0][1]
            end   = segment_s[-1][1]
            snp_count = len(segment_s)
            pattern = self.grep_pattern(segment_s[0])
            #print chrom, start, end, snp_count, pattern
            out.write('{0}\n'.format('\t'.join([chrom, start, end, str(snp_count), pattern])))
        out.close()

    def iter_segment(self):
        segments = list()
        pre_pattern = ''
        for line in open(self.in_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#CHROM']:
                print 'sample order = {0}'.format(self.grep_pattern(items))
                continue
            pattern = self.grep_pattern(items)
            #
            if not pattern in [pre_pattern]:
                yield segments
                segments = [items]
            else:
                segments.append(items)
            pre_pattern = pattern
        yield segments

    def grep_pattern(self, items):
        patterns = []
        for idx in self.idxs:
            origin = items[idx].split('=')[-1]
            patterns.append(origin)
        return '@'.join(patterns)

def main(args):
    ds = DEFIND_SEGMENT(args)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_fn', default='TBD171076_Glycine2_DiffHomo.1.1-HWANGGEUM_vs_2-DAEPUNG.v2.xls')
    parser.add_argument('--out_fn', default='TBD171076_Glycine2_DiffHomo.1.1-HWANGGEUM_vs_2-DAEPUNG.v2.segment.bed')
    parser.add_argument('--idxs', nargs='+', default=[16,17,18,19,20,21])
    args = parser.parse_args()
    main(args)
