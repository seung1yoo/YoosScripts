import os

def main(args):
    for filename in args.inputfiles:
        print filename
        filenames = os.path.splitext(filename)
        up_outfile = '{0}.{1}.UP{2}'.format(filenames[0], args.degcriteria, filenames[1])
        dn_outfile = '{0}.{1}.DOWN{2}'.format(filenames[0], args.degcriteria, filenames[1])
        #
        up_outhandle = open(up_outfile, 'w')
        dn_outhandle = open(dn_outfile, 'w')
        titles = []
        for line in open(filename):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['miRNA']:
                idxDic = dict()
                for idx, item in enumerate(items):
                    if item in ['m.value', 'p.value', 'q.value', 'estimatedDEG']:
                        idxDic.setdefault(item, idx)
                up_outhandle.write(line)
                dn_outhandle.write(line)
                continue
            m_value = float(items[idxDic['m.value']])
            p_value = float(items[idxDic['p.value']])
            q_value = float(items[idxDic['q.value']])
            estimatedDEG = int(items[idxDic['estimatedDEG']])
            deg_is_tag = 0
            deg_ud_tag = 'updown'
            if args.degcriteria in ['TCC']:
                if estimatedDEG in [1] and m_value > 0:
                    deg_is_tag = 1
                    deg_ud_tag = 'up'
                elif estimatedDEG in [1] and m_value < 0:
                    deg_is_tag = 1
                    deg_ud_tag = 'down'
                else:
                    deg_is_tag = 0
                    deg_ud_tag = 'updown'
            elif args.degcriteria in ['FC2P005']:
                if p_value < 0.05 and m_value >= 1:
                    deg_is_tag = 1
                    deg_ud_tag = 'up'
                elif p_value < 0.05 and m_value <= -1:
                    deg_is_tag = 1
                    deg_ud_tag = 'down'
                else:
                    deg_is_tag = 0
                    deg_ud_tag = 'updown'
            #
            if deg_is_tag in [1] and deg_ud_tag in ['up']:
                up_outhandle.write(line)
            elif deg_is_tag in [1] and deg_ud_tag in ['down']:
                dn_outhandle.write(line)
        up_outhandle.close()
        dn_outhandle.close()

if __name__=='__main__':
    import glob
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ifs', '--inputfiles', nargs='+',
            default=glob.glob('/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD161022-TBD161075-TBD170259-TBD170428_HUMC_Mouse_SmallRNA_Report_20171218/differential_expression/DE.DEG_*.targetscan.xls'))
    parser.add_argument('-dc', '--degcriteria', choices=('FC2P005', 'TCC'),
            default='FC2P005')
    args = parser.parse_args()
    main(args)
