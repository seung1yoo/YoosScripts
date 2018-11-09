




def table_to_dic(infn):
    a_dic = dict()
    for line in open(infn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['ID']:
            idx_dic = dict()
            for idx, item in enumerate(items):
                idx_dic.setdefault(item, idx)
            print items
            print idx_dic
            continue
        _id = items[idx_dic['ID']]
        _strand = items[idx_dic['']]
        for colname, idx in idx_dic.iteritems():
            if colname in ['N_BASE']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_BASE.ALL.{0}'.format(_strand), items[idx_dic['N_BASE']])
            elif colname in ['N_BASE_Q30']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_BASE.Q30.{0}'.format(_strand), items[idx_dic['N_BASE_Q30']])
            elif colname in ['N_BASE_Q20']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_BASE.Q20.{0}'.format(_strand), items[idx_dic['N_BASE_Q20']])
            elif colname in ['N_READ']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_READ.ALL.{0}'.format(_strand), items[idx_dic['N_READ']])
            elif colname in ['N_READ_Q30']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_READ.Q30.{0}'.format(_strand), items[idx_dic['N_READ_Q30']])
            elif colname in ['N_READ_Q20']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_READ.Q20.{0}'.format(_strand), items[idx_dic['N_READ_Q20']])
            elif colname in ['N_BASE_GC']:
                a_dic.setdefault(_id, {}).setdefault('RAW.GC_BASE.{0}'.format(_strand), items[idx_dic['N_BASE_GC']])
            elif colname in ['N_BASE_AT']:
                a_dic.setdefault(_id, {}).setdefault('RAW.AT_BASE.{0}'.format(_strand), items[idx_dic['N_BASE_AT']])
            elif colname in ['N_BASE_NN']:
                a_dic.setdefault(_id, {}).setdefault('RAW.N_BASE.{0}'.format(_strand), items[idx_dic['N_BASE_NN']])
            else:
                continue
    return a_dic



def main(args):
    a_dic = table_to_dic(args.infn)

    out = open(args.outfn, 'w')
    for _id, key_dic in a_dic.iteritems():
        for key, value in key_dic.iteritems():
            out.write('{0}\n'.format(','.join([_id, key, value])))
    out.close()

    print '''Do like this'''
    print '''cd [PATH of Rine]/db'''
    print '''sqlite3 rine.db'''
    print '''sqlite> .headers on'''
    print '''sqlite> .separator ","'''
    print '''sqlite> .import {0} samples_info'''.format(args.outfn)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infn', default='stat_fastq.summary')
    parser.add_argument('--outfn', default='stat_fastq.summary.csv')
    args = parser.parse_args()
    main(args)
