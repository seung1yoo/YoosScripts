
def fileFinder(path, fileName):
    cmd = "find {0} -iname '{1}'".format(path, fileName)
    files = subprocess.check_output(cmd, shell=True)
    fileDic = dict()
    for afile in files.split('\n'):
        if not afile:
            continue
        sample_name = afile.split('/')[-2]
        fileDic.setdefault(sample_name, afile)
    return fileDic

def align_summary_parser(afile):
    infoDic = dict()
    for line in open(afile):
        if not line.strip():
            break
        items = line.strip().split()
        if items[1] in ['reads:']:
            tag = items[0]
            infoDic.setdefault(tag, {})
        elif items[1] in [':']:
            strand = items[0]
            value = items[2]
            infoDic[tag].setdefault(strand, value)
        elif items[1] in ['these:']:
            if len(items) in [10]:
                hit2 = items[2]
                hit20 = items[7].lstrip('(')
                infoDic[tag].setdefault('HitOver2', hit2)
                infoDic[tag].setdefault('HitOver20', hit20)
            elif len(items) in [11]:
                hit2 = items[2]
                hit20 = items[8].lstrip('(')
                infoDic[tag].setdefault('HitOver2', hit2)
                infoDic[tag].setdefault('HitOver20', hit20)
        elif items[1] in ['overall']:
            overall = items[0]

    for tag, typeDic in infoDic.iteritems():
        rate = int(typeDic['Mapped']) / float(typeDic['Input']) * 100
        typeDic.setdefault('Mapped_Rate', str(round(rate,2)))
        rate = int(typeDic['HitOver2']) / float(typeDic['Input']) * 100
        typeDic.setdefault('HitOver2_Rate', str(round(rate,2)))
        rate = int(typeDic['HitOver20']) / float(typeDic['Input']) * 100
        typeDic.setdefault('HitOver20_Rate', str(round(rate,2)))

    return infoDic, overall


def main(args):
    fileDic = fileFinder(args.path, 'align_summary.txt')
    '''
    tags = ['Left', 'Right', 'Unpaired']
    types = ['Input', 'Mapped', 'Mapped_Rate', 'HitOver2', 'HitOver2_Rate', 'HitOver20', 'HitOver20_Rate']
    titles = ['Sample']
    titles.extend(['{0}_{1}'.format(tag, type) for tag in tags for type in types])
    titles.append('Overall_rate')
    out = open(args.outFile , 'w')
    out.write('{0}\n'.format('\t'.join(titles)))
    '''
    out = open(args.outFile , 'w')
    tags = ['Left', 'Right', 'Unpaired']
    types = ['Input', 'Mapped', 'Mapped_Rate', 'HitOver2', 'HitOver2_Rate', 'HitOver20', 'HitOver20_Rate']
    titles_1 = ['']
    for tag in tags:
        titles_1.append(tag)
        for n in range(len(types)-1):
            titles_1.append('')
    out.write('{0}\n'.format('\t'.join(titles_1)))
    titles = ['Sample']
    titles.extend(['{0}'.format(typev) for tag in tags for typev in types])
    titles.append('Overall_rate')
    out.write('{0}\n'.format('\t'.join(titles)))
    for sample_name, afile in fileDic.iteritems():
        print sample_name, afile
        infoDic, overall = align_summary_parser(afile)
        out_items = [sample_name]
        for tag in tags:
            for type in types:
                out_items.append(infoDic[tag][type])
        out_items.append(overall)
        out.write('{0}\n'.format('\t'.join(out_items)))
    out.close()

if __name__=='__main__':
    import subprocess
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', default='/BiO/BioProjects/TBD160530-Centipede-RNAref-20160909/06.mapping')
    parser.add_argument('-o', '--outFile', default='align_summary_parser.py.result')
    args = parser.parse_args()
    main(args)
