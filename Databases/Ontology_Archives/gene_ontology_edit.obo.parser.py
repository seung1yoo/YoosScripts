def term_iter(input):
    block = []
    for line in open(input):
        if line.startswith('[Term]'):
            yield block
            block = [line.strip()]
        else:
            block.append(line.strip())
    yield block

def main(args):
    goDic = dict()
    term_count = 0
    for block in term_iter(args.input):
        term_count += 1
        alt_count = 0
        if term_count in [1]:
            continue
        for line in block:
            if line.startswith('[Term]'):
                continue
            elif line.startswith('[Typedef]'):
                continue
            elif line.startswith('id:'):
                acc = line.split('id:')[-1].strip()
                goDic.setdefault(acc, {})
            elif line.startswith('name:'):
                goDic[acc].setdefault('name', line.split('name:')[-1].strip())
            elif line.startswith('namespace:'):
                goDic[acc].setdefault('namespace', line.split('namespace:')[-1].strip())
            elif line.startswith('alt_id:'):
                alt_count += 1
                acc_alt = line.split('alt_id:')[-1].strip()
                goDic.setdefault(acc_alt, {})
                goDic[acc_alt].setdefault('name', '{0}_alt{1}'.format(goDic[acc]['name'], alt_count))
                goDic[acc_alt].setdefault('namespace', '{0}'.format(goDic[acc]['namespace']))
                continue
            elif line.startswith('def:'):
                #goDic[acc].setdefault('def', line.split('def:')[-1].strip())
                continue
            elif line.startswith('synonym:'):
                continue
            elif line.startswith('is_a:'):
                continue
            elif line.startswith('subset:'):
                continue
            elif line.startswith('xref:'):
                continue
            elif line.startswith('comment:'):
                continue
            elif line.startswith('is_obsolete:'):
                continue
            elif line.startswith('consider:'):
                continue
            elif line.startswith('relationship:'):
                continue
            elif line.startswith('replaced_by:'):
                continue
            elif line.startswith('transitive_over:'):
                continue
            elif line.startswith('is_metadata_tag:'):
                continue
            elif line.startswith('is_class_level:'):
                continue
            elif line.startswith('is_transitive:'):
                continue
            elif line.startswith('holds_over_chain:'):
                continue
            elif not line.strip():
                continue
            else:
                print (line)
                import sys
                sys.exit()
    print (len(goDic))

    out = open(args.output, 'w')
    for go, infoDic in goDic.items():
        out.write('{0}\n'.format('\t'.join([go, infoDic['namespace'], infoDic['name']])))
    out.close()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='ftp://ftp.geneontology.org/go/ontology-archive//gene_ontology_edit.obo.2016-06-01.gz',
                        default='gene_ontology_edit.obo.2020-08-01')
    parser.add_argument('-o', '--output', default='gene_ontology_edit.obo.2020-08-01.table')
    args = parser.parse_args()
    main(args)
