import os

def term_checker(category, termDic):
    for term in termDic:
        if term in category:
            return 1
    return 0

def qtl_parser(originBED, listFile, trait, sub_trait, release):
    termDic = dict()
    for line in open(listFile):
        termDic.setdefault(line.strip(), '')
    #
    cateDic = dict()
    out = open('cattle-QTLdb.{0}.{1}.{2}.bed'.format(release, trait, sub_trait), 'w')
    for line in open(originBED):
        if line.startswith('#'):
            continue
        items = line.rstrip().split('\t')
        name = items[3]
        category = '('.join(name.split('(')[:-1]).strip()
        #
        start = items[1]
        end = items[2]
        if not start or not end:
            continue
        if start in ['0'] and end in ['100']:
            continue
        #
        if term_checker(category, termDic) and category.endswith('QTL'):
            out.write(line)
        #
        cateDic.setdefault(category, {}).setdefault(name, 0)
        cateDic[category][name] += 1
    out.close()
    cmd = 'bedtools sort -i cattle-QTLdb.{0}.{1}.{2}.bed > cattle-QTLdb.{0}.{1}.{2}.sort.bed'.format(release, trait, sub_trait)
    print cmd
    os.system(cmd)
    #
    out = open('cattle-QTLdb.{0}.{1}.{2}.stats.xls'.format(release, trait, sub_trait), 'w')
    for category, nameDic in cateDic.iteritems():
        if term_checker(category, termDic) and category.endswith('QTL'):
            out.write('{0}\n'.format('\t'.join([category, '@', str(len(nameDic))])))
    out.close()
    #

def main(args):
    
    qtl_parser(args.originBED, args.traitListFn, args.trait, args.sub_trait, args.release)

    # # Regacy command for DEV.
    #qtl_parser(originBED, 'Reproduction_Traits.All.QTLterms', 'Reproduction', 'All')
    #qtl_parser(originBED, 'Reproduction_Traits.Fertility.QTLterms', 'Reproduction', 'Fertility')
    #qtl_parser(originBED, 'Reproduction_Traits.General.QTLterms', 'Reproduction', 'General')
    #qtl_parser(originBED, 'Reproduction_Traits.Reproductive-hormone-level.QTLterms', 'Reproduction', 'Reproductive-hormone-level')
    #qtl_parser(originBED, 'Reproduction_Traits.Semen-quality.QTLterms', 'Reproduction', 'Semen-quality')
    # #
    #qtl_parser(originBED, 'Milk_Traits.All.QTLterms', 'Milk', 'All')
    #qtl_parser(originBED, 'Milk_Traits.Milk_composition_fat.QTLterms', 'Milk', 'Milk_composition_fat')
    #qtl_parser(originBED, 'Milk_Traits.Milk_composition_other.QTLterms', 'Milk', 'Milk_composition_other')
    #qtl_parser(originBED, 'Milk_Traits.Milk_composition_protein.QTLterms', 'Milk', 'Milk_composition_protein')
    #qtl_parser(originBED, 'Milk_Traits.Milk_processing_trait.QTLterms', 'Milk', 'Milk_processing_trait')
    #qtl_parser(originBED, 'Milk_Traits.Milk_yield.QTLterms', 'Milk', 'Milk_yield')
    # #

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ob', '--originBED', help='cattle QTL information as BED format',
            default='cattle-QTLdb.Release33.bed')
    parser.add_argument('-tl', '--traitListFn', help='trait list file name',
            default='Milk_Traits.All.QTLterms')
    parser.add_argument('-r', '--release', default='Release33')
    parser.add_argument('-t', '--trait', help='Trait category',
            default='Milk')
    parser.add_argument('-st', '--sub-trait', help='Sub-Trait name',
            default='All')
    args = parser.parse_args()
    main(args)

