## Parse line started with "D"
def D_line(line) :
    lines = line.split(';')
    Entry_Name = lines[1].strip().split() 
    return get_result(lines, Entry_Name)

## Parse line started with "E" 
def E_line(line) :
    if ';' in line :
        lines = line.split(';')
        Entry_Name = lines[1].strip().split()
        return get_result(lines, Entry_Name)


## Get result from line 
def get_result(lines, Entry_Name):
    result_ = {}
    result_['GeneID'] = lines[0].split()[1]
    if len(lines) == 2:
        result_['Entry'] = Entry_Name[0]; result_['Name'] = " ".join(Entry_Name[1:])
        result_['Def'] = '-'; result_['EC'] = '-'
    
    elif len(lines) == 3 :
        if len(Entry_Name) == 2 :
            result_['Entry'] = Entry_Name[0]; result_['Name'] = Entry_Name[1]
        else :
            result_['Entry'] = Entry_Name[0]; result_['Name'] = " ".join(Entry_Name[1:])
    
        if '[EC:' in lines[2] :
            Def_EC = re.search(r'([\w\W\s]+)\s+(\[[\w\W\s]+\])',lines[2].strip())
            result_['Def'] = Def_EC.group(1); result_['EC'] = Def_EC.group(2)
        else :
            Def = re.search(r'[\w\W\s]+', lines[2].strip()).group()
            result_['Def'] = Def; result_['EC'] = '-'
    return result_

## Parse q00001.keg and q00003.keg files
def parser_q1_q3(q1fn, fo) :
    check_A = 'off'; check_B = 'off'; check_C = 'off';
    for line in open(q1fn) :
        if line.startswith('+') :
            Keg = line.strip().split()[1]
            if Keg == 'KO' :
                Keg = 'Orthology'
            elif Keg == 'RModule' :
                Keg = 'ReModule'
            continue
        elif line.startswith('#') or line.startswith('!') or len(line.strip()) == 1 :
            continue

        if line.startswith('A') :
            check_A = 'on'
            cat_A = line.strip()[1:]
        elif check_A == 'on' and line.startswith('B') :
            check_B = 'on'
            cat_B = re.search(r'B\s+([\w\W\s]+)', line.strip()).group(1)
        elif check_A == 'on' and check_B == 'on' and line.startswith('C') :
            check_C = 'on'
            cat_C = re.search(r'C\s+([\w\W\s]+)', line.strip()).group(1)
        elif check_A == 'on' and check_B == 'on' and check_C == 'on' and line.startswith('D') :
            cat_D = '-'
            result = D_line(line.strip())
            fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(Keg, result['GeneID'], result['Entry'], result['Name'], result['Def'], result['EC'], cat_A.replace("<b>", "").replace("</b>", ""), cat_B, cat_C, cat_D))

## Parse q00002.keg file
def parser_q2(q2fn, fo) :
    check_A = 'off'; check_B = 'off'; check_C = 'off'; check_D = 'off';
    for line in open(q2fn) :
        if line.startswith('+') :
            Keg = line.strip().split()[1]
            if Keg == 'Module' :
                Keg = 'Modules'
            continue
        elif line.startswith('#') or line.startswith('!') or len(line.strip()) == 1 :
            continue

        if line.startswith('A') :
            check_A = 'on'
            cat_A = line.strip()[1:]
        elif check_A == 'on' and line.startswith('B') :
            check_B = 'on'
            cat_B = re.search(r'B\s+([\w\W\s]+)', line.strip()).group(1)
        elif check_A == 'on' and check_B == 'on' and line.startswith('C') :
            check_C = 'on'
            cat_C = re.search(r'C\s+([\w\W\s]+)', line.strip()).group(1)
        elif check_A == 'on' and check_B == 'on' and check_C == 'on' and line.startswith('D') :
            check_D = 'on'
            cat_D = re.search(r'D\s+([\w\W\s]+)', line.strip()).group(1)
        elif check_A == 'on' and check_B == 'on' and check_C == 'on' and check_D == 'on' and line.startswith('E') :
            result = E_line(line.strip())
            if result :
                fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(Keg, result['GeneID'], result['Entry'], result['Name'], result['Def'], result['EC'],cat_A.replace("<b>", "").replace("</b>", ""), cat_B.replace("<b>", "").replace("</b>", ""), cat_C, cat_D))

## Find genes not existed in keg files
def Genes_not_in_KegFile(fo) :
    gene_in_fasta = []
    gene_in_keg = []

    header = os.popen("grep '>' C_sinensis_plantkingdomgdb_1.pep.fa").readlines()
    gene_q = os.popen("grep '\.1' q0000* | awk '{print substr($2,1,11)}'").readlines()
    for each in header :
        gene = re.search(r'>([\w.]+)', each.strip()).group(1)
        gene_in_fasta.append(gene)
    for each in gene_q :
        gene_in_keg.append(each.strip())
    for i in gene_in_fasta :
        if i not in gene_in_keg :
            fo.write("-\t%s\t-\t-\t-\t-\t-\t-\t-\t-\n" %(i))
 
## Sort result file by column, "GeneID" 
def file_sort(fn) :
    fi = pd.read_table(fn, sep='\t')
    fi_sort = fi.sort_values(['GeneID'])
    fi_sort.to_csv(fn+".sorted.xls", sep='\t', index=False)


def main(args) :
    fo = open(args.outfn, 'a')
    fo.write("KegFileName\tGeneID\tEntry\tName\tDefinition\tEC:Num\tCategory_A\tCategory_B\tCategory_C\tCategory_D\n")
    Genes_not_in_KegFile(fo)
    parser_q1_q3(args.q3fn, fo)
    parser_q1_q3(args.q1fn, fo)
    parser_q2(args.q2fn, fo)
    fo.close()
    file_sort(args.outfn)



if __name__ == '__main__' :
    import argparse
    import re
    import os
    import pandas as pd
    parser = argparse.ArgumentParser()
    parser.add_argument('-q1', '--q1fn', help = 'q1 file name')
    parser.add_argument('-q2', '--q2fn', help = 'q2 file name')
    parser.add_argument('-q3', '--q3fn', help = 'q3 file name')
    parser.add_argument('-o', '--outfn', help = 'output file name')
    args = parser.parse_args()
    main(args)
