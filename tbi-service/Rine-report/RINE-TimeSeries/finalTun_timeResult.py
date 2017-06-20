import glob, pdb

Files = glob.glob('pre_*.xls')

for File in Files :
    
    fr = open('%s' % File, 'r')
        
    lst = []

    for line in fr.xreadlines() :
        words       = line.rstrip('\n').split('\t')
        geneId      = words[1]
        geneAcc     = words[2]
        geneName    = words[3]
        geneDesc    = words[10]
        exp_1       = words[12]
        exp_2       = words[13]
        exp_3       = words[14]
        deg_1_pv = words[16]
        deg_1_fc = words[18]
        deg_2_pv = words[24]
        deg_2_fc = words[26]
        go = words[-1]

        sub_lst = [geneId, geneAcc, geneName, geneDesc,
                   exp_1, exp_2, exp_3,
                   deg_1_pv, deg_1_fc,
                   deg_2_pv, deg_2_fc,
                   go]

        lst.append(sub_lst)

    fw = open('%s' % '_'.join(File.split('_')[1:]), 'w')

    for item in lst :
        fw.write('%s\n' % '\t'.join(item))
    
    fw.close()

