import glob, pdb

Files = glob.glob('pre_*.xls')

for File in Files :
    
    fr = open('%s' % File, 'r')
        
    lst = []

    for line in fr.xreadlines() :
        words       = line.rstrip('\n').split('\t')
        #
        infos       = words[1:11] # GeneId ~ Desc
        #
        exps_1      = words[18:21] # Timeseries 1 
        exps_2      = words[15:18] # Timeseries 2
        exps_3      = words[12:15] # Timeseries 3
        #exps_4      = words[]
        #
        degs_1 = words[21:25] # DEG 1
        degs_2 = words[29:33] # DEG 2
        #
        gos = [words[-1]]

        sub_lst = infos
        sub_lst.extend(exps_1)
        sub_lst.extend(exps_2)
        sub_lst.extend(exps_3)
        sub_lst.extend(degs_1)
        sub_lst.extend(degs_2)
        sub_lst.extend(gos)

        lst.append(sub_lst)

    fw = open('%s' % '_'.join(File.split('_')[1:]), 'w')

    for item in lst :
        fw.write('%s\n' % '\t'.join(item))
    
    fw.close()

