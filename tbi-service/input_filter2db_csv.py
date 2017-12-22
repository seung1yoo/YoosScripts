import os,sys

try:
        input_stat=sys.argv[1]

except:
        print 'python input_filter2db_csv.py input/filter2db > output.txt'
        exit(1)

for item in open(input_stat).read().split('//')[1:]:
	buf=list()
	name=''
	for rec in item.strip().split('\n'):
		if rec.startswith('NAME'):
			name=rec.strip().split()[2]
			pass
		elif rec.startswith('RAW.N_READ'):
			buf.append(('RAW.N_READ',rec.strip().split()[1]))
			pass
		elif rec.startswith('RAW.N_SELECT'):
			buf.append(('RAW.N_SELECT',rec.strip().split()[1]))
			pass
		elif rec.startswith('RAW.N_UNSELECT'):
			buf.append(('RAW.N_UNSELECT',rec.strip().split()[1]))
			pass
		elif rec.startswith('RAW.N_BASE'):
			buf.append(('RAW.N_BASE',rec.strip().split()[1]))
			pass

		for title,value in buf:
			print "%s,%s,%s"%(name,title,value)
		
		buf=list()
'''
for line in open(input_stat).split('//'):
	if line.startswith('NAME'):
	elif 
        if line.startswith('ID'):
                header=line.strip().split('\t')
        else:
                rec=line.strip().split('\t')
                sID,sNAME,sPAIR,sBASE,sBASE_Q30,sBASE_Q20,sREAD,sREAD_Q30,sREAD_Q20,sBASE_GC,sBASE_AT,sBASE_NN=rec

                if sPAIR=='READ_1':
                        print "%s,RAW.N_BASE.ALL.READ_1,%s"%(sID,sBASE)
                        print "%s,RAW.N_BASE.Q30.READ_1,%s"%(sID,sBASE_Q30)
                        print "%s,RAW.N_BASE.Q20.READ_1,%s"%(sID,sBASE_Q20)
                        print "%s,RAW.N_READ.ALL.READ_1,%s"%(sID,sREAD)
                        print "%s,RAW.N_READ.Q30.READ_1,%s"%(sID,sREAD_Q30)
                        print "%s,RAW.N_READ.Q20.READ_1,%s"%(sID,sREAD_Q20)
                        print "%s,RAW.GC_BASE.READ_1,%s"%(sID,sBASE_GC)
                        print "%s,RAW.AT_BASE.READ_1,%s"%(sID,sBASE_AT)
                        print "%s,RAW.N_BASE.READ_1,%s"%(sID,sBASE_NN)
                        READ_LEN=int(sBASE)/int(sREAD)
                        print "%s,RAW.LEN_AVG.READ_1,%s"%(sID,READ_LEN)
                        print "%s,RAW.LEN_MIN.READ_1,%s"%(sID,READ_LEN)
                        print "%s,RAW.LEN_MAX.READ_1,%s"%(sID,READ_LEN)
                else:
                        print "%s,RAW.N_BASE.ALL.READ_2,%s"%(sID,sBASE)
                        print "%s,RAW.N_BASE.Q30.READ_2,%s"%(sID,sBASE_Q30)
                        print "%s,RAW.N_BASE.Q20.READ_2,%s"%(sID,sBASE_Q20)
                        print "%s,RAW.N_READ.ALL.READ_2,%s"%(sID,sREAD)
                        print "%s,RAW.N_READ.Q30.READ_2,%s"%(sID,sREAD_Q30)
                        print "%s,RAW.N_READ.Q20.READ_2,%s"%(sID,sREAD_Q20)
                        print "%s,RAW.GC_BASE.READ_2,%s"%(sID,sBASE_GC)
                        print "%s,RAW.AT_BASE.READ_2,%s"%(sID,sBASE_AT)
                        print "%s,RAW.N_BASE.READ_2,%s"%(sID,sBASE_NN)
                        READ_LEN=int(sBASE)/int(sREAD)
                        print "%s,RAW.LEN_AVG.READ_2,%s"%(sID,READ_LEN)
                        print "%s,RAW.LEN_MIN.READ_2,%s"%(sID,READ_LEN)
                        print "%s,RAW.LEN_MAX.READ_2,%s"%(sID,READ_LEN)
'''