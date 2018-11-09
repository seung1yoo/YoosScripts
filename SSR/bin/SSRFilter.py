import sys, os

if len(sys.argv) != 3 :
	print "python %s <in.ssr> <out.ssr>" % sys.argv[0]
	sys.exit()

f = open(sys.argv[1], 'r')
w = open(sys.argv[2], 'w')
w.write("ID\tSSR_Number\tMotif_Length\tMotif\tRepeats\tStart\tEnd\tSeq_Length\n")

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if len(words[3]) == 4 :
		if words[3][:2] == words[3][2:] :
			continue
		else :
			w.write(lines)
	elif len(words[3]) == 6 :
		if words[3][:2] == words[3][2:4] and words[3][:2] == words[3][4:] :
			continue
		elif words[3][:3] == words[3][3:] :
			continue
		else :
			w.write(lines)
	elif len(words[3]) == 8 :
		if words[3][:2] == words[3][2:4] and words[3][:2] == words[3][4:6] and words[3][:2] == words[3][6:] :
			continue
		elif words[3][:4] == words[3][4:] :
			continue
		else :
			w.write(lines)
	elif len(words[3]) == 9 :
		if words[3][:3] == words[3][3:6] and words[3][:3] == words[3][6:] :
			continue
		else :
			w.write(lines)
	elif len(words[3]) == 10 :
		if words[3][:2] == words[3][2:4] and words[3][:2] == words[3][4:6] and words[3][:2] == words[3][6:8] and words[3][:2] == words[3][8:] :
			continue
		elif words[3][:5] == words[3][5:] :
			continue
		else :
			w.write(lines)
	else :
		w.write(lines)

f.close()
w.close()
