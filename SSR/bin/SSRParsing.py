import sys, os

if len(sys.argv) != 5 :
	print "Usage : python %s <ssr.result> <reference.stat> <SSR.Mer.Freq.outname> <SSR.Motif.Freq.outname>" % sys.argv[0]
	sys.exit()

SSR_Class = {}

f = open(sys.argv[2], 'r')

TotalSize = 0

for lines in f.xreadlines() :
	if lines[0] == '#' :
		continue
	else :
		words = lines.rstrip().split("\t")
		TotalSize = int(words[1])

f.close()

f = open(sys.argv[1], 'r')

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if words[0] == 'ID' :
		continue
	if SSR_Class.has_key(int(words[2])) :
		if SSR_Class[int(words[2])].has_key(words[3].upper()) :
			SSR_Class[int(words[2])][words[3].upper()] += 1
		else :
			SSR_Class[int(words[2])][words[3].upper()] = 1
	else :
		SSR_Class[int(words[2])] = {}
		SSR_Class[int(words[2])][words[3].upper()] = 1

f.close()

w1 = open(sys.argv[3], 'w')
w2 = open(sys.argv[4], 'w')

w1.write("#Repeat Type\tFrequency\tFrequency per million\n")
w2.write("#Motif\tFrequency\tFrequency per million\n")
for MER in sorted(SSR_Class.keys()) :
	TotalFreq = 0
	for SSR in sorted(SSR_Class[MER].keys()) :
		TotalFreq += SSR_Class[MER][SSR]
		w2.write("%s\t%d\t%.2f\n" % (SSR, SSR_Class[MER][SSR], (SSR_Class[MER][SSR] * 1000000.0 / TotalSize)))
	w1.write("%d\t%d\t%.2f\n" % (MER, TotalFreq, (TotalFreq * 1000000.0 / TotalSize)))

w1.close()
w2.close()
