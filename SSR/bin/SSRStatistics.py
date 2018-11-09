import sys

if len(sys.argv) != 3 :
	print "python %s <ssr.freq> <di-tri>" % sys.argv[0]
	sys.exit()

f = open (sys.argv[1], 'r')

len2Total = 0
len3Total = 0
AC = 0
AG = 0
AT = 0
CG = 0
AAT = 0
AAG = 0
AAC = 0
ATG = 0
AGT = 0
AGG = 0
AGC = 0
ACG = 0
ACC = 0
GGC = 0

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if lines[0] == '#' or words[0] == 'ID':
		continue
	if len(words[0]) == 2 :
		len2Total += int(words[1])
		if words[0] == 'AC' or words[0] == 'CA' or words[0] == 'TG' or words[0] == 'GT' :
			AC += int(words[1])
		elif words[0] == 'AG' or words[0] == 'GA' or words[0] == 'TC' or words[0] == 'CT' :
			AG += int(words[1])
		elif words[0] == 'AT' or words[0] == 'TA' :
			AT += int(words[1])
		elif words[0] == 'CG' or words[0] == 'GC' :
			CG += int(words[1])
	elif len(words[0]) == 3 :
		len3Total += int(words[1])
		if words[0] == 'AAT' or words[0] == 'ATA' or words[0] == 'TAA' or words[0] == 'ATT' or words[0] == 'TTA' or words[0] == 'TAT' :
			AAT += int(words[1])
		elif words[0] == 'AAG' or words[0] == 'AGA' or words[0] == 'GAA' or words[0] == 'CTT' or words[0] == 'TTC' or words[0] == 'TCT' :
			AAG += int(words[1])
		elif words[0] == 'AAC' or words[0] == 'ACA' or words[0] == 'CAA' or words[0] == 'GTT' or words[0] == 'TTG' or words[0] == 'TGT' :
			AAC += int(words[1])
		elif words[0] == 'ATG' or words[0] == 'TGA' or words[0] == 'GAT' or words[0] == 'CAT' or words[0] == 'ATC' or words[0] == 'TCA' :
			ATG += int(words[1])
		elif words[0] == 'AGT' or words[0] == 'GTA' or words[0] == 'TAG' or words[0] == 'ACT' or words[0] == 'CTA' or words[0] == 'TAC' :
			AGT += int(words[1])
		elif words[0] == 'AGG' or words[0] == 'GGA' or words[0] == 'GAG' or words[0] == 'CCT' or words[0] == 'CTC' or words[0] == 'TCC' :
			AGG += int(words[1])
		elif words[0] == 'AGC' or words[0] == 'GCA' or words[0] == 'CAG' or words[0] == 'GCT' or words[0] == 'CTG' or words[0] == 'TGC' :
			AGC += int(words[1])
		elif words[0] == 'ACG' or words[0] == 'CGA' or words[0] == 'GAC' or words[0] == 'CGT' or words[0] == 'GTC' or words[0] == 'TCG' :
			ACG += int(words[1])
		elif words[0] == 'ACC' or words[0] == 'CCA' or words[0] == 'CAC' or words[0] == 'GGT' or words[0] == 'GTG' or words[0] == 'TGG' :
			ACC += int(words[1])
		elif words[0] == 'GGC' or words[0] == 'GCG' or words[0] == 'CGG' or words[0] == 'GCC' or words[0] == 'CCG' or words[0] == 'CGC' :
			GGC += int(words[1])

f.close()

w = open(sys.argv[2], 'w')

w.write("Motif\tFrequency\t%%\n")
w.write("AC (AC/CA/TG/GT)\t%d\t%.2f%%\n" % (AC, (AC * 100.0 / len2Total)))
w.write("AG (AG/GA/TC/CT)\t%d\t%.2f%%\n" % (AG, (AG * 100.0 / len2Total)))
w.write("AT (AT/TA)\t%d\t%.2f%%\n" % (AT, (AT * 100.0 / len2Total)))
w.write("CG (CG/GC)\t%d\t%.2f%%\n" % (CG, (CG * 100.0 / len2Total)))
w.write("AAT (AAT/ATA/TAA/ATT/TTA/TAT)\t%d\t%.2f%%\n" % (AAT, (AAT * 100.0 / len3Total)))
w.write("AAG (AAG/AGA/GAA/CTT/TTC/TCT)\t%d\t%.2f%%\n" % (AAG, (AAG * 100.0 / len3Total)))
w.write("AAC (AAC/ACA/CAA/GTT/TTG/TGT)\t%d\t%.2f%%\n" % (AAC, (AAC * 100.0 / len3Total)))
w.write("ATG (ATG/TGA/GAT/CAT/ATC/TCA)\t%d\t%.2f%%\n" % (ATG, (ATG * 100.0 / len3Total)))
w.write("AGT (AGT/GTA/TAG/ACT/CTA/TAC)\t%d\t%.2f%%\n" % (AGT, (AGT * 100.0 / len3Total)))
w.write("AGG (AGG/GGA/GAG/CCT/CTC/TCC)\t%d\t%.2f%%\n" % (AGG, (AGG * 100.0 / len3Total)))
w.write("AGC (AGC/GCA/CAG/GCT/CTG/TGC)\t%d\t%.2f%%\n" % (AGC, (AGC * 100.0 / len3Total)))
w.write("ACG (ACG/CGA/GAC/CGT/GTC/TCG)\t%d\t%.2f%%\n" % (ACG, (ACG * 100.0 / len3Total)))
w.write("ACC (ACC/CCA/CAC/GGT/GTG/TGG)\t%d\t%.2f%%\n" % (ACC, (ACC * 100.0 / len3Total)))
w.write("GGC (GGC/GCG/CGG/GCC/CCG/CGC)\t%d\t%.2f%%\n" % (GGC, (GGC * 100.0 / len3Total)))
w.write("Di-mer Total\t%d\n" % len2Total)
w.write("Tri-mer Total\t%d\n" % len3Total)

w.close()
