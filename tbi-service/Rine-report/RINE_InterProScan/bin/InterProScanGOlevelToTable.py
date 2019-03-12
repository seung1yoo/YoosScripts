import sys

if len(sys.argv) != 3 :
	print "Usage : python %s <FileList.txt> <OutPrefix>" % sys.argv[0]
	sys.exit()


def GOLevel2(_go) :
	GOLIST = _go
	LEV2 = []
	while 1 :
		TMPLIST = []
		for GOID in GOLIST :
			if GOUPPER.has_key(GOID) :
				for GOV1 in GOUPPER[GOID] :
					if not GOUPPER.has_key(GOV1) :
						if not GOID in LEV2 :
							LEV2.append(GOID)
					if not GOV1 in TMPLIST :
						TMPLIST.append(GOV1)
		GOLIST = TMPLIST
		if len(GOLIST) == 0 :
			break

	return LEV2

def GOLevel4(_go) :
	GOLIST = _go
	LEV4 = []
	while 1 :
		TMPLIST = []
		for GOID in GOLIST :
			if GOUPPER.has_key(GOID) :
				for GOV3 in GOUPPER[GOID] :
					if GOUPPER.has_key(GOV3) :
						for GOV2 in GOUPPER[GOV3] :
							if GOUPPER.has_key(GOV2) :
								for GOV1 in GOUPPER[GOV2] :
									if not GOUPPER.has_key(GOV1) :
										if not GOID in LEV4 :
											LEV4.append(GOID)
					if not GOV3 in TMPLIST :
						TMPLIST.append(GOV3)
		GOLIST = TMPLIST
		if len(GOLIST) == 0 :
			break

	return LEV4

def GOLevel6(_go) :
	GOLIST = _go
	LEV6 = []
	while 1 :
		TMPLIST = []
		for GOID in GOLIST :
			if GOUPPER.has_key(GOID) :
				for GOV5 in GOUPPER[GOID] :
					if GOUPPER.has_key(GOV5) :
						for GOV4 in GOUPPER[GOV5] :
							if GOUPPER.has_key(GOV4) :
								for GOV3 in GOUPPER[GOV4] :
									if GOUPPER.has_key(GOV3) :
										for GOV2 in GOUPPER[GOV3] :
											if GOUPPER.has_key(GOV2) :
												for GOV1 in GOUPPER[GOV2] :
													if not GOUPPER.has_key(GOV1) :
														if not GOID in LEV6 :
															LEV6.append(GOID)
					if not GOV5 in TMPLIST :
						TMPLIST.append(GOV5)
		GOLIST = TMPLIST
		if len(GOLIST) == 0 :
			break

	return LEV6

def GOLevel8(_go) :
	GOLIST = _go
	LEV8 = []
	while 1 :
		TMPLIST = []
		for GOID in GOLIST :
			if GOUPPER.has_key(GOID) :
				for GOV7 in GOUPPER[GOID] :
					if GOUPPER.has_key(GOV7) :
						for GOV6 in GOUPPER[GOV7] :
							if GOUPPER.has_key(GOV6) :
								for GOV5 in GOUPPER[GOV6] :
									if GOUPPER.has_key(GOV5) :
										for GOV4 in GOUPPER[GOV5] :
											if GOUPPER.has_key(GOV4) :
												for GOV3 in GOUPPER[GOV4] :
													if GOUPPER.has_key(GOV3) :
														for GOV2 in GOUPPER[GOV3] :
															if GOUPPER.has_key(GOV2) :
																for GOV1 in GOUPPER[GOV2] :
																	if not GOUPPER.has_key(GOV1) :
																		if not GOID in LEV8 :
																			LEV8.append(GOID)
					if not GOV7 in TMPLIST :
						TMPLIST.append(GOV7)
		GOLIST = TMPLIST
		if len(GOLIST) == 0 :
			break

	return LEV8

def GOLevel10(_go) :
	GOLIST = _go
	LEV10 = []
	while 1 :
		TMPLIST = []
		for GOID in GOLIST :
			if GOUPPER.has_key(GOID) :
				for GOV9 in GOUPPER[GOID] :
					if GOUPPER.has_key(GOV9) :
						for GOV8 in GOUPPER[GOV9] :
							if GOUPPER.has_key(GOV8) :
								for GOV7 in GOUPPER[GOV8] :
									if GOUPPER.has_key(GOV7) :
										for GOV6 in GOUPPER[GOV7] :
											if GOUPPER.has_key(GOV6) :
												for GOV5 in GOUPPER[GOV6] :
													if GOUPPER.has_key(GOV5) :
														for GOV4 in GOUPPER[GOV5] :
															if GOUPPER.has_key(GOV4) :
																for GOV3 in GOUPPER[GOV4] :
																	if GOUPPER.has_key(GOV3) :
																		for GOV2 in GOUPPER[GOV3] :
																			if GOUPPER.has_key(GOV2) :
																				for GOV1 in GOUPPER[GOV2] :
																					if not GOUPPER.has_key(GOV1) :
																						if not GOID in LEV10 :
																							LEV10.append(GOID)
					if not GOV9 in TMPLIST :
						TMPLIST.append(GOV9)
		GOLIST = TMPLIST
		if len(GOLIST) == 0 :
			break

	return LEV10



FILELIST = sys.argv[1]
GOLEV2COUNTTABLE = sys.argv[2] + ".GOLev2.count.xls"
GOLEV2GENESTABLE = sys.argv[2] + ".GOLev2.genes.xls"
GOLEV4COUNTTABLE = sys.argv[2] + ".GOLev4.count.xls"
GOLEV4GENESTABLE = sys.argv[2] + ".GOLev4.genes.xls"
GOLEV6COUNTTABLE = sys.argv[2] + ".GOLev6.count.xls"
GOLEV6GENESTABLE = sys.argv[2] + ".GOLev6.genes.xls"
GOLEV8COUNTTABLE = sys.argv[2] + ".GOLev8.count.xls"
GOLEV8GENESTABLE = sys.argv[2] + ".GOLev8.genes.xls"
GOLEV10COUNTTABLE = sys.argv[2] + ".GOLev10.count.xls"
GOLEV10GENESTABLE = sys.argv[2] + ".GOLev10.genes.xls"
GOEDGEFILE = "/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180756-JNU-Haliotis-RNAref-SmallReport-FC2P005-20181116/addTask/20181128_InterProScan/si_scripts/bin/go.edge.txt"

GOUPPER = {}
GOPROCESS = {}
GODESCIPTION = {}

f = open(GOEDGEFILE, 'r')

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if not GOUPPER.has_key(words[0]) :
		GOUPPER[words[0]] = []
	GOUPPER[words[0]].append(words[1])
	GOPROCESS[words[0]] = words[3]
	GODESCIPTION[words[0]] = words[2]

f.close()

f = open(FILELIST, 'r')

SPECIESLIST = []
SPFILES = {}

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if lines[0] == '#' :
		continue
	SPECIESLIST.append(words[1])
	SPFILES[words[1]] = words[0]

f.close()

LV2COUNTS = {}
LV2GENES = {}
LV4COUNTS = {}
LV4GENES = {}
LV6COUNTS = {}
LV6GENES = {}
LV8COUNTS = {}
LV8GENES = {}
LV10COUNTS = {}
LV10GENES = {}

for SP in SPECIESLIST :
	f = open(SPFILES[SP], 'r')
	for lines in f.xreadlines() :
		if lines[0] == '#' :
			continue
		else :
			words = lines.rstrip().split("\t")
			if len(words) < 12 :
				continue
			if words[11][:2] == 'GO' :
				name = words[11].split("|")
				GOLV2LIST = GOLevel2(name)
				GOLV4LIST = GOLevel4(name)
				GOLV6LIST = GOLevel6(name)
				GOLV8LIST = GOLevel8(name)
				GOLV10LIST = GOLevel10(name)
				for GOI in GOLV2LIST :
					if not LV2GENES.has_key(GOI) :
						LV2COUNTS[GOI] = {}
						LV2GENES[GOI] = {}
					if not LV2GENES[GOI].has_key(SP) :
						LV2COUNTS[GOI][SP] = 0
						LV2GENES[GOI][SP] = []
					if not words[0] in LV2GENES[GOI][SP] :
						LV2GENES[GOI][SP].append(words[0])
						LV2COUNTS[GOI][SP] += 1
				for GOI in GOLV4LIST :
					if not LV4GENES.has_key(GOI) :
						LV4COUNTS[GOI] = {}
						LV4GENES[GOI] = {}
					if not LV4GENES[GOI].has_key(SP) :
						LV4COUNTS[GOI][SP] = 0
						LV4GENES[GOI][SP] = []
					if not words[0] in LV4GENES[GOI][SP] :
						LV4GENES[GOI][SP].append(words[0])
						LV4COUNTS[GOI][SP] += 1
				for GOI in GOLV6LIST :
					if not LV6GENES.has_key(GOI) :
						LV6COUNTS[GOI] = {}
						LV6GENES[GOI] = {}
					if not LV6GENES[GOI].has_key(SP) :
						LV6COUNTS[GOI][SP] = 0
						LV6GENES[GOI][SP] = []
					if not words[0] in LV6GENES[GOI][SP] :
						LV6GENES[GOI][SP].append(words[0])
						LV6COUNTS[GOI][SP] += 1
				for GOI in GOLV8LIST :
					if not LV8GENES.has_key(GOI) :
						LV8COUNTS[GOI] = {}
						LV8GENES[GOI] = {}
					if not LV8GENES[GOI].has_key(SP) :
						LV8COUNTS[GOI][SP] = 0
						LV8GENES[GOI][SP] = []
					if not words[0] in LV8GENES[GOI][SP] :
						LV8GENES[GOI][SP].append(words[0])
						LV8COUNTS[GOI][SP] += 1
				for GOI in GOLV10LIST :
					if not LV10GENES.has_key(GOI) :
						LV10COUNTS[GOI] = {}
						LV10GENES[GOI] = {}
					if not LV10GENES[GOI].has_key(SP) :
						LV10COUNTS[GOI][SP] = 0
						LV10GENES[GOI][SP] = []
					if not words[0] in LV10GENES[GOI][SP] :
						LV10GENES[GOI][SP].append(words[0])
						LV10COUNTS[GOI][SP] += 1
	f.close()


w1 = open(GOLEV2COUNTTABLE, 'w')
w2 = open(GOLEV2GENESTABLE, 'w')

w1.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))

for GOID in sorted(LV2GENES.keys()) :
	w1.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	w2.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	for SP in SPECIESLIST :
		if LV2GENES[GOID].has_key(SP) :
			w1.write("\t%d" % LV2COUNTS[GOID][SP])
			w2.write("\t%s" % (" ".join(LV2GENES[GOID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()


w1 = open(GOLEV4COUNTTABLE, 'w')
w2 = open(GOLEV4GENESTABLE, 'w')

w1.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))

for GOID in sorted(LV4GENES.keys()) :
	w1.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	w2.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	for SP in SPECIESLIST :
		if LV4GENES[GOID].has_key(SP) :
			w1.write("\t%d" % LV4COUNTS[GOID][SP])
			w2.write("\t%s" % (" ".join(LV4GENES[GOID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()

w1 = open(GOLEV6COUNTTABLE, 'w')
w2 = open(GOLEV6GENESTABLE, 'w')

w1.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))

for GOID in sorted(LV6GENES.keys()) :
	w1.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	w2.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	for SP in SPECIESLIST :
		if LV6GENES[GOID].has_key(SP) :
			w1.write("\t%d" % LV6COUNTS[GOID][SP])
			w2.write("\t%s" % (" ".join(LV6GENES[GOID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()

w1 = open(GOLEV8COUNTTABLE, 'w')
w2 = open(GOLEV8GENESTABLE, 'w')

w1.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))

for GOID in sorted(LV8GENES.keys()) :
	w1.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	w2.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	for SP in SPECIESLIST :
		if LV8GENES[GOID].has_key(SP) :
			w1.write("\t%d" % LV8COUNTS[GOID][SP])
			w2.write("\t%s" % (" ".join(LV8GENES[GOID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()

w1 = open(GOLEV10COUNTTABLE, 'w')
w2 = open(GOLEV10GENESTABLE, 'w')

w1.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#GO_Term\tOntology\tDesc\t%s\n" % ("\t".join(SPECIESLIST)))

for GOID in sorted(LV10GENES.keys()) :
	w1.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	w2.write("%s\t%s\t%s"% (GOID, GOPROCESS[GOID], GODESCIPTION[GOID]))
	for SP in SPECIESLIST :
		if LV10GENES[GOID].has_key(SP) :
			w1.write("\t%d" % LV10COUNTS[GOID][SP])
			w2.write("\t%s" % (" ".join(LV10GENES[GOID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()

