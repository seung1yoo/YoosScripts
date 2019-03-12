import sys

if len(sys.argv) != 3 :
	print "Usage : python %s <FileList.txt> <OutPrefix>" % sys.argv[0]
	sys.exit()

KEGGINFO = '/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180756-JNU-Haliotis-RNAref-SmallReport-FC2P005-20181116/addTask/20181128_InterProScan/si_scripts/bin/KEGG.map.xls'
ECINFO   = '/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180756-JNU-Haliotis-RNAref-SmallReport-FC2P005-20181116/addTask/20181128_InterProScan/si_scripts/bin/EC_number.xls'

KEGGNAME = {}

f = open(KEGGINFO, 'r')

for lines in f.xreadlines() :
	if lines[0] == '#' :
		continue
	else :
		words = lines.rstrip().split("\t")
		KEGGNAME[words[0]] = words[1]

f.close()

ECNAME = {}

f = open(ECINFO, 'r')

for lines in f.xreadlines() :
	if lines[0] == '#' :
		continue
	else :
		words = lines.rstrip().split("\t")
		ECNAME[words[0]] = words[1]

f.close()

FILELIST = sys.argv[1]
COUNTTABLE = sys.argv[2] + ".count.xls"
GENESTABLE = sys.argv[2] + ".genes.xls"

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

COUNTS = {}
GENES = {}

for SP in SPECIESLIST :
	f = open(SPFILES[SP], 'r')
	for lines in f.xreadlines() :
		if lines[0] == '#' :
			continue
		else :
			words = lines.rstrip().split("\t")
			if len(words) == 13 :
				names = words[12].split("|")
				for name in names :
					PATHWAY = name.split(" ")[0]
					if PATHWAY == 'KEGG:' :
						MAPID = name.split(" ")[1].split("+")[0] + "\t" + name.split(" ")[1].split("+")[1]
						if not GENES.has_key(MAPID) :
							COUNTS[MAPID] = {}
							GENES[MAPID] = {}
						if not GENES[MAPID].has_key(SP) :
							COUNTS[MAPID][SP] = 0
							GENES[MAPID][SP] = []
						if not words[0] in GENES[MAPID][SP] :
							GENES[MAPID][SP].append(words[0])
							COUNTS[MAPID][SP] += 1
	f.close()

w1 = open(COUNTTABLE, 'w')
w2 = open(GENESTABLE, 'w')

w1.write("#KEGG_MapID\tKEGG_Name\tKEGG_GeneID\tGene_Desc\t%s\n" % ("\t".join(SPECIESLIST)))
w2.write("#KEGG_MapID\tKEGG_Name\tKEGG_GeneID\tGene_Desc\t%s\n" % ("\t".join(SPECIESLIST)))

for KEGGID in sorted(GENES.keys()) :
	KEGGID1 = "map" + KEGGID.split("\t")[0]
	GENENAME = KEGGID.split("\t")[1]
	w1.write("%s\t%s\t%s\t%s"% (KEGGID1, KEGGNAME[KEGGID1], GENENAME, ECNAME[GENENAME]))
	w2.write("%s\t%s\t%s\t%s"% (KEGGID1, KEGGNAME[KEGGID1], GENENAME, ECNAME[GENENAME]))
	for SP in SPECIESLIST :
		if GENES[KEGGID].has_key(SP) :
			w1.write("\t%d" % COUNTS[KEGGID][SP])
			w2.write("\t%s" % (" ".join(GENES[KEGGID][SP])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()

