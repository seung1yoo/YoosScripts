import sys

if len(sys.argv) != 3 :
	print "Usage : python %s <InterProScan.out.table.xls> <Out.gene.to.GO>" % sys.argv[0]
	sys.exit()


TABLEFILE = sys.argv[1]
GENETOGOFILE = sys.argv[2] 
#GOEDGEFILE = "./go.edge.txt"
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

f = open(TABLEFILE, 'r')
w = open(GENETOGOFILE, 'w')
w.write("#GeneID\tGO_Term\tOntology\tDesc\n")

INTABLEVALUE = {}

for lines in f.xreadlines() :
	if lines[0] == '#' :
		continue
	else :
		words = lines.rstrip().split("\t")
		if len(words) < 12 :
			continue
		if words[11][:2] == 'GO' :
			name = words[11].split("|")
			for GOID in name :
				if GOPROCESS.has_key(GOID) :
					DATA = "%s\t%s\t%s\t%s\n" % (words[0], GOID, GOPROCESS[GOID], GODESCIPTION[GOID])
					if not INTABLEVALUE.has_key(DATA) :
						w.write(DATA)
						INTABLEVALUE[DATA] = 1

f.close()
w.close()

