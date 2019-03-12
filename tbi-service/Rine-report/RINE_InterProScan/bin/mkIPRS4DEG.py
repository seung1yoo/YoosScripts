#!/usr/bin/python

import sys
import glob

def mkGenLst(file) :
	fr = open(file, "r")

	degName = file.split(".")[1]
	degLst  = list()

	for line in fr.xreadlines() :
		if line.startswith("GeneAcc") :
			continue

		else :
			splitted = line.strip().split("\t")
			genId    = splitted[0]

			degLst.append(genId)

	fr.close()
	return degName, degLst

def main() :
	files = glob.glob("GenLst*txt")

	for file in files :

		(degName, degLst) = mkGenLst(file)

		fr = open(sys.argv[1], "r")
		fw = open("InterproScanAnnotation.{0}.xls".format(degName), "w")

		print file

		for line in fr.xreadlines() :
			if line.startswith("#Query") :
				fw.write(line)

			else :
				splitted = line.strip().split("\t")
				genId    = splitted[0].split("_")[0]
				others   = splitted[1:]

				if genId in degLst :
					fw.write("{0}\t{1}\n".format(genId, "\t".join(others)))

		fr.close()
		fw.close()

if __name__ == "__main__" :
	if len(sys.argv) != 2 :
		print "Usage : {0} <InterproScanAnnotation.xls>".format(sys.argv[0])
		sys.exit()

	else :
		main()
