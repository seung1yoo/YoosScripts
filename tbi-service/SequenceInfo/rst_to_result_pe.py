
import json

def main(args):
	sample_dic = json.load(open(args.conf_json))

	f2 = open(args.outfn, 'w')
	f2.write('SampleID\tTBI_ID\tTotalReads\tTotalBases\tTotalBases(Gb)\tGC_Count\tGC_Rate\tN_ZeroReads\tN_ZeroReadsRate\tN5_LessReads\tN5_LessReadsRate\tN_Count\tN_Rate\tQ30_MoreBases\tQ30_MoreBasesRate\tQ20_MoreBases\tQ20_MoreBasesRate\n')

	TotalReadCnt = 0
	TotalLength  = 0
	TotalGCCnt   = 0
	NZeroReadCnt = 0
	N5ReadCnt    = 0
	TotalNCnt    = 0
	TotalQ30     = 0
	TotalQ20     = 0

	for statFile, infoDic in sample_dic.items():
		outline = []
		for stats in open(statFile):
			if stats.startswith('#') :
				continue
			else :
				words        = stats.strip('\n').split("\t")
				TotalReadCnt = int(words[0])
				TotalLength  = int(words[1])
				TotalGCCnt   = int(words[2])
				NZeroReadCnt = int(words[3])
				N5ReadCnt    = int(words[4])
				TotalNCnt    = int(words[5])
				TotalQ30     = (int(words[8])  + int(words[9]))
				TotalQ20     = (int(words[12]) + int(words[13]))
				TotalLengthGB = "%.2f Gb" % (TotalLength * 1.0 / 1000000000)
				GCRate        = TotalGCCnt   * 100.0 / TotalLength
				NzRate        = NZeroReadCnt * 100.0 / TotalReadCnt
				N5Rate        = N5ReadCnt    * 100.0 / TotalReadCnt
				TotalNRate    = TotalNCnt    * 100.0 / TotalLength
				TotalQ30Rate  = TotalQ30     * 100.0 / TotalLength
				TotalQ20Rate  = TotalQ20     * 100.0 / TotalLength

				outline.append(infoDic['cst_id'])
				outline.append(infoDic['tbi_id'])
				outline.append(str(TotalReadCnt))
				outline.append(str(TotalLength))
				outline.append(str(TotalLengthGB))
				outline.append(str(TotalGCCnt))
				outline.append("%.2f%%" % GCRate)
				outline.append(str(NZeroReadCnt))
				outline.append("%.2f%%" % NzRate)
				outline.append(str(N5ReadCnt))
				outline.append("%.2f%%" % N5Rate)
				outline.append(str(TotalNCnt))
				outline.append("%.2f%%" % TotalNRate)
				outline.append(str(TotalQ30))
				outline.append("%.2f%%" % TotalQ30Rate)
				outline.append(str(TotalQ20))
				outline.append("%.2f%%" % TotalQ20Rate)

		f2.write('\t'.join(outline)+'\n')
		outline = []

if __name__ == "__main__" :
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--conf-json')
	parser.add_argument('--outfn')
	args = parser.parse_args()
	main(args)
