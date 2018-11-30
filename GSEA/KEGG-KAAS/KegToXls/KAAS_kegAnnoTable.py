class KEGANNOTABLE:
	def __init__(self, args):
		self.infiles = args.kegfiles
		self.pepfile = args.pepfile
		self.outfile = args.outfile
		self.pep_geneList = self.genelist()
		self.keg_geneList, self.infoDic = self.contentDic()

	def filenameDic(self, afile):
			if 'q00001.keg' in [afile]:
				return 'Orthology'
			if 'q00002.keg' in [afile]:
				return 'Modules'
			if 'q00003.keg' in [afile]:
				return 'ReModules'

	def genelist(self):
		pep_geneList = list()
		for line in open(self.pepfile):
			if line.startswith('>'):
				gene = line.lstrip('>').split(' ')[0]
				pep_geneList.append(gene)
		return pep_geneList

	def struc_division(self, letter):
		divList = list()
		for alpha in range(ord('A'), ord(letter)):
			divList.append(chr(alpha))
		return divList

	def reverseMapDic(self, dic):
		categorys = list()
		for key in sorted(dic):
			categorys.append(dic[key])
		if 'D' not in dic.keys():
			categorys.append('-')
		return categorys

	def contentsParser(self, line, infokey):
		div = line[0]
		if div in [infokey]:
			infos = line[1:].strip()
			infos = infos.replace(';', ' ').replace('[EC', ' [EC').split('  ')
			letter = infos[1]
			if letter.startswith('M'):
				infos.append('-')
			if 'EC:' not in infos[-1]:
				infos.append('-')
			return infos

	def identiDic(self):
		for afile in self.infiles:
			divDic = dict()
			for line in open(afile):
				if line.startswith('+'):
					infoKey = line.lstrip('+').split('\t')[0]
				if line.startswith('#') or line.startswith('!'):
					continue
				if line[0] in self.struc_division(infoKey):
					div = line[0]
					contents = line[1:].strip().strip('<b>').strip('</b>')
					if contents in ['']:
						continue
					if div not in [infoKey]:
						divDic[div] = contents
				yield line, afile, divDic, infoKey

	def contentDic(self):
		infoDic = dict()
		keg_geneList = set()
		for line, afile, divDic, infoKey in self.identiDic():
			infos = self.contentsParser(line, infoKey)
			if infos in [None]:
				continue
			else:
				infos.extend(self.reverseMapDic(divDic))
				geneId = infos[0]
				items = infos[1:]
				keg_geneList.add(geneId)
				infoDic.setdefault(afile, {}).setdefault(geneId, [])
				infoDic[afile][geneId].append(items)
		return keg_geneList, infoDic

	def outfileMaker(self):
		infoDic = self.infoDic
		outFile = open(self.outfile, 'w')
		header = ['KegFileName', 'GeneID', 'Entry', 'Name',
							'Definition', 'EC:Num', 'Category_A', 'Category_B', 'Category_C', 'Category_D']
		outFile.write('{0}\n'.format('\t'.join(header)))
		for geneId in self.pep_geneList:
			if geneId not in self.keg_geneList:
					text = ['-'] * len(header)
					text[1] = geneId
					outFile.write('\t'.join(text) + '\n')
			if geneId in self.keg_geneList:
				for afile, geneDic in infoDic.iteritems():
					if geneId in geneDic.keys():
						for infos in geneDic[geneId]:
							text = '{0}\t{1}\t{2}\n'.format(self.filenameDic(afile), geneId, '\t'.join(infos))
							outFile.write(text)

def main(args):
	result = KEGANNOTABLE(args)
	result.outfileMaker()

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--kegfiles', nargs='+',
						default='q00001.keg' 'q00002.keg' 'q00003.keg')
	parser.add_argument('--pepfile',
						default='C_sinensis_plantkingdomgdb_1.pep.fa')
	parser.add_argument('--outfile',
						default='KegAnnoTable.xls')
	args = parser.parse_args()
	main(args)
