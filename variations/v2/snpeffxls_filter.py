

def positions_from_vcf(vcf_fn):
    positions = list()
    for line in open(vcf_fn):
        if line.startswith('#'):
            continue
        items = line.rstrip('\n').split('\t')
        chrom = items[0]
        pos = items[1]
        positions.append('{0}_{1}'.format(chrom, pos))
    return positions

def select_only_positions_from_snpefftxt(positions, inxls, outxls):
    out = open(outxls,'w')
    for line in open(inxls):
        if line.startswith('#'):
            out.write(line)
            continue
        items = line.rstrip('\n').split('\t')
        chrom = items[0]
        pos = items[1]
        #
        if '{0}_{1}'.format(chrom, pos) in positions:
            out.write(line)
    out.close()




def main():
    positions = positions_from_vcf('L_edodes_newgeneset.SNP.bi.good.30-70.vcf.VarComparing')
    select_only_positions_from_snpefftxt(positions, 'multisample.snpeff.xls', 'L_edodes_newgeneset.SNP.bi.good.30-70.vcf.VarComparing.xls')
    #
    positions = positions_from_vcf('L_edodes_newgeneset.SNP.bi.good.20-80.vcf.VarComparing')
    select_only_positions_from_snpefftxt(positions, 'multisample.snpeff.xls', 'L_edodes_newgeneset.SNP.bi.good.20-80.vcf.VarComparing.xls')
    #
    positions = positions_from_vcf('L_edodes_newgeneset.SNP.bi.vcf.VarComparing')
    select_only_positions_from_snpefftxt(positions, 'multisample.snpeff.xls', 'L_edodes_newgeneset.SNP.bi.vcf.VarComparing.xls')


if __name__=='__main__':
    main()
