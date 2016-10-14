
def main(blastxml, outtable, parameters):
    out = open(outtable, 'w')
    titles = ['queryNum','hit_num','hsp_num','queryID','e-value','bit-score','queryLen','targetLen',
              'targetAcc','targetDesc','queryCov','targetCov','alignLen','alignSim','queryStr','queryEnd',
              'targetStr','targetEnd','gapInfo','frameInfo','querySeq','comparison','targetSeq','targetDB']
    print '{0}'.format('\t'.join(titles))
    out.write('{0}\n'.format('\t'.join(titles)))
    parameters = [float(x) for x in parameters]
    hit_num_cut, hsp_num_cut, evalue_cut = parameters
    query_num = 0
    for blast_record in NCBIXML.parse(open(blastxml)):
        query_num += 1
        alignment_num = 0
        for alignment in blast_record.alignments:
            alignment_num += 1
            hsp_num = 0
            for hsp in alignment.hsps:           
                hsp_num += 1
                mismatch = hsp.match.count('+')+hsp.match.count(' ')
                alignRate = (hsp.align_length-mismatch)/float(hsp.align_length)*100
                query_mapRegions = [hsp.query_start, hsp.query_end]
                query_mapRegions.sort()
                queryRate = (query_mapRegions[1]-query_mapRegions[0]+1)/float(blast_record.query_length)*100
                sbj_mapRegions = [hsp.sbjct_start, hsp.sbjct_end]
                sbj_mapRegions.sort()
                sbjctRate = (sbj_mapRegions[1]-sbj_mapRegions[0]+1)/float(alignment.length)*100
                '''
                try:
                    subjectframe, queryframe = hsp.frame
                except:
                    subjectframe, queryframe = str(hsp.frame).strip(')').strip('(').split(',')
                    if not subjectframe:
                        subjectframe = 0
                    if not queryframe:
                        queryframe = 0

                if int(queryframe) < 0:
                    queryRate = (hsp.query_start-hsp.query_end+1)/float(blast_record.query_length)*100
                else:
                    queryRate = (hsp.query_end-hsp.query_start+1)/float(blast_record.query_length)*100
                if int(subjectframe) < 0:
                    sbjctRate = (hsp.sbjct_start-hsp.sbjct_end+1)/float(alignment.length)*100
                else:
                    sbjctRate = (hsp.sbjct_end-hsp.sbjct_start+1)/float(alignment.length)*100
                '''
                lineUnits = [query_num, alignment_num, hsp_num,
                             blast_record.query, hsp.expect, hsp.bits, 
                             blast_record.query_length, alignment.length,
                             alignment.accession, alignment.hit_def, queryRate, sbjctRate,
                             hsp.align_length, alignRate, hsp.query_start, hsp.query_end,
                             hsp.sbjct_start, hsp.sbjct_end, hsp.gaps, hsp.frame,
                             hsp.query, hsp.match, hsp.sbjct, alignment.hit_id]
                lineUnits = [str(x) for x in lineUnits]
                if alignment_num <= hit_num_cut and hsp_num <= hsp_num_cut and hsp.expect <= evalue_cut:
                    #print '{0}'.format('\t'.join(lineUnits))
                    out.write('{0}\n'.format('\t'.join(lineUnits)))
                else:
                    print '{0}'.format('\t'.join(lineUnits))

if __name__=='__main__':
    from Bio.Blast import NCBIXML
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--blastxml',
            default='SoftCoral_vs_nr.blastx.xml')
    parser.add_argument('-o', '--outtable',
            default='SoftCoral_vs_nr.blastx.table')
    parser.add_argument('-p', '--parameters', nargs='+', help='[hit_num, hsp_num, evalue]',
            default=[10, 1, 1e-5])
    args = parser.parse_args()
    main(args.blastxml, args.outtable, args.parameters)

