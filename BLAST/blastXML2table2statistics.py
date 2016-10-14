import math
import sys

def makeStats(table):
    queryDic = dict()
    targetDic = dict()
    for line in open(table):
        items = line.strip('\n').split('\t')
        if line.startswith('queryNum'):
            continue
        queryDic.setdefault(items[3], 0)
        queryDic[items[3]] += 1
        targetDic.setdefault(items[8], 0)
        targetDic[items[8]] += 1
    out = open('{0}.stats'.format(table), 'w')
    out.write('query unique : {0}\n'.format(len(queryDic.keys())))
    out.write('target unique : {0}\n'.format(len(targetDic.keys())))
    out.close()

def targetOrder(table, excepts):
    targetDic = dict()
    memDic = dict()
    for line in open(table):
        items = line.strip('\n').split('\t')
        if line.startswith('queryNum'):
            continue
        query = items[3]
        targetID = items[8]
        targetDesc = items[9]

        exceptCount = 0
        for exceptChar in excepts:
            if targetDesc.startswith(exceptChar):
                exceptCount += 1
        if not exceptCount in [0]:
            continue

        targetDic.setdefault(targetID, targetDesc)
        memDic.setdefault(targetID, {}).setdefault(query, 0)
        memDic[targetID][query] += 1

    orderDic = dict()
    for targetID, queryDic in memDic.iteritems():
        orderDic.setdefault(targetID, len(queryDic.keys()))

    out = open('{0}.targetOrder'.format(table), 'w')
    for targetID, count in sorted(orderDic.iteritems(), key=lambda (k,v):(v,k), reverse=True):
        out.write('{0}\n'.format('\t'.join([targetID, targetDic[targetID], str(count), ','.join(memDic[targetID].keys())])))
    out.close()

def speciesOrder(table, excepts):
    import re
    speciesDic = dict()
    for line in open(table):
        items = line.strip('\n').split('\t')
        if line.startswith('queryNum'):
            continue
        targetDesc = items[9]

        ## species term find 1
        p = re.compile('[\[][a-zA-Z]+ [a-zA-Z]+[\]]')
        l = p.findall(targetDesc)
        if len(l) in [0]:
            species = 'Cannot search species name'
        elif len(l) in [1]:
            species = l[0]
        else:
            species = l[0]
        ##

        ##species term find 2
        #exceptCount = 0
        #for exceptChar in excepts:
        #    if targetDesc.startswith(exceptChar):
        #        exceptCount += 1
        #if not exceptCount in [0]:
        #    continue
        #species = ' '.join(targetDesc.split()[:2])
        ##

        speciesDic.setdefault(species, 0)
        speciesDic[species] += 1
    
    out = open('{0}.speciesOrder'.format(table), 'w')
    for species, count in sorted(speciesDic.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        out.write('{0}\n'.format('\t'.join([species, str(count)])))
    out.close()
    
def hspCoverageOrder(table):
    covDic = dict()
    for line in open(table):
        items = line.strip('\n').split('\t')
        if line.startswith('queryNum'):
            continue

        qCov = float(items[10])
        qCov_k = int(math.ceil(int(qCov)*0.1)*10)
        covDic.setdefault('q', {}).setdefault(qCov_k, 0)
        covDic['q'][qCov_k] += 1

        tCov = float(items[11])
        tCov_k = int(math.ceil(int(tCov)*0.1)*10)
        covDic.setdefault('t', {}).setdefault(tCov_k, 0)
        covDic['t'][tCov_k] += 1
        
    out = open('{0}.hspCoverage'.format(table), 'w')
    for tag, rangeDic in covDic.iteritems():
        for range, count in sorted(rangeDic.iteritems(), key=lambda (k,v):(v,k), reverse=True):
            out.write('{0}\n'.format('\t'.join([tag, str(range), str(count)])))
    out.close()

def bestTargetCoverage(table):
    covDic = dict()
    for line in open(table):
        items = line.strip('\n').split('\t')
        if line.startswith('queryNum'):
            continue
        q_id = items[3]
        t_id = items[9]
        t_s, t_e = [int(x) for x in sorted(items[16:18])]
        t_len = int(items[7])
        #print q_id, t_id, t_s, t_e, t_len
        for i in range(t_s, t_e+1):
            covDic.setdefault(t_id, {}).setdefault('pos', {}).setdefault(i, 0)
            covDic.setdefault(t_id, {}).setdefault('unit', {}).setdefault(q_id, 0)
            covDic.setdefault(t_id, {}).setdefault('len', t_len)

    bestDic = dict()
    for t_id, tagDic in covDic.iteritems():
        #print t_id, tagDic['len'], tagDic['unit'].keys(), tagDic['pos'].keys()
        t_cov = len(tagDic['pos'].keys())/float(tagDic['len'])*100
        bestDic.setdefault(t_id, t_cov)

    out = open('{0}.bestTargetCoverage'.format(table), 'w')
    out.write('{0}\n'.format('\t'.join(['Target_coverage', 'Target', 'covered_base(bp)', 'Target_length', 'queries', '#queries'])))
    for t_id, t_cov in sorted(bestDic.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        #print '{0}'.format('\t'.join([str(x) for x in [t_cov, t_id, ','.join(covDic[t_id]['unit'].keys()), len(covDic[t_id]['pos'].keys()), covDic[t_id]['len']]]))
        out.write('{0}\n'.format('\t'.join([str(x) for x in [t_cov, t_id, len(covDic[t_id]['pos'].keys()), covDic[t_id]['len'], ','.join(covDic[t_id]['unit'].keys()), len(covDic[t_id]['unit'].keys())]])))
    out.close()


def main(args):
    print args
    makeStats(args.table)
    targetOrder(args.table, args.excepts)
    speciesOrder(args.table, args.excepts)
    hspCoverageOrder(args.table)
    bestTargetCoverage(args.table)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--table',
            default='ContamSeqContigs_vs_Mitochondria_20140311.fna.1e-1.xml.100.1000.1e-1.table')
    parser.add_argument('-e', '--excepts',
            default=['Uncharacterized', 'Putative', 'Predicted', 'PREDICTED', 'TPA'])
    ## TPA : Third Party Annotation ; submitter-provided annotation (http://www.ncbi.nlm.nih.gov/genbank/tpa)
    args = parser.parse_args()
    main(args)

