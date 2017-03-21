#!/usr/bin/env python
# AUTHOR: Seung Jae Noh
# Modified by : Seung-il Yoo
# FileName: 
# Description: deg_parser_ee.py
import string, os, sys, datetime
from math import *
#import numpy as np
#import scipy.misc

#===================================================================================
#       LOG_SUMMATION
#
#       for easy calculation of log(N!)
#       cf. log(0) --> error!, summation start FROM value 1 not FROM 0

def log_sum(end):
    sum = 0.0
    k = 1
    while k <= end:
        sum += log(k, 10)
        k += 1
    return sum

#===================================================================================
#       Stirling's formula : N! ~= sqrt(2*pi*n)(n/e)**n
#                        log(N!) = 1/2*(log(2)+log(pi)+log(n)) + n*(log(n)-log(e))

def stirling(n):
    if n == 0 or n == 1:
        return 0.0
    elif n >= 2 and n < 100:
        return log_sum(n)
    elif n >= 100 and n < 10000:
        return n*(log(n, 10)-log(e, 10)) + 0.5*(log(2, 10)+log(pi, 10)+log(n, 10))
    elif n >= 10000:
        return n*(log(n, 10)-log(e, 10))

#=========================================================================================
def eePvalue_calc(x, eg, p):
    pvalue = 0.0
    eg = int(eg)
    for m in range(eg-x+1):
        num = x + m
        logP = stirling(eg) - stirling(eg-x) - stirling(x) + num*log(p, 10) + (eg-num)*log((1-p), 10)
        pvalue += pow(e, logP)
    return pvalue
#=========================================================================================

def annoParser(file):
    annotDic = {}
    for line in open(file):
        if line.strip() == '':
            continue
        items = line.strip().split('\t')
        trId = items[0]
        annotDic.setdefault(trId, items[1:])
    print ("annotDic: %d" % len(annotDic))
    return annotDic

def pListMaker(file, tissueList):
    libSumList = [0] * len(tissueList) # [s1, s2, s3,,,,s6]: si= total no. reads in each lib
    total = 0
    for line in open(file):
        if line.strip() == '':
            continue
        lineTmp = line.strip().split('\t')
        if lineTmp[0] == "TrID":
            continue
        for i in range(len(tissueList)):
            try:
                readNum = int(lineTmp[2+i])
            except:
                readNum = int(round(float(lineTmp[2+i]), 0))
            libSumList[i] += readNum
            total += readNum

    pList = [] # pi = si/sum(si)
    for i in range(len(tissueList)):
        p = float(libSumList[i]/float(total))
        pList.append(p)
        print ("%s:\t%d,\tp=%.3f" % (tissueList[i], libSumList[i], p))
    return pList

def iterator(file, tissueLen):
    for line in open(file):
        if line.strip() == '' or line.startswith('TrID'):
            continue
        lineTmp = line.strip().split('\t') 
        trId = lineTmp[0]
        trLen = lineTmp[1]
        readValues = '\t'.join(lineTmp[2:2+tissueLen])  
        fpkmValues = '\t'.join(lineTmp[2+tissueLen:2+tissueLen*2])

        ##### make readList using readValues
        readList = []
        eg = 0 # total reads of one gene from 6 tissues
        for j in readValues.split('\t'):
            try:
                readNum = int(j)
            except:
                readNum = int(round(float(j)))
            readList.append(readNum)
            eg += float(readNum)
        yield trId, trLen, readValues, fpkmValues, readList, eg

def eeMaker(trId, readList, eg, pList, tissueList):
    for k, x in enumerate(readList): # k = index, x = readList[index]
        p = pList[k]
        fi = eg*p
        try:
            ee = x/fi
        except Exception as e:
            print (trId, readList, '--> %s' % e)
            continue
        eePvalue = eePvalue_calc(x, eg, p)
        tissue = tissueList[k]
        if eePvalue >= 0.0001:
            eePval = "%.4f" % eePvalue
        else:
            eePval = "%.3e" % eePvalue
        yield k, ee, eePvalue

def dataFileMaker(filelist, tissueList):
    Trlist = {}
    for file in filelist:
        for line in file:
            sample, trId, trLen, readValue, fpkmValue = line.strip().split('\t')
            trId = '%s_%s' % (trId, trLen)
            if sample == 'Sample':
                continue
            if Trlist.has_key(trId):
                Trlist[trId]['readValue'].append(readValue)
                Trlist[trId]['fpkmValue'].append(fpkmValue)
            else:
                Trlist.update({trId:{'trLen':trLen, 'readValue':[readValue], 'fpkmValue':[fpkmValue]}})

    s = str(datetime.datetime.now()).replace(' ', '_')
    dataFile = '%s.dataFile' % s
    out = open(dataFile, 'w')
    out.write('TrID\tTrLen\t%s\t%s\n' % ('\t'.join(tissueList), '\t'.join(tissueList)))
    for trId, values in Trlist.items():
        trLen = values['trLen']
        lowExpCount = 0
        for fpkmValue in values['fpkmValue']:
            if float(fpkmValue) < 0.3:
                lowExpCount += 1
        if lowExpCount == len(values['fpkmValue']):
            continue
        outline = '%s\t%s\t%s\t%s\n' % (trId, trLen, '\t'.join(values['readValue']), '\t'.join(values['fpkmValue']))
        out.write(outline)
    out.close()
    return dataFile

def main(dataFile, rawList, annoFile, outFile, exp_cutoff, ee_cutoff, eePval_cutoff, tissueList):
    annotDic = {}
    if annoFile:
        annotDic = annoParser(annoFile)
    else:
        print ('You did not specify the annotation file; Please check it again!!!')

    if not dataFile:
        dataFile = dataFileMaker(rawList, tissueList)

    pList = pListMaker(dataFile, tissueList) # [s1, s2, s3,,,,s6]: si= total no. reads in each lib

    #======================================================================================================    
    #DEG
    degCount, count = [0] * len(tissueList), 0
    tissues = '\t'.join(tissueList[:])

    fw = open(outFile, 'w')
    fw.write("trId\ttrLen\t%s\t%s\tTissue\tEE\tPvalue_ee\tannotations\n" % (tissues, tissues))
    for trId, trLen, readValues, fpkmValues, readList, eg in iterator(dataFile, len(tissueList)):
        if trId == "TrID":
            continue
        count += 1
        if trId in annotDic:
            annos = annotDic[trId]
        else:
            annos = []
        
        eeList, eePvalList = [], []
        for k, ee, eePvalue in eeMaker(trId, readList, eg, pList, tissueList):
            eeList.append("%.2f" % ee)
            eePvalList.append(eePvalue)
            if ee >= ee_cutoff and eePvalue <= eePval_cutoff:
                newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (trId, trLen, readValues, fpkmValues, tissueList[k], ee, eePvalue, '\t'.join(annos))
                fw.write(newLine)
                degCount[k] += 1
        print (count, trId, readValues, eeList, eePvalList)
    fw.close()
    for i, tissue in enumerate(tissueList):
        print ("%s: %d" % (tissue, degCount[i]))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Parse DEG based on the expression enrichment')
    parser.add_argument('-d', '--dataFile',
                        help='Specify the database file, It should contain gene name, gene length, mapped read count, rpkm seperated by tab and title line should be start with \'TrID\'')
    parser.add_argument('-r', '--rawdatalist', nargs='+', type=argparse.FileType('r'),
                        help='Specify the raw database file which contain each library. \
                              The order of input file should be same as tissueList')
    parser.add_argument('-a', '--annoFile',
                        help='If you have the file which contain annotation, then specify that file')
    parser.add_argument('-o', '--outFile',
                        help='Specify the output file name')
    parser.add_argument('-e', '--ee_cutoff', type=float,
                        help='Specify the expression enrichment cutoff, default is 4.0',
                        default=2.0)
    parser.add_argument('-f', '--exp_cutoff', type=float,
                        help='Specify the cutoff of FPKM, default is 0.3',
                        default=0.3)
    parser.add_argument('-p', '--eePval_cutoff', type=float,
                        help='Specify the pvalue cutoff of expression enrichment, default is 1e-3',
                        default=1e-3)
    parser.add_argument('-t', '--tissueList', nargs='+', 
                        help='Specify the library names, the name order should be same as the one in the dataFile',
                        default=[])
    args = parser.parse_args()
    main(args.dataFile, args.rawdatalist, args.annoFile, args.outFile, args.exp_cutoff, args.ee_cutoff, args.eePval_cutoff, args.tissueList)
    
