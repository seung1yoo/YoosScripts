import glob
import subprocess
import os

def fileFinder(path, extension, firstTag, secondTag):
    fileDic = dict()
    #
    cmd = "find {0} -iname '*.{1}'".format(path, extension)
    filesBites = subprocess.check_output(cmd, shell=True)
    for afile in filesBites.decode("utf-8").split('\n'):
        if not afile:
            continue
        fileName = afile.split('/')[-1].rstrip(extension)
        if fileName.endswith(firstTag):
            sampleName = fileName.strip(firstTag)
            fileDic.setdefault(sampleName, {}).setdefault('R1', os.path.abspath(afile))
        elif not secondTag:
            continue
        elif fileName.endswith(secondTag):
            sampleName = fileName.strip(secondTag)
            fileDic.setdefault(sampleName, {}).setdefault('R2', os.path.abspath(afile))
        else:
            print("ERROR: Check the Tag arguments")
            import sys
            sys.exit()
    #
    return fileDic

def exe_FasterFastqStatistics(script, rs):
    rstFile = '{0}.rst.xls'.format(rs[0])
    if os.path.isfile(rstFile):
        pass
    else:
        cmd = '{0} {1}'.format(script, ' '.join(rs))
        print(cmd)
        os.system(cmd)
    return rstFile

def make_summary(rstDic, outFile):
    out = open(outFile, 'w')
    titles = ['#sample']
    for sampleName, rstFile in rstDic.items():
        for line in open(rstFile):
            if line.startswith('#totalReadCnt'):
                if len(titles) in [1]:
                    titles.extend(line.rstrip('\n').split('\t'))
                    out.write('{0}\n'.format('\t'.join(titles)))
                continue
            items = line.rstrip('\n').split('\t')
            out.write('{0}\t{1}\n'.format(sampleName, '\t'.join(items)))
        #
    out.close()

def main(args):
    print(args)
    fileDic = fileFinder(args.path, args.extension, args.firstTag, args.secondTag)
    #
    rstDic = dict()
    for sampleName, tagDic in fileDic.items():
        print (sampleName, tagDic)
        if 'R1' in tagDic and 'R2' in tagDic:
            r1 = tagDic['R1']
            r2 = tagDic['R2']
            rs = [r1, r2]
            rstFile =exe_FasterFastqStatistics(args.script, rs)
            rstDic.setdefault(sampleName, rstFile)
        elif 'R1' in tagDic and not 'R2' in tagDic:
            r1 = tagDic['R1']
            rs = [r1]
            rstFile = exe_FasterFastqStatistics(args.script, rs)
            rstDic.setdefault(sampleName, rstFile)
        else:
            print("ERROR: Check the Tag arguments")
            import sys
            sys.exit()
        #
    #
    make_summary(rstDic, args.outFile)



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--script', help='FasterFastqStatistics',
            default='~/YoosScripts/tbi-service/SequenceInfo/FasterFastqStatistics')
    parser.add_argument('-p', '--path', help='The PATH for file search',
            default='/BiO/BioProjects/TBD170409-SCHU-Fungi-smallRNA-20170818/Rawdata')
    parser.add_argument('-e', '--extension', help='File extension to search', default='fq.gz')
    parser.add_argument('-1', '--firstTag', help='first read tag keyword', default='1')
    parser.add_argument('-2', '--secondTag', help='second read tag keyword', default='')
    parser.add_argument('-o', '--outFile', help='Summary file', default='SequenceInfo.fastq.Report.xls')
    args = parser.parse_args()
    main(args)
