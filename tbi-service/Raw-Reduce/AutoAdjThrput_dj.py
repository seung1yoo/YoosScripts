import os, sys

adjThrScript = './FastqThrAdjustment.py'
pipeManager  = './pipelineManager.v03.py'

def getAdjThrScript(parameter) :
    
    f2 = open(sys.argv[2], 'w')

    f2.write('#fork\n')

    for File in sampLst :
        f2.write('\tpython %s %s_1.fq.gz adjThrFiles/%s_1.fastq %s\n' % (adjThrScript, File, File, float(parameter)))
        f2.write('\tpython %s %s_2.fq.gz adjThrFiles/%s_2.fastq %s\n' % (adjThrScript, File, File, float(parameter)))

    f2.write('#join')

if __name__ == "__main__" :
    if len(sys.argv) != 4 :
        print "Usage : python %s <inputFileList.txt> <tmp_adjThrFileList.txt> <parameter>" % sys.argv[0]
        sys.exit()
    
    else :
        os.system('mkdir ./adjThrFiles')
        
        f1 = open(sys.argv[1], 'r')

        sampLst = []

        for line in f1.xreadlines() :
            sampleId = line.strip()
            sampLst.append(sampleId)
        
        parameter = sys.argv[3]

        getAdjThrScript(parameter)

        os.system('%s SA < %s' % (pipeManager, sys.argv[2]))
        os.system('rm %s' % sys.argv[2])
