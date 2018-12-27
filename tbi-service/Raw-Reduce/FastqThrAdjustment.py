import sys
import gzip
import commands

if len(sys.argv) != 4 :
    print "Usage : %s <oriFqFile> <adjFqFile> <parameter>" % sys.argv[0]
    sys.exit()

oriFqFile = sys.argv[1]
adjFqFile = sys.argv[2]
parameter = sys.argv[3]

totalLineNo = commands.getoutput('zcat %s | wc -l' % oriFqFile)

f1 = gzip.open(oriFqFile, 'r')
f2 = open(adjFqFile, 'w')

lineNo = 0

for line in f1.readlines() :
    
    lineNo += 1

    if int(lineNo) >= (int(totalLineNo) * float(parameter)) :
        if lineNo % 4 == 0 :
            f2.write(line)
            sys.exit()

    f2.write(line)
