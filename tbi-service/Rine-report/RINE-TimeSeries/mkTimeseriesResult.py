import sys

if len(sys.argv) != 6 :
    print "Usage : %s <genes.xls> <colname> <prefix> <log2FC_cut> <p_cut>" % sys.argv[0]
    sys.exit()

f2 = open(sys.argv[2], 'r')
PREFIX = str(sys.argv[3])
FCFILT   = float(sys.argv[4]) # log
p_cut = float(sys.argv[5]) # p-value

COLNAMES = []
for line in f2.xreadlines() :
    if line.rstrip() == '' or line[0] == '#' :
        continue
    COLNAMES.append(line.rstrip())


f1 = open(sys.argv[1], 'r')

COLNUM   = {}
TIMELIST = {}
header   = ''

for line in f1.xreadlines() :
    if line.startswith('Order') :
        words = line.rstrip().split("\t")
        header = line.rstrip()
        for num in range(0, len(words)) :
            if words[num] in COLNAMES :
                COLNUM[words[num]] = num

    else :
        pattern = ''
        CNT     = 0
        words = line.rstrip().split("\t")
        for NAME in COLNAMES :
            #
            p_value = words[COLNUM[NAME]+1]
            if p_value in ['-']: # p-value
                p_value = 1.0
            else:
                p_value = float(p_value)

            if p_value > p_cut:
                CNT += 1
            #
            log2fc_value = 0.0
            if words[COLNUM[NAME]+3] in ['-']: # Log2FC value
                log2fc_value = 0.0
            else:
                log2fc_value = float(words[COLNUM[NAME]+3])

            if float(log2fc_value) in [0.0]:
                CNT += 1
            #
            #####################
            '''
            if words[COLNUM[NAME]] == 'Y' :
                if float(log2fc_value) <= float(-(FCFILT)) :
                    pattern += 'D'
                elif float(log2fc_value) >= float(FCFILT) :
                    pattern += 'U'
            else :
                pattern += 'F'
            '''
            if -(FCFILT) < float(log2fc_value) < FCFILT:
                pattern += 'F'
            elif float(log2fc_value) <= -(FCFILT) :
                pattern += 'D'
            elif float(log2fc_value) >= FCFILT :
                pattern += 'U'
            else:
                print 'error: check the fold change'
                import sys
                sys.exit()
            #'''
            #####################

        if not CNT == 0 :
            #pass
            continue

        if not TIMELIST.has_key(pattern) :
            TIMELIST[pattern] = []
        TIMELIST[pattern].append("%s" % "\t".join(words))

for PAT in TIMELIST.keys() :
    w = open("pre_%s_%s.xls" % (PREFIX, PAT), 'w')
    w.write("%s\n" % header)
    w.write("%s" % "\n".join(TIMELIST[PAT]))
    w.close()

