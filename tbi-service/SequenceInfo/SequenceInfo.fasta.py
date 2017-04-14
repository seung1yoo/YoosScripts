
def run_clc_sequence_info(file, cutoff):
    outfile = "{0}.{1}.SeqInfo".format(file, cutoff)
    if not os.path.isfile(outfile):
        pass
        #cmd = "~/YoosScripts/BioTools/CLC_AC_4.4.2_linux64/"\
        #cmd = "~/YoosScripts/BioTools/CLC_AC_4.4.2_mac64/"\
        #cmd = "/BiO/BioPeople/siyoo/00.Tools/clc-assembly-cell-4.4.2-linux_64/"\
        cmd = "/TBI/People/tbi/siyoo/YoosScripts/BioTools/CLC_AC_4.4.2_linux64/"\
              "clc_sequence_info -r -n -c {2} {0} > {1}".format(file, outfile, cutoff)
        print (cmd)
        os.system(cmd)
    else:
        print ('file is exist')
    return outfile

def make_infoDict(seqinfo_file, infoDict):
    for line in open(seqinfo_file):
        line = line.strip()
        if line == '':
            continue
        items = line.split() 
        line_tag = items[0] 
        line_value = items[-1]
        if line_tag == "Number" :
            if not "%" in line_value:
                infoDict[seqinfo_file].setdefault('No_Reads', line_value)
            base_tag = items[2] 
            if "A" in base_tag:
                a = items[3]
            elif "C" in base_tag:
                c = items[3]
            elif "G" in base_tag:
                g = items[3]
            elif "T" in base_tag:
                t = items[3]
            elif "N" in base_tag:
                n = items[3]

        elif line_tag == "Total" :
            infoDict[seqinfo_file].setdefault("Residues", line_value)
            total = line_value
        elif line_tag == "Minimum" : 
            infoDict[seqinfo_file].setdefault("Minimum", line_value)
        elif line_tag == 'Maximum' : 
            infoDict[seqinfo_file].setdefault("Maximum", line_value)
        elif line_tag == 'Average' : 
            infoDict[seqinfo_file].setdefault("Average", line_value)
        elif line_tag == 'N50' : 
            infoDict[seqinfo_file].setdefault("N50", line_value)
        else :
            pass 

    gc = (int(g)+int(c))/float(total)*100
    npct = int(n)/float(total)*100
    infoDict[seqinfo_file].setdefault("Npct", npct)
    infoDict[seqinfo_file].setdefault("GC", gc)

    return infoDict
    
def makeFileList(path, extensions):
    files = []
    for extension in extensions:
        newfiles = glob.glob('{0}/*.{1}'.format(path, extension))
        files.extend(newfiles)
    return files

def number_format(n):
    s = str(n)
    if s.find('.') < 0:
        e = re.compile(r"(\d)(\d\d\d)$")
        s, cnt = re.subn(e, r"\1,\2", s)
    e = re.compile(r"(\d)(\d\d\d([.,]))")
    while 1:
        s, cnt = re.subn(e, r"\1,\2", s)
        if not cnt: break
    return s

def main(path, cutoff, extensions):
    files = makeFileList(path, extensions)

    infoDict = dict() 
    for file in files:
        seqinfo_file = run_clc_sequence_info(file, cutoff)
        infoDict.setdefault(seqinfo_file, {})
        infoDict = make_infoDict(seqinfo_file, infoDict)

    outfile = open("{0}/SequenceInfo.CutOff.{1}.Report".format(path, cutoff), "w")
    outfile.write('||{0}||{1}||{2}||{3}||{4}||{5}||{6}||{7}||{8}||\n'.format(\
            'File name', 'No. Contigs', 'Residues', 'Average',\
            'Minimum', 'Maximum', 'N50', 'Npct', 'GC'
            ))
    for info_file, ainfoDict in sorted(infoDict.items()):
        outfile.write('||{0}||{1}||{2}||{3}||{4}||{5}||{6}||{7}||{8}||\n'.format(\
                info_file.split('/')[-1],
                number_format(ainfoDict['No_Reads']),
                number_format(ainfoDict['Residues']),
                number_format(ainfoDict['Average']), 
                number_format(ainfoDict['Minimum']), 
                number_format(ainfoDict['Maximum']),
                number_format(ainfoDict['N50']),
                round(ainfoDict['Npct'], 1), round(ainfoDict['GC'], 1)
                ))

if __name__=='__main__':
    import glob
    import os
    import sys
    import re
    import argparse
    parser = argparse.ArgumentParser(description='Run clc_sequence_info and Merge')
    parser.add_argument('-p', '--path', type=str, help='Input the target PATH')
    parser.add_argument('-c', '--cutoff', type=str, default='1', help='Input the minimum cut off value')
    parser.add_argument('-e', '--extensions', nargs='+', help='Input the filename extensions as list')
    args = parser.parse_args()
    main(args.path, args.cutoff, args.extensions)
