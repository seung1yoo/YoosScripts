#!/usr/bin/python

import os
import sys
import gzip

class Foo :
    pass

class FastQ :
    def __init__(self,fn) :
        self._fn=fn
    def __iter__(self) :
        _in=None
        if self._fn.endswith(".gz") :
            _in=gzip.open(self._fn)
        else :
            _in=file(self._fn)
        while True :
            f=Foo()
            f.title=_in.readline().strip()[1:]
            f.seq=_in.readline().strip()
            f.comment=_in.readline().strip()
            f.score=_in.readline().strip()
            if f.title=='' or f.title==None :
                break
            yield f
        _in.close()

def reverse(seq) :
    buf=[]
    for nt in seq.upper() :
        if nt=='A' :
            buf.append('T')
        elif nt=='T' :
            buf.append('A')
        elif nt=='G' :
            buf.append('C')
        elif nt=='C' :
            buf.append('G')
        elif nt=='N' :
            buf.append('N')
        elif nt=='X' :
            buf.append('X')
        else :
            sys.stderr.write("ERROR: Unknown NT (%s)\n"%nt)
    buf.reverse()
    return "".join(buf)

def usage(out) :
    out.write("python %s [forward|reverse] [IN_FQ] [OUT_FA]\n"%__file__)

def main() :
    if len(sys.argv)<3 :
        usage(sys.stderr)
        return

    is_forward=True
    if sys.argv[1]=='forward' :
        is_forward=True
    elif sys.argv[1]=='reverse' :
        is_forward=False
    else :
        usage(sys.stderr)
        return

    in_fn=sys.argv[2]
    out_fn=sys.argv[3]

    out=file(out_fn,'wt')
    for read in FastQ(in_fn) :
        out.write(">%s\n"%read.title)
        if is_forward :
            out.write("%s\n"%read.seq)
        else :
            out.write("%s\n"%reverse(read.seq))
    out.close()


if __name__=='__main__' :
    main()

