#!/usr/bin/python

import os
import sys
import copy
import unittest

#GFF_READ="/BiO/BioTools/Rine/Tools/Cufflinks/current/gffread"
GFF_READ="gffread"

TYPE_GENE_CODING=10
TYPE_GENE_NONCODING=20

TYPE_MOL_UNKN=0
TYPE_MOL_GENE=10
TYPE_MOL_TRAN=20
TYPE_MOL_EXON=30
TYPE_MOL_CDS=40

SELECT_S=("X","Y","W","Z","MT","Mt","M","PT","Pt","P")

class Model :
    def __init__(self):
        self.gene_name=None
        self._type_mol=TYPE_MOL_UNKN
        self._coding_type=None
        self.tran=None
        self.gene=None
        self._score=None
        self._orf=None
    def id(self):
        return self._id
    def start(self):
        return self._start
    def end(self):
        return self._end
    def strand(self):
        return self._strand
    def __len__(self):
        return abs(self.end()-self.start()+1)
    def get_length(self):
        _sum=0
        if self.exon_s!=None and len(self.exon_s)>0 :
            for exon in self.exon_s :
                _sum+=len(exon)
            return _sum
        if self.cds_s!=None and len(self.cds_s)>0 :
            for exon in self.cds_s :
                _sum+=len(exon)
            return _sum
    def to_gtf(self,to_exon=True,exon_number=None):
        buf=[]
        buf.append(self._chr.name)
        buf.append(self.source)
        if self._type_mol==TYPE_MOL_GENE:
            buf.append('gene')
        elif self._type_mol==TYPE_MOL_TRAN:
            buf.append('transcript')
        elif to_exon :
            buf.append('exon')
        else :
            buf.append(self.feature)
        buf.append("%s"%(self.start()))
        buf.append("%s"%(self.end()))
        if self._score==None :
            buf.append(".")
        else :
            buf.append(self._score)
        buf.append(self._strand)
        if self._orf==None :
            buf.append(".")
        else :
            buf.append(self._orf)

        sub=[]
        #
        sub.append('gene_id "%s";'%(self.gene_id))
        if self.tran!=None and self.tran.gene._coding_type==TYPE_GENE_CODING :
            sub.append('gene_biotype "protein_coding";')
        elif self.gene!=None and self.gene._coding_type==TYPE_GENE_CODING :
            sub.append('gene_biotype "protein_coding";')
        elif self._type_mol==TYPE_MOL_GENE and self._coding_type==TYPE_GENE_CODING:
            sub.append('gene_biotype "protein_coding";')
        else :
            sub.append('gene_biotype "misc_RNA";')
        if self.gene_name!=None :
            sub.append('gene_name "%s";'%(self.gene_name))
        #
        if self._type_mol>TYPE_MOL_GENE :
            sub.append('transcript_id "%s";'%(self.tran_id))
            if self.tran!=None and self.tran._coding_type==TYPE_GENE_CODING :
                sub.append('transcript_biotype "protein_coding";')
            elif self._coding_type==TYPE_GENE_CODING :
                sub.append('transcript_biotype "protein_coding";')
            else :
                sub.append('transcript_biotype "misc_RNA";')
        #
        if self._type_mol>TYPE_MOL_TRAN :
            if exon_number!=None :
                sub.append('exon_number "%d";'%exon_number)
        buf.append(" ".join(sub))
        return "\t".join(buf)
    def __eq__(self,other) :
        if other==None :
            return False
        if self._strand!=other.strand :
            return False
        if self.start()!=other.start() :
            return False
        if self.end()!=other.end() :
            return False
        if self.id()!=None and self.id()!=other.id() :
            return False
        return True
    def __cmp__(self,other) :
        if other==None :
            return -1

        if self._chr<other._chr :
            return -1
        if self._chr>other._chr :
            return 1

        if self._type_mol in (TYPE_MOL_GENE,TYPE_MOL_TRAN) :
            if self.start()<other.start() :
                return -1
            if self.start()>other.start() :
                return 1
            if self.end()<other.end() :
                return -1
            if self.end()>other.end() :
                return 1
            if self.strand()=='+' and other.strand()=='-' :
                return -1
            if self.strand()=='-' and other.strand()=='+' :
                return 1
            return 0

        if self.strand()=='-' :
            if self.start()<other.start() :
                return 1
            if self.start()>other.start() :
                return -1
            if self.end()<other.end() :
                return 1
            if self.end()>other.end() :
                return -1
            return 0

        if self.start()<other.start() :
            return -1
        if self.start()>other.start() :
            return 1
        if self.end()<other.end() :
            return -1
        if self.end()>other.end() :
            return 1
        return 0



class DNA :
    def __init__(self,name):
        self.name=name
        self._int_name=None
        self.type=None
        if name.lower().startswith("chr") :
            self.int_name()
    def int_name(self):
        if self._int_name==None :
            try :
                self._int_name=int(self.name)
            except ValueError:
                pass
        return self._int_name
    def roman_name(self):
        return ConvertRomanNumeralToNumber(self.name)
    def is_ref(self):
        if self.type=='REF' :
            return True
        return False
    def is_chromosome(self):
        if self.name.lower().startswith("chr") :
            return True
        if self.name.lower().startswith("scaffold") :
            return False
        if self.name.lower().startswith("contig") :
            return False
        return True
    def has_int_name(self):
        if self._int_name!=None :
            return True
        try :
            self._int_name=int(self.name)
        except ValueError:
            return False
        return True
    def has_roman_name(self):
        self._roman_name=ConvertRomanNumeralToNumber(self.name)
        if self._roman_name==None :
            return False
        return True
    def has_sel_name(self):
        if self.name in SELECT_S :
            return True
        return False
    def __eq__(self,other):
        if other==None :
            return False
        if self.name!=other.name :
            return False
        return True
    def __cmp__(self,other) :
        if not(self.is_chromosome() and  other.is_chromosome()) :
            if self.is_chromosome() :
                return -1
            if other.is_chromosome() :
                return 1

        if self.has_int_name() and other.has_int_name() :
            if self.int_name()<other.int_name() :
                return -1
            if self.int_name()>other.int_name() :
                return 1
            return 0
        #
        if self.has_int_name() :
            return -1
        if other.has_int_name() :
            return 1
        #
        if self.has_roman_name() and other.has_roman_name() :
            if self.roman_name()<other.roman_name() :
                return -1
            if self.roman_name()>other.roman_name() :
                return 1
            return 0
        #
        if self.has_sel_name() and other.has_sel_name() :
            idx_1=SELECT_S.index(self.name)
            idx_2=SELECT_S.index(other.name)
            if idx_1<idx_2 :
                return -1
            if idx_1>idx_2 :
                return 1
            return 0
        #
        if self.has_sel_name() :
            return -1
        if other.has_sel_name() :
            return 1
        #
        if self.name<other.name :
            return -1
        if self.name>other.name :
            return 1
        return 0

def ConvertRomanNumeralToNumber(roman_numeral):
    '''
    http://pythonicprose.blogspot.kr/2011/02/python-roman-numeral-to-number.html
    '''
    is_init=False
    number_result = 0
    roman_numerals = {
        1:"I", 4:"IV", 5:"V", 9:"IX",
        10:"X", 40:"XL", 50:"L",
        90:"XC", 100:"C", 400:"CD",
        500:"D", 900:"CM", 1000:"M"}

    for numeral_value in sorted(roman_numerals,
            key=lambda roman: len(roman_numerals[roman]),
            reverse=True):
        keep_converting = True
        while keep_converting:
            if roman_numeral.find(roman_numerals[numeral_value]) != -1:
                number_result += numeral_value
                is_init=True
                roman_numeral = roman_numeral.replace(\
                        roman_numerals[numeral_value], "", 1)
            else:
                keep_converting = False
    if not is_init :
        return
    if len(roman_numeral)>0 :
        return
    return number_result

def parse_attr_s(line) :
    attr_s={}
    buf=[]
    is_quot=False
    for ch in line :
        if ch=='"' :
            is_quot=not is_quot
        if ch==';' and (not is_quot) :
            attr=("".join(buf)).strip()
            buf=[]
            yield attr
        else :
            buf.append(ch)

def read_gtf_line(line,feature=None) :
    unit_s=line.strip().split("\t",8)
    if feature!=None and unit_s[2]!=feature :
        return

    foo=Model()
    foo._chr=DNA(unit_s[0])
    foo.source=unit_s[1]
    foo.feature=unit_s[2]
    foo._start=int(unit_s[3])
    foo._end=int(unit_s[4])
    foo._score=unit_s[5]
    foo._strand=unit_s[6]
    foo._orf=unit_s[7]

    foo.gene_id=None
    foo.gene_name=None

    foo.tran_id=None
    foo.exon_number=None

    if foo.feature.lower()=='exon'  :
        foo._type_mol=TYPE_MOL_EXON
    elif foo.feature.lower()=='cds'  :
        foo._type_mol=TYPE_MOL_CDS

    for attr in parse_attr_s(unit_s[8]) :
        attr=attr.strip()
        if attr=='' :
            continue
        key,val=attr.split(None,1)
        if key=='gene_id' :
            foo.gene_id=val[1:-1]
        elif key=='gene_name' :
            foo.gene_name=val[1:-1]
        elif key=='transcript_id' :
            foo.tran_id=val[1:-1]
        elif key=='transcript_accession' :
            foo._name=val[1:-1]
        elif key=='exon_number' :
            foo.exon_number=int(val[1:-1])

    if foo.gene_id==None :
        foo.gene_id="%sG"%foo.tran_id

    return foo

def _update_transcript_position(tran) :
    if len(tran.exon_s)>0 :
        tran._start=tran.exon_s[0]._start
        tran._end=tran.exon_s[0]._end
        tran._length=0
        for exon in tran.exon_s :
            if tran._start > exon._start :
                tran._start=exon._start
            if tran._end < exon._end :
                tran._end=exon._end
            tran._length+=len(exon)
    elif len(tran.cds_s)>0 :
        tran._start=tran.cds_s[0]._start
        tran._end=tran.cds_s[0]._end
        tran._length=0
        for cds in tran.cds_s :
            if tran._start > cds._start :
                tran._start=cds._start
            if tran._end < cds._end :
                tran._end=cds._end
            tran._length+=len(cds)
    else :
        print "ERROR: No EXON or CDS infor."

def _read_gtf_transcript(_in) :
    tran=None
    while True :
        line=_in.readline()
        if line=='' :
            break
        if line.startswith("#") :
            continue
        #
        exon=read_gtf_line(line)
        if exon==None :
            continue
        #
        if tran!=None and tran.tran_id==exon.tran_id  :
            if exon.feature=='exon' :
                tran.exon_s.append(exon)
            elif exon.feature=='CDS' :
                tran.cds_s.append(exon)
            exon.tran=tran
            continue

        if tran!=None :
            _update_transcript_position(tran)
            _update_transcript_coding_type(tran)
            yield tran
        #
        tran=copy.copy(exon)
        tran._id=tran.tran_id
        tran.exon_s=[]
        tran.cds_s=[]
        if exon.feature=='exon' :
            tran.exon_s=[exon]
        elif exon.feature=='CDS' :
            tran.cds_s=[exon]
        else :
            print "ERROR: Unknown FEATURE: %s"%(tran.feature)
        tran._type_mol=TYPE_MOL_TRAN
        exon.tran=tran

    if tran!=None :
        _update_transcript_position(tran)
        _update_transcript_coding_type(tran)
        yield tran

def read_gtf_transcript(fn) :
    return _read_gtf_transcript(file(fn))

def _update_gene_position(gene) :
    min_start=None
    max_end=None
    for tran in gene.tran_s :
        if min_start==None :
            min_start=tran.start()
            max_end=tran.end()
        else :
            min_start=min(min_start,tran.start())
            max_end=max(max_end,tran.end())
    gene._start=min_start
    gene._end=max_end

    gene.tran_s.sort()

def _update_transcript_coding_type(tran) :
    if tran.source==None :
        tran._coding_type=TYPE_GENE_CODING
        return
    if tran.source=='protein_coding' :
        tran._coding_type=TYPE_GENE_CODING
    elif tran.source=='tRNA' or tran.source=='MT_tRNA' :
        tran._coding_type=TYPE_GENE_NONCODING_TRNA
    elif tran.source=='rRNA' or tran.source=='MT_rRNA' :
        tran._coding_type=TYPE_GENE_NONCODING_RRNA
    elif tran.source.startswith("IG_") or  tran.source.startswith("TR_")  :
        has_coding=True
        tran._coding_type=TYPE_GENE_CODING
    else :
        tran._coding_type=TYPE_GENE_CODING


def _update_gene_coding_type(gene) :
    has_coding=False
    has_linc=False
    has_nc=False
    has_rrna=False
    has_trna=False

    for tran in gene.tran_s :
        if tran.source==None :
            continue
        if tran.source=='protein_coding' :
            has_coding=True
            tran._coding_type=TYPE_GENE_CODING
            break
        elif tran.source=='lincRNA' :
            has_linc=True
            tran._coding_type=TYPE_GENE_NONCODING_LINC
        elif tran.source=='tRNA' :
            has_trna=True
            tran._coding_type=TYPE_GENE_NONCODING_TRNA
        elif tran.source=='rRNA' :
            has_rrna=True
            tran._coding_type=TYPE_GENE_NONCODING_RRNA
        else :
            has_coding=True
            tran._coding_type=TYPE_GENE_CODING

    if has_coding :
        gene._coding_type=TYPE_GENE_CODING
    elif has_linc :
        gene._coding_type=TYPE_GENE_NONCODING_LINC
    elif has_rrna :
        gene._coding_type=TYPE_GENE_NONCODING_RRNA
    elif has_trna :
        gene._coding_type=TYPE_GENE_NONCODING_TRNA
    elif has_nc :
        gene._coding_type=TYPE_GENE_NONCODING
    else :
        print "UNKNOWN TYPE: %s"%(gene.id())


def read_gtf_gene(fn) :
    _in=file(fn)
    gene=None
    for tran in _read_gtf_transcript(_in) :
        if gene!=None and tran.gene_id==gene.gene_id :
            gene.tran_s.append(tran)
            tran.gene=gene
            continue
        #
        if gene!=None :
            _update_gene_position(gene)
            _update_gene_coding_type(gene)
            yield gene
        #
        gene=Model()
        gene._id=tran.gene_id
        gene.gene_id=tran.gene_id
        gene._name=tran.gene_name
        gene.source=tran.source
        gene._type_mol=TYPE_MOL_GENE
        #
        gene._chr=tran._chr
        gene._strand=tran._strand
        #
        gene._coding_type=TYPE_GENE_CODING
        gene.tran_s=[tran]
        #
        tran.gene=gene

    if gene!=None :
        _update_gene_position(gene)
        _update_gene_coding_type(gene)
        yield gene


def usage() :
    print "python %s [GFF_FN] [GTF_FN]"%(os.path.split(__file__)[1])


def main() :
    if len(sys.argv)<3 :
        usage()
        return

    gff_fn=sys.argv[1]
    gtf_fn=sys.argv[2]

    tmp_gtf_fn=gff_fn

    if gff_fn.endswith(".gff") or gff_fn.endswith(".gff3") :
        tmp_gtf_fn="%s.raw.gtf"%(gff_fn[:-5])
        print "%s >> %s >> %s"%(gff_fn,tmp_gtf_fn,gtf_fn)
        if not os.path.exists(tmp_gtf_fn) :
            os.system("%s %s -o %s -T"%(GFF_READ,gff_fn,tmp_gtf_fn))
        else :
            print "PASS: FOUND RAW_GTF_FN %s"%(tmp_gtf_fn)
        if (not os.path.exists(tmp_gtf_fn)) or os.path.getsize(tmp_gtf_fn)<1500 :
            print "ERROR: INCOMPLETE GTF: %s"%(tmp_gtf_fn)
            print "%s %s -o %s -T"%(GFF_READ,gff_fn,tmp_gtf_fn)
            return
    else :
        print "%s >> %s"%(gff_fn,gtf_fn)

    print "STAGE\tstarting\tto\tread GTF_FN(%s)"%tmp_gtf_fn
    gene_id_s=set()
    gene_s={}
    for gene in read_gtf_gene(tmp_gtf_fn) :
        if gene.id() in gene_id_s :
            print "WARNING\tREDUNDANT GENE_ID_S: %s"%(gene.id())
            gene_org=gene_s[gene.id()]
            for tran in gene.tran_s :
                has_org=False
                for tran_org in gene_org.tran_s :
                    if tran.id()==tran_org.id() :
                        for exon in tran.exon_s :
                            exon.tran=tran_org
                            tran_org.exon_s.append(exon)
                        for cds in tran.cds_s :
                            cds.tran=tran_org
                            tran_org.cds_s.append(cds)
                        has_org=True
                        print "\tMERGED ISOFORM: %s: %d/ %d"%\
                                (tran.id(),len(tran.exon_s),len(tran_org.exon_s))
                        _update_transcript_position(tran_org)
                        break
                if not has_org :
                    print "\tADD NOVEL ISOFORM: %s"%(tran.id())
                    gene_org.tran_s.append(tran)
                    tran.gene=gene_org
                    _update_gene_position(gene_org)
                    _update_gene_coding_type(gene_org)
        else :
            gene_id_s.add(gene.id())
            gene_s[gene.id()]=gene

    gene_s=gene_s.values()
    gene_s.sort()

    for gene in gene_s :
        gene.tran_s.sort()
        for tran in gene.tran_s :
            tran.exon_s.sort()
            tran.cds_s.sort()

    print "Found %s Genes"%len(gene_s)
    out=file(gtf_fn,'wt')
    for gene in gene_s :
        out.write("%s\n"%gene.to_gtf())
        for tran in gene.tran_s :
            out.write("%s\n"%tran.to_gtf())
            if len(tran.exon_s)>0 :
                for i_exon,exon in enumerate(tran.exon_s) :
                    out.write(exon.to_gtf(exon_number=i_exon+1))
                    out.write("\n")
            elif len(tran.cds_s)> 0 :
                for i_cds,cds in enumerate(tran.cds_s) :
                    out.write(cds.to_gtf(to_exon=True,exon_number=i_cds+1))
                    out.write("\n")
            else :
                print "ERROR: No EXONS and CDSs"
    out.close()

    out_rep=file("%s.repr.xls"%(gtf_fn[:-4]),'wt')
    for gene in gene_s :
        max_tran=None
        for tran in gene.tran_s :
            if max_tran==None or max_tran.get_length()<tran.get_length() :
                max_tran=tran
        out_rep.write("%s\t%d\t%s\t%d\n"%\
                (gene.id(),len(gene.tran_s),tran.id(),tran.get_length()))
    out_rep.close()


class TestSequenceFunctions(unittest.TestCase):
    def setUp(self):
        pass
    def test_sort(self):
        chr_1=DNA("CM001265.1")
        chr_2=DNA("CM001266.1")
        chr_3=DNA("CM001272.1")
        chr_s=[chr_1,chr_2,chr_3]
        chr_s.sort()
        self.assertEqual(chr_s[0].name,"CM001265.1")
        self.assertEqual(chr_s[1].name,"CM001266.1")
        self.assertEqual(chr_s[2].name,"CM001272.1")
    def test_sort_2(self):
        chr_3=DNA("CM001265.1")
        chr_2=DNA("CM001266.1")
        chr_1=DNA("CM001272.1")
        chr_s=[chr_1,chr_2,chr_3]
        chr_s.sort()
        self.assertEqual(chr_s[0].name,"CM001265.1")
        self.assertEqual(chr_s[1].name,"CM001266.1")
        self.assertEqual(chr_s[2].name,"CM001272.1")


class TestConvertRoman(unittest.TestCase):
    def setUp(self):
        pass
    def test_run_non_roman(self):
        rst=ConvertRomanNumeralToNumber("CM001272.1")
        self.assertEqual(rst,None)
    def test_run_roman_2(self):
        rst=ConvertRomanNumeralToNumber("I3")
        self.assertEqual(rst,None)
    def test_run_roman(self):
        rst=ConvertRomanNumeralToNumber("I")
        self.assertEqual(rst,1)


if __name__=='__main__' :
    main()
    #unittest.main()
