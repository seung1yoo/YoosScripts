#!/usr/bin/python


import sys
import os
import copy

import Queue
import threading

import StringIO
import sqlite3 as lite
import socket
import ftplib

import unittest

CODING_S=('protein_coding','nonsense_mediated_decay','non_stop_decay',\
        'IG_V_gene','IG_J_gene','IG_C_gene','IG_D_gene',\
        'TR_V_gene','TR_C_gene','TR_J_gene','TR_D_gene',\
        'TEC','ambiguous_orf',\
        '3prime_overlapping_ncrna',\
        'sense_overlapping',\
        'processed_transcript',\
        )

NONCODING_S=('rRNA','tRNA','pseudogene','miRNA','ncRNA',\
        'transposable_element_gene','snoRNA','snRNA',\
        'misc_RNA','retrotransposed',\
        'retained_intron',\
        'Mt_tRNA','Mt_rRNA','lincRNA','non_coding',\
        'unprocessed_pseudogene','processed_pseudogene','polymorphic_pseudogene',\
        'transcribed_processed_pseudogene','transcribed_unprocessed_pseudogene',\
        'unitary_pseudogene',\
        'TR_V_pseudogene','TR_J_pseudogene',\
        'IG_V_pseudogene','IG_C_pseudogene','IG_J_pseudogene',\
        'antisense',\
        'sense_intronic',\
        )

T_R_RNA_S=('rRNA','tRNA','Mt_tRNA','Mt_rRNA')
LINC_RNA_S=('lincRNA',)
MI_RNA_S=('miRNA',)

MT_name_s = ('MT','Mt','mt','Mito','MtDNA')
PT_name_s = ('PT','Pt','pt')

FTP_ENS_ANIMAL_ULR="ftp.ensembl.org"
FTP_ENS_ENSGENOME_ULR="ftp.ensemblgenomes.org"

##################################################
# Gene Structures
TYPE_GS_DEFAULT = 0
TYPE_GS_KNOWN = 10
TYPE_GS_NOVEL = 20
TYPE_GS_MOD   = 30
TYPE_GS_MISSASSEMBLY =110
TYPE_GS_KNOWN_S=(TYPE_GS_KNOWN,TYPE_GS_MOD)
#
TYPE_GS_ALT_EXON_SKIPPING = 2
TYPE_GS_ALT_INTRON_RETENTION = 4
TYPE_GS_ALT_ALT_5 = 6
TYPE_GS_ALT_ALT_3 = 8
#
GS_GENE_TYPES=(TYPE_GS_DEFAULT,TYPE_GS_KNOWN,TYPE_GS_NOVEL,TYPE_GS_MOD)
GS_TRAN_TYPES=(TYPE_GS_DEFAULT,TYPE_GS_KNOWN,TYPE_GS_NOVEL)
#
##################################################
# Gene TYpe
TOOLS_HOME="/BiO/BioTools"
#
_TYPE_CODING_BASE           =10
TYPE_CODING_BASE            =10
TYPE_GENE_CODING            = 0
TYPE_GENE_CODING_COMPLETE   = 0
TYPE_GENE_CODING_5PARTIAL   = 1
TYPE_GENE_CODING_3PARTIAL   = 2
TYPE_GENE_CODING_INTERNAL   = 4
#
#
################################################################################
# External Programs
TOOLS_HOME="/BiO/BioTools"
#
SAMTOOLS_EXEC="%s/samtools/samtools-0.1.19/samtools"%TOOLS_HOME
BWA_EXEC="%s/bwa/bwa-0.7.10/bwa"%TOOLS_HOME
BOWTIE1_BUILD_EXEC="%s/bowtie/bowtie-1.1.1/bowtie-build"%TOOLS_HOME
BOWTIE2_BUILD_EXEC="%s/bowtie/bowtie2-2.2.3/bowtie2-build"%TOOLS_HOME
#
STAR_N_CPU=16
STAR_EXEC="%s/STAR/STAR_2.3.0e/STAR"%TOOLS_HOME
STAR_PATH="staridx"
#STAR_LIMITRAM="25000000000"
STAR_LIMITRAM="85359006421"
#
MONO_EXEC="%s/mono/mono-2.10.9/bin/mono"%TOOLS_HOME
FUSIONMAP_PATH="FusionMap"
FUSIONMAP_EXEC="%s/FusionMap/FusionMap_2014-01-01/bin/FusionMap.exe"%TOOLS_HOME
#
PICARD_CREATESEQUENCEDICT_EXEC="%s/Picard/picard-tools-1.119/CreateSequenceDictionary.jar"%TOOLS_HOME
#
BLAST_EXEC="%s/ncbi-blast-plus/ncbi-blast-2.2.29+/bin/makeblastdb" %TOOLS_HOME
################################################################################


class Model :
    def __init__(self):
        self._id=None
        self._name=None
        #
        self._chr=None
        self._strand=None
        self._start=None
        self._end=None
        #
        self._length=None
        self.biotype=None
        #
        self.kegg_s=None
        self.ncbi_gene_s=None
        self.go_with_anch_id_s=None
        self.go_id_s=None
        self.interpro_s=None
    def id(self):
        return self._id
    #
    def chr(self):
        return self._chr
    def strand(self):
        return self._strand
    def start(self):
        return self._start
    def end(self):
        return self._end
    # 
    def __len__(self):
        if self._length!=None and self._length!=0 :
            return self._length
        return abs(self.end()-self.start()+1)
    # 
    def repr_tran(self):
        tran_s=[]
        for tran in self.tran_s :
            tran_s.append( (len(tran),tran) )
        tran_s.sort()
        tran_s.reverse()

        for _len, tran in tran_s  :
            if tran.transl_id==None :
                #sys.stdout.write("INFO: %s;%s is non-coding\n"%(self.gene_id,tran.trans_name))
                continue
            return tran

        return tran_s[0][1]


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
    return number_result
 

class Foo :
    def __init__(self):
        self.tran_id=None
        self._id=None
        self.transl_id=None
        self.symbol='-'
        #
        self.go_id_s=None
        self.interpro_s=None
        self.refseq_s=None
        self.ncbi_gene_s=None
        self.protein=None
    def id(self):
        if self._id==None :
            return self.tran_id
        return self._id
    #
    def chr(self) :
        return self._chr
    def strand(self) :
        return self._strand
    def start(self):
        return self._start
    def end(self):
        return self._end
    def __len__(self):
        return abs(self._end-self._start+1)

class Runner(threading.Thread) :
    def __init__(self,job_que,rst_que) :
        threading.Thread.__init__(self)
        self.job_que=job_que
        self.rst_que=rst_que
    def run(self):
        while True :
            job=self.job_que.get()
            if job=='END' :
                break
            i_cmd,cmd,stderr,stdout,work_home=job
            old_cwd=os.getcwd()
            if work_home!=None :
                os.chdir(work_home)
            p=None
            if stderr!=None and stdout!=None :
                p=Popen(cmd,stderr=stderr,stdout=stdout)
            elif stderr!=None :
                p=Popen(cmd,stderr=stderr)
            elif stdout!=None :
                p=Popen(cmd,stdout=stdout)
            else :
                p=Popen(cmd)
            out,err=p.communicate()
            if work_home!=None :
                os.chdir(old_cwd)
            self.rst_que.put( (i_cmd,cmd,out,err) )
        self.rst_que.put("END")

class Reporter(threading.Thread) :
    def __init__(self,n_proc,rst_que) :
        threading.Thread.__init__(self)
        self.n_proc=n_proc
        self.rst_que=rst_que
        self.rst_s=[]
    def run(self):
        n_finished=0
        while True :
            if n_finished==self.n_proc :
                break
            job=self.rst_que.get()
            if job=='END' :
                n_finished+=1
                continue
            i_cmd,cmd,out,err=job
            self.rst_s.append((i_cmd,cmd,out,err))

class Executor :
    def __init__(self,exec_type,n_cpu,path_finder=None):
        self.exec_type=exec_type
        self._n_cpu=n_cpu
        self._que=[]
        self._log=[]
        if path_finder==None :
            self._path="%s/scripts"%os.getcwd()
            self._log_path="%s/logs"%os.getcwd()
        else :
            self._path="%s/scripts"%(path_finder.home())
            self._log_path="%s/logs"%(path_finder.home())
        self._logger=None
    def init(self):
        self._que=[]
        self._log=[]
    def n_cpu(self):
        return self._n_cpu
    def add(self,cmd,stderr=None,stdout=None,work_home=None) :
        if not os.path.exists(self._path) :
            os.makedirs(self._path)
        if not os.path.exists(self._log_path) :
            os.makedirs(self._log_path)
        self._que.append( (cmd,stderr,stdout,work_home) )
        return len(self._que)-1
    def run(self) :
        self._log=[]
        if self.exec_type==EXEC_TYPE_CLI :
            self._run_cli()
        elif self.exec_type==EXEC_TYPE_THREADS :
            self._run_threads()
        else :
            return False
        self._que=[]
        return True
    def _run_cli(self) :
        i_cmd=-1
        for cmd,stderr,stdout,work_home in self._que :
            i_cmd+=1
            old_cwd=os.getcwd()
            if work_home!=None :
                os.chdir(work_home)
            p=None
            if stderr!=None and stdout!=None :
                p=Popen(cmd,stderr=stderr,stdout=stdout)
            elif stderr!=None :
                p=Popen(cmd,stderr=stderr)
            elif stdout!=None :
                p=Popen(cmd,stdout=stdout)
            else :
                p=Popen(cmd)
            out,err=p.communicate()
            if work_home!=None :
                os.chdir(old_cwd)
            self._log.append( (i_cmd,cmd,out,err) )
        return True

    def _run_threads(self):
        job_que=Queue.Queue(self._n_cpu*2)
        rst_que=Queue.Queue(self._n_cpu*2)
        #
        runner_s=[]
        for i_thread in range(self._n_cpu) :
            runner=Runner(job_que,rst_que)
            runner.start()
            runner_s.append(runner)
        #
        reporter=Reporter(self._n_cpu,rst_que)
        reporter.start()
        #
        i_cmd=-1
        for cmd,stderr,stdout,work_home in self._que :
            i_cmd+=1
            job_que.put( (i_cmd,cmd,stderr,stdout,work_home) )
        for i_thread in range(self._n_cpu) :
            job_que.put('END')
        #
        for runner in runner_s :
            runner.join()
        reporter.join()

        self._log.extend(reporter.rst_s)
            
        return True
    def wait(self) :
        return self._log

class DnaGroup(list):
    def ref_dna_s(self):
        has_ref=False
        has_chr=False
        for dna in self :
            if dna.is_ref() :
                has_ref=True
                break
        for dna in self :
            if dna.is_chromosome() :
                has_chr=True
                break
        #
        if has_ref and has_chr :
            dna_s=[]
            for dna in self :
                if dna.is_ref() and dna.is_chromosome() :
                    dna_s.append(dna)
            return dna_s
        elif has_ref :
            dna_s=[]
            for dna in self :
                if dna.is_ref() :
                    dna_s.append(dna)
            return dna_s
        else :
            return self
    def major_s(self):
        has_chr=False
        has_ref=False
        for dna in self :
            if dna.is_chromosome() :
                has_chr=True
                break
        for dna in self :
            if dna.is_ref() :
                has_ref=True
                break
        sel=[]
        for dna in self :
            if dna.is_ref()!=has_ref :
                continue
            if dna.is_chromosome()!=has_chr :
                continue
            sel.append(dna)
        return sel

class DNA :
    def __init__(self,name,level):
        self.name=name
        self.level=level
        self._int_name=None
        self.type=None
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
        if self.level=='chromosome' :
            return True
        return False
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
    
    if unit_s[2]=='gene' :
        return

    foo=Foo()
    foo._chr=unit_s[0]
    foo.source=unit_s[1]
    foo._start=int(unit_s[3])
    foo._end=int(unit_s[4])
    foo._strand=unit_s[6]

    foo.gene_id=None
    foo.gene_name=None
    foo.gene_biotype=None
    foo.transcript_biotype="Unknown"
    foo._symbol=None

    foo.tran_id=None
    foo.exon_number=None
    
    foo.trans_name=None
    foo.transl_name=None

    foo.attr_s={}
    for attr in parse_attr_s(unit_s[8]) :
        attr=attr.strip()
        if attr=='' :
            continue
        key,val=attr.split(None,1)
        val=val[1:-1]
        if key=='gene_id' :
            foo.gene_id=val
        elif key=='gene_accession' :
            foo.gene_name=val
        elif key=='gene_name' :
            foo._symbol=val
        elif key=='transcript_id' :
            foo.tran_id=val
        elif key=='transcript_accession' :
            foo._name=val
        elif key=='exon_number' :
            foo.exon_number=int(val)
        elif key=='gene_biotype' :
            foo.gene_biotype=val
        elif key=='transcript_biotype' :
            foo.transcript_biotype=val
        foo.attr_s[key]=val

    if foo.gene_biotype==None \
            and (foo.source in CODING_S or foo.source in NONCODING_S) :
        foo.gene_biotype=foo.source
    if foo.transcript_biotype=="Unknown" \
            and (foo.source in CODING_S or foo.source in NONCODING_S) :
        foo.transcript_biotype=foo.source

    foo.line=line

    return foo

def _update_transcript_position(tran) :
    tran._start=tran.exon_s[0]._start
    tran._end=tran.exon_s[0]._end
    tran._length=0
    for exon in tran.exon_s :
        if tran._start > exon._start :
            tran._start=exon._start
        if tran._end < exon._end :
            tran._end=exon._end
        tran._length+=len(exon)

def _read_gtf_transcript(_in) :
    tran=None
    while True :
        line=_in.readline()
        if line=='' :
            break
        if line.startswith("#") :
            continue
        #
        exon=read_gtf_line(line,'exon')
        if exon==None :
            continue
        #
        if tran!=None and tran.tran_id==exon.tran_id  :
            tran.exon_s.append(exon)
            continue
        if tran!=None :
            _update_transcript_position(tran)
            _update_transcript_coding_type(tran)
            yield tran
        #
        tran=copy.copy(exon)
        tran.source=exon.source
        tran._id=tran.tran_id
        tran._alt_id=tran.tran_id
        tran.protein=None
        tran.exon_s=[exon]
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
        return
    if tran.source=='protein_coding' :
        tran._coding_type=TYPE_GENE_CODING
    elif tran.source=='lincRNA' :
        tran._coding_type=TYPE_GENE_NONCODING_LINC
    elif tran.source=='tRNA' or tran.source=='MT_tRNA' :
        tran._coding_type=TYPE_GENE_NONCODING_TRNA
    elif tran.source=='rRNA' or tran.source=='MT_rRNA' :
        tran._coding_type=TYPE_GENE_NONCODING_RRNA
    elif tran.source.find("pseudogene")>-1 \
            or tran.source in NONCODING_S :
        tran._coding_type=TYPE_GENE_NONCODING
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
        elif tran.source.find("pseudogene")>-1 \
                or tran.source in NONCODING_S :
            has_nc=True
            tran._coding_type=TYPE_GENE_NONCODING
        elif tran.source.startswith("IG_") or  tran.source.startswith("TR_")  :
            has_coding=True
            tran._coding_type=TYPE_GENE_CODING
            break
            
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




def read_gtf_gene(fn) :
    _in=file(fn)
    gene=None
    for tran in _read_gtf_transcript(_in) :
        if gene!=None and tran.gene_id==gene.gene_id :
            gene.tran_s.append(tran)
            continue
        #
        if gene!=None :
            _update_gene_position(gene)
            _update_gene_coding_type(gene)
            yield gene
        #
        gene=Model()
        gene._id=tran.gene_id
        gene._alt_id=tran.gene_id
        gene.source=tran.source
        gene.gene_biotype=tran.gene_biotype
        #
        gene._name=tran.gene_name
        if gene._name==None :
            gene._name=tran.gene_id
        gene._symbol=tran._symbol
        #
        gene._chr=tran._chr
        gene._strand=tran._strand
        gene.symbol=None
        gene.status=None
        #
        gene._mol_type=TYPE_GS_DEFAULT
        gene._coding_type=TYPE_GENE_CODING
        gene.gene_id=tran.gene_id
        gene.tran_s=[tran]
        gene.desc=None
        gene.go_id_s=[]
        gene.interpro_s=[]
        gene.ncbi_gene_s=None
        #
    if gene!=None :
        _update_gene_position(gene)
        _update_gene_coding_type(gene)
        yield gene


class DnaGroupTest (unittest.TestCase):
    def test_group(self):
        pass

    pass

if __name__=='__main__' :
    main()

