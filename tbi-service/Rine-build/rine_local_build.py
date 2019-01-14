#!/usr/bin/python

import os
import sys
import gzip
import time
import unittest
import MySQLdb
import copy
import getopt

import pdb

from rine_liblocal import *

#REFERENCE_HOME="/BiOfs/BioResources/References"
#REFERENCE_HOME="/BiO/References"

REF_TYPE="ensembl"

SPECIES=['animal','PLANTS','BACTERIA','FUNGI']

INFO_DOWN_FN="info.down"
INFO_DESC_FN="info.desc"
INFO_REPEAT_FN="info.repeat"
INFO_IDX_FN="info.index"
INFO_FUSION="info.fusion"

DB_PATH="db"

CHR_PATH="chr.ensembl"
CHR_OUT_PATH="chr"

TXT_s=["gene.txt.gz","transcript.txt.gz","translation.txt.gz",\
        "exon.txt.gz","exon_transcript.txt.gz","protein_feature.txt.gz",\
        "external_db.txt.gz","interpro.txt.gz","object_xref.txt.gz",\
        "ontology_xref.txt.gz","xref.txt.gz","analysis.txt.gz",\
        "repeat_consensus.txt.gz","repeat_feature.txt.gz","seq_region.txt.gz",\
        "meta.txt.gz"]

DB_EXTRA_s=["exon_transcript","external_db","interpro",\
        "object_xref","ontology_xref","xref","analysis",\
        "repeat_consensus","repeat_feature","seq_region","meta"]

PROTEIN_FEATURE_S=["pfam","superfamily","Smart","scanprosite",\
        "uniprot_mammal","uniprot_vertebrate","uniprot_non_vertebrate",\
        "prints","tmhmm","seg"]

BOWTIE_PATH="bowtie_idx"

ENS_GENE2TRANSCRIPT_FN="ens_gene2transcript.txt"
ENS_GENE_DESC_FN="ens_gene.desc.txt"
ENS_GOA_GENE_FN="ens_go.gene.txt"
ENS_GOA_TRANSCRIPT_FN="ens_go.transcript.txt"

SELECT_S=("X","Y","W","Z","MT","Mt","M","PT","Pt","P")

GO_DB_HOST="wolf.myomics.net"
GO_DB_USER="giomics"
GO_DB_PS="giomics"
GO_DB_NAME="GO_201312"

NCBI_GENE_DB_HOST="wolf.myomics.net"
NCBI_GENE_DB_NAME="NCBI_GENE_20131216"

KEGG_DB="/BiO/BioResources/DBs/KeggPathway/KEGG.sqlite.current"

ENS_GENOMES=('PLANTS','BACTERIA','FUNGI')

class PathFinder :
    def __init__(self,species,release,name_s) :
        self.species=species
        self.release=release
        self.name_s=name_s
        #
        self.db_s=[]
        #
        self.dna_s=[]
        self.core_db=None
    #
    def home(self):
        formed_name=self.ens_formatted_species()
        new_name="%s%s"%(formed_name[0].upper(),formed_name[1:])
        return "%s/%s/%s"%(self._ref_home,new_name,self.short_name())
    def to_rel_path(self,path):
        return path.replace(self.home(),"[PATH]")
    def fusionmap_path(self):
        return "%s/%s"%(self.home(),FUSIONMAP_PATH)
    def ens_formatted_species(self):
        buf=[]
        for name in self.name_s :
            buf.append(name.lower().replace("-","_"))
        return "_".join(buf)
    #
    def meta_path(self):
        return "%s/meta"%(self._ref_home)
    def meta_list_fn(self):
        return "%s/ens_genomes_%s"%(self.meta_path(),self.release)
    #
    def collection_path(self) :
        buf=[]
        for tag in self.core_db.split("_") :
            if tag=='core' :
                break
            buf.append(tag)
        return "_".join(buf)
    #
    def staridx_path(self):
        return "%s/%s"%(self.home(),STAR_PATH)
    def chr_path(self):
        return "%s/%s"%(self.home(),CHR_OUT_PATH)
    def chr_source_path(self):
        return "%s/%s"%(self.home(),CHR_PATH)
    #
    def formal_full_name(self):
        return "%s%s"%(self.species[0].upper(),self.species[1:])
    def short_name(self):
        buf=[]
        buf.append(self.species[0].upper())
        buf.append("_")
        buf.append(self.species.split("_",1)[1])
        buf.append("_")
        buf.append(self._source)
        buf.append("_")
        buf.append(self.release)
        short_name="".join(buf)
        return short_name
    def _prefix(self):
        return "%s/%s"%(self.home(),self.short_name())
    def prefix(self,tag=None):
        if tag==None :
            return "%s/%s.chr"%(self.home(),self.short_name())
        return "%s/%s.chr.%s"%(self.home(),self.short_name(),tag)
    #
    def fa_list_fn(self,fa_fn) :
        return "%s.list"%fa_fn
    #
    def chr_fn(self) :
        return "%s/%s.chr.fa"%(self.home(),self.short_name())
    def rna_fn(self) :
        return "%s/%s.rna.fa"%(self.home(),self.short_name())
    def rna_coding_fn(self) :
        return "%s/%s.rna.coding.fa"%(self.home(),self.short_name())
    def ncrna_fn(self) :
        return "%s/%s.rna.noncoding.fa"%(self.home(),self.short_name())
    def pep_fn(self) :
        return "%s/%s.pep.fa"%(self.home(),self.short_name())
    #
    def gtf_fn(self,tag=None):
        return "%s.gtf"%(self.prefix(tag))
    def non_gtf_fn(self,tag=None):
        return "%s.removed.gtf"%(self.prefix(tag))
    def repr_gtf(self,tag=None):
        return "%s.repr.gtf"%self.prefix(tag)
    def coding_gtf(self,tag=None):
        return "%s.coding.gtf"%self.prefix(tag)
    def coding_repr_gtf(self,tag=None):
        return "%s.coding.repr.gtf"%self.prefix(tag)
    def non_coding_all_gtf(self,tag=None):
        return "%s.ncRNA.all.gtf"%self.prefix(tag)
    def non_coding_repr_gtf(self,tag=None):
        return "%s.ncRNA.repr.gtf"%self.prefix(tag)
    def t_r_rna_gtf(self,tag=None):
        return "%s.tnrRNA.gtf"%self.prefix(tag)
    def mt_rna_gtf(self,tag=None):
        return "%s.mtRNA.gtf"%self.prefix(tag)
    def pt_rna_gtf(self,tag=None):
        return "%s.ptRNA.gtf"%self.prefix(tag)
    def t_r_mt_pt_rna_gtf(self,tag=None):
        return "%s.TnRnMTnPTRNA.gtf"%self.prefix(tag)
    def mirna_gtf(self,tag=None):
        return "%s.miRNA.gtf"%self.prefix(tag)
    def linc_gtf(self,tag=None):
        return "%s.lincRNA.gtf"%self.prefix(tag)
    def mask_gtf(self,tag=None):
        return "%s.mask.gtf"%self.prefix(tag)
    def t_r_mt_pt_rna_nc_gtf(self,tag=None):
        return "%s.TnRnMTnPTRNAnNC.gtf"%self.prefix(tag)
    #
    def repeat_bed_fn(self):
        return "%s.repeat.bed"%(self.prefix())
    #
    def gene_id_info(self):
        return "%s/info.id"%(self.home())
    def gene_id_coding_info(self):
        return "%s/info.id.coding"%(self.home())
    #
    def gene_tran_map(self):
        return "%s.map"%self.prefix()
    def gene_to_tran(self):
        return "%s/%s"%(self.home(),ENS_GENE2TRANSCRIPT_FN)
    def gene_info(self):
        return "%s.genes.info"%self._prefix()
    def gene_feature(self):
        return "%s.genes.feature"%self._prefix()
    #
    def tran_info(self):
        return "%s.transcripts.info"%self._prefix()
    def tran_feature(self):
        return "%s.transcripts.feature"%self._prefix()
    #
    def info_down_fn(self):
        return "%s/%s"%(self.home(),INFO_DOWN_FN)
    def info_desc_fn(self):
        return "%s/%s"%(self.home(),INFO_DESC_FN)
    def info_repeat_fn(self):
        return "%s/%s"%(self.home(),INFO_REPEAT_FN)
    def info_idx_fn(self):
        return "%s/%s"%(self.home(),INFO_IDX_FN)
    def info_fusion_fn(self):
        return "%s/%s"%(self.home(),INFO_FUSION)
    def info_geneset(self):
        return "%s/info.geneset"%(self.home())
    def info_geneset_exonic(self):
        return "%s/info.geneset.exonic"%(self.home())
    def info(self):
        return "%s/ref_%s"%(self.home(),self.short_name())
    #
    def info_build_fn(self):
        return "%s/info.build"%(self.home())
    #
    def db_path(self):
        return "%s/db"%(self.home())
    def db_fn(self):
        return "%s/%s.ensembl.db"%(self.db_path(),self.species)
    def db_done(self):
        return "%s/%s.ensembl.db.done"%(self.db_path(),self.species)

class Run :
    def __init__(self) :
        self.species=None

class NcbiGeneDatabase :
    def __init__(self):
        self.db=None
        self.cursor=None
        self.db = MySQLdb.connect(db=GO_DB_NAME,user=GO_DB_USER,passwd=GO_DB_PS,host=GO_DB_HOST)
        self.cursor = self.db.cursor()
    def __del__(self):
        if self.cursor!=None :
            self.cursor.close()
        if self.db!=None :
            self.db.close()
    def search_protein_acc(self,acc):
        buf=[]
        buf.append("SELECT")
        buf.append("ga.GeneId, ga.rna_acc, ga.protein_acc")
        buf.append(", gi.Symbol, gi.description")
        buf.append("FROM gene2accession ga")
        buf.append("LEFT JOIN gene_info AS gi")
        buf.append("ON gi.GeneId = ga.GeneId")
        buf.append("WHERE ga.protein_acc=%s")
        buf.append("GROUP BY ga.GeneId")
        sql=" ".join(buf)
        #
        rst=self.cursor.execute(sql,(acc,))
        if rst==0 :
            return
        rst=self.cursor.fetchone()
        #
        foo=Foo()
        foo._db_id=rst[0]
        foo.rna_acc=rst[1]
        foo.acc=rst[2]
        foo.protein_acc=rst[2]
        foo.symbol=rst[3]
        foo.desc=rst[4]
        return foo

class GoDB :
    def __init__(self):
        self.db=None
        self.cursor=None
        self.db = MySQLdb.connect(db=GO_DB_NAME,user=GO_DB_USER,passwd=GO_DB_PS,host=GO_DB_HOST)
        self.cursor = self.db.cursor()
    def __del__(self):
        if self.cursor!=None :
            self.cursor.close()
        if self.db!=None :
            self.db.close()
    def search_ancestor(self,acc_id) :
        buf=[]
        buf.append("SELECT DISTINCT")
        buf.append("ancestor.acc")
        buf.append(", graph_path.distance")
        buf.append("FROM term")
        buf.append("INNER JOIN graph_path ON (term.id=graph_path.term2_id)")
        buf.append("INNER JOIN term AS ancestor")
        buf.append("ON (ancestor.id=graph_path.term1_id)")
        buf.append("WHERE term.acc='%s'"%acc_id)
        sql=" ".join(buf)
        rst=self.cursor.execute(sql)
        if rst==0 :
            return []
        go_id_s=set()
        for rst in self.cursor.fetchall() :
            go_id,dist=rst
            if go_id not in go_id_s :
                go_id_s.add(go_id)
        go_id_s=list(go_id_s)
        go_id_s.sort()
        return go_id_s

class BufferedGoDB :
    def __init__(self):
        self.go_db=GoDB()
        self._local_id_s=set()
        self._local={}
        self._local_anch_id_s=set()
        self._local_anch={}
    def category(self):
        return self.go_db.category()
    def search(self,go_id) :
        if go_id in self._local_id_s :
            return copy.copy(self._local[go_id])
        go=self.go_db.search(go_id)
        self._local_id_s.add(go_id)
        self._local[go_id]=go
        return go
    def search_ancestor(self,go_id) :
        if go_id in self._local_anch_id_s :
            return copy.copy(self._local_anch[go_id])
        go_id_s=self.go_db.search_ancestor(go_id)
        self._local_anch_id_s.add(go_id)
        self._local_anch[go_id]=go_id_s
        return go_id_s


class Downloads :
    def __init__(self):
        self.ncrna_fn=None
    def uncompress(self) :
        if self.dna_fn.endswith(".gz") :
            if os.path.exists(self.dna_fn) :
                os.system("gzip -d %s"%(self.dna_fn))
            self.dna_fn=self.dna_fn[:-3]
        #
        if self.cdna_fn.endswith(".gz") :
            if os.path.exists(self.cdna_fn) :
                os.system("gzip -d %s"%(self.cdna_fn))
            self.cdna_fn=self.cdna_fn[:-3]
        #
        if self.ens_gtf_fn.endswith(".gz") :
            if os.path.exists(self.ens_gtf_fn) :
                os.system("gzip -d %s"%(self.ens_gtf_fn))
            self.ens_gtf_fn=self.ens_gtf_fn[:-3]


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


class Meta :
    def __init__(self,path_finder):
        self._path_finder=path_finder
    def get_coding_id_s(self):
        id_s=[]
        for line in file(self._path_finder.gene_id_coding_info()) :
            id_s.append(line.strip().split()[0])
        return id_s
    def get_coding_tran_id_s(self):
        id_s=[]
        for line in file(self._path_finder.gene_id_coding_info()) :
            id_s.extend(line.strip().split()[2].split(","))
        return id_s


def prepare_path(path_finder) :
    if not os.path.exists(path_finder.home()) :
        os.makedirs(path_finder.home())
    if not os.path.exists(path_finder.chr_path()) :
        os.makedirs(path_finder.chr_path())
    if not os.path.exists(path_finder.chr_path()) :
        sys.stderr.write("ERROR: Failed to make out_path.\n")
        return False
    return True


def _gen_cdna_fn(gtf_fn,chr_fn,rna_fn) :
    os.system("gtf_to_fasta %s %s %s.raw"%(gtf_fn,chr_fn,rna_fn))
    out=file(rna_fn,'wt')
    for line in file("%s.raw"%rna_fn) :
        if line.startswith(">") :
            out.write(">%s"%(line.split(None,1)[1]))
        else :
            out.write(line)
    out.close()


def rebuild_cdna_fn(path_finder) :
    if os.path.exists(path_finder.rna_fn()) \
            and os.path.getsize(path_finder.rna_fn())>1000 :
        sys.stdout.write("Found RNA_FN(%s).\n"%(path_finder.rna_fn()))
        return True
    sys.stdout.write("Starting to rebuild RNA_FN.\n")
    #
    _gen_cdna_fn(path_finder.gtf_fn(),path_finder.chr_fn(),\
            path_finder.rna_fn())
    _gen_cdna_fn(path_finder.coding_gtf(),path_finder.chr_fn(),\
            path_finder.rna_coding_fn())
    _gen_cdna_fn(path_finder.non_coding_all_gtf(),path_finder.chr_fn(),\
            path_finder.ncrna_fn())
    #
    return True

# ------------------------------------------------------------------------------
# Download
#
def download_meta(ftp,path_finder,species) :
    sys.stdout.write("INFO: Downloading META RELEASE FILE\n")
    path=ftp.path(rel_top=True)
    ftp.download("%s/species.txt"%path,path_finder.meta_list_fn())

def download_db(ftp,path_finder,species) :
    local_path="%s/db"%(path_finder.home())
    if not os.path.exists(local_path) :
        os.makedirs(local_path)

    path=ftp.path("mysql")
    sel_path=None
    if path_finder.type in ENS_GENOMES :
        sel_path="%s/%s"%(path,path_finder.core_db)
    else :
        for line in ftp.list(path) :
            fn=line.strip().split()[-1]
            if fn.startswith("%s_core"%(species.lower())) :
                sel_path="%s/%s"%(path,fn)
                break

    if sel_path==None :
        sys.stderr.write("ERROR\tNo selected path for mysql (%s,%s)\n"%(species,path))
        return False
    #
    sys.stdout.write("LIST %s\n"%sel_path)
    sql_selected=None
    for line in ftp.list(sel_path) :
        fn=line.strip().split()[-1]
        if fn.endswith(".sql.gz") :
            sql_selected=fn
            break
    if sql_selected==None :
        sys.stderr.write("ERROR: No SQL_FN.\n")
        return False
    #
    sql_fn="%s/%s"%(local_path,sql_selected)
    sys.stdout.write("Download SQL %s\n"%sql_fn)
    ftp.download("%s/%s"%(sel_path,sql_selected),sql_fn)
    path_finder.db_s.append(("SQL",sql_fn))
    #
    for txt_fn in TXT_s :
        local_fn="%s/%s"%(local_path,txt_fn)
        path_finder.db_s.append((txt_fn,local_fn))
        if os.path.exists(local_fn) :
            msg="INFO: Found previous downloaded file (%s)\n"%local_fn
            sys.stdout.write(msg)
            continue
        sys.stdout.write("INFO: Starting to download %s.\n"%txt_fn)
        ftp.download("%s/%s"%(sel_path,txt_fn),local_fn)

def download_gtf(ftp,path_finder,species) :
    path=ftp.path("gtf")
    sel_path=None
    #
    if path_finder.type=='BACTERIA' :
        if path_finder.core_db!=None :
            path="%s/%s"%(path,path_finder.collection_path())
        sys.stdout.write("INFO: GTF: BACTERIA: %s (%s)\n"%(path,path_finder.core_db))
        sub_path=None
        name_s=species.split("_")
        prefix=path_finder.ens_formatted_species()
        sel_path="%s/%s"%(path,path_finder.ens_formatted_species())
    elif path_finder.type in ENS_GENOMES :
        if path_finder.core_db!=None :
            path="%s/%s"%(path,path_finder.division)
        sys.stdout.write("INFO: GTF: %s (%s)\n"%(path,path_finder.core_db))
        sub_path=None
        name_s=species.split("_")
        prefix=path_finder.ens_formatted_species()
        sel_path=path
    else :
        for line in ftp.list(path) :
            fn=line.strip().split()[-1]
            if fn.startswith(species.lower()) :
                sel_path="%s/%s"%(path,fn)
                break
    if sel_path==None :
        sys.stderr.write("ERROR: %s: gtf\n"%species)
        return False
    #
    sys.stdout.write("LIST %s\n"%sel_path)
    gtf_selected=None
    for line in ftp.list(sel_path) :
        fn=line.strip().split()[-1]
        if fn.endswith(".gtf.gz") :
            gtf_selected=fn
            break
    if gtf_selected==None :
        sys.stderr.write("ERROR: No GTF_FN.\n")
        return False
    #
    gtf_fn="%s/%s"%(path_finder.home(),gtf_selected)
    if os.path.exists(gtf_fn) :
        sys.stdout.write("PASS: Found previous downloaded files (%s)\n"%gtf_fn)
        path_finder.ens_gtf_fn=gtf_fn
    elif gtf_fn.endswith(".gz") and os.path.exists(gtf_fn[:-3]) :
        sys.stdout.write("PASS: Found previous downloaded files (%s)\n"%gtf_fn)
        path_finder.ens_gtf_fn=gtf_fn[:-3]
    else :
        sys.stdout.write("Downloading GTF %s\n"%gtf_fn)
        ftp.download("%s/%s"%(sel_path,gtf_selected),gtf_fn)
        path_finder.ens_gtf_fn=gtf_fn

    return True

def download_fasta_dna(ftp,path_finder,species,release) :
    path=ftp.path("fasta")
    sel_path=None
    if path_finder.type=='BACTERIA' :
        if path_finder.core_db!=None :
            path="%s/%s"%(path,path_finder.collection_path())
        sys.stdout.write("INFO: BACTERIA: %s\n"%path)
        sub_path=None
        name_s=species.split("_")
        prefix=path_finder.ens_formatted_species()
        sel_path="%s/%s/dna"%(path,path_finder.ens_formatted_species())
    elif path_finder.type in ENS_GENOMES :
        sel_path="%s/%s/dna"%(path,path_finder.division)
    else :
        for line in ftp.list(path) :
            fn=line.strip().split()[-1]
            if fn.startswith(path_finder.species.lower()) :
                sel_path="%s/%s/dna"%(path,fn)
                break
    if sel_path==None :
        sys.stderr.write("ERROR: %s: fasta_dna\n"%\
                (path_finder.species))
        return False
    #
    sys.stdout.write("FASTA_DNA: LIST %s\n"%sel_path)
    fn_s=[]
    for line in ftp.list(sel_path) :
        fn=line.strip().split()[-1]
        if fn.find(".%s.dna."%release)>-1 :
            fn_s.append(fn)
    #
    if not os.path.exists(CHR_PATH) :
        os.makedirs(CHR_PATH)
    #
    for fn in fn_s :
        out_fn="%s/%s"%(path_finder.chr_source_path(),fn)
        if fn.endswith("toplevel.fa.gz") :
            out_fn="%s/%s"%(path_finder.home(),fn)
        else :
            continue
        #
        if os.path.exists(out_fn) :
            sys.stdout.write("PASS\tFound previous files (%s)\n"%out_fn)
            path_finder.dna_s.append((fn,out_fn))
        elif out_fn.endswith(".gz") and os.path.exists(out_fn[:-3]) :
            out_fn=out_fn[:-3]
            sys.stdout.write("PASS\tFound previous files (%s)\n"%out_fn)
            path_finder.dna_s.append((fn,out_fn))
        else :
            sys.stdout.write("INFO\tStarting to download %s.\n"%fn)
            ftp.download("%s/%s"%(sel_path,fn),out_fn)
            path_finder.dna_s.append((fn,out_fn))
    return True

def download_fasta_cdna(ftp,path_finder,species,release) :
    path=ftp.path("fasta")
    sel_path=None
    if path_finder.type=='BACTERIA' :
        if path_finder.core_db!=None :
            path="%s/%s"%(path,path_finder.collection_path())
        sys.stdout.write("INFO: BACTERIA: %s\n"%path)
        sub_path=None
        name_s=species.split("_")
        prefix=path_finder.ens_formatted_species()
        sel_path="%s/%s/cdna"%(path,path_finder.ens_formatted_species())
    elif path_finder.type in ENS_GENOMES :
        sel_path="%s/%s/cdna"%(path,path_finder.division)
    else :
        for line in ftp.list(path) :
            fn=line.strip().split()[-1]
            if fn.startswith(path_finder.species.lower()) :
                sel_path="%s/%s/cdna"%(path,fn)
                break
    if sel_path==None :
        sys.stderr.write("ERROR: %s: fasta_cdna\n"%\
                (path_finder.species))
        return False
    #
    sys.stdout.write("LIST %s\n"%sel_path)
    selected=None
    for line in ftp.list(sel_path) :
        fn=line.strip().split()[-1]
        if fn.endswith(".cdna.all.fa.gz") :
            selected=fn
            break
    if selected==None :
        sys.stderr.write("ERROR: No CDNA_FN.\n")
        return False
    #
    out_fn="%s/%s"%(path_finder.home(),selected)
    if os.path.exists(out_fn) :
        sys.stdout.write("PASS\tFound previous files (%s)\n"%out_fn)
        path_finder.ens_cdna_fn=out_fn
    elif out_fn.endswith(".gz") and os.path.exists(out_fn[:-3]) :
        out_fn=out_fn[:-3]
        sys.stdout.write("PASS\tFound previous files (%s)\n"%out_fn)
        path_finder.ens_cdna_fn=out_fn
    else :
        sys.stdout.write("Downloading CDNA %s\n"%selected)
        ftp.download("%s/%s"%(sel_path,selected),out_fn)
        path_finder.ens_cdna_fn=out_fn

    return True

def download_fasta_ncrna(ftp,path_finder,species,release) :
    path_finder.ens_ncrna_fn=None
    #
    path=ftp.path("fasta")
    sel_path=None
    if path_finder.type=='BACTERIA' :
        if path_finder.core_db!=None :
            path="%s/%s"%(path,path_finder.collection_path())
        sys.stdout.write("INFO: BACTERIA: %s\n"%path)
        sub_path=None
        name_s=species.split("_")
        prefix=path_finder.ens_formatted_species()
        sel_path="%s/%s/ncrna"%(path,path_finder.ens_formatted_species())
    elif path_finder.type in ENS_GENOMES :
        sel_path="%s/%s/ncrna"%(path,path_finder.division)
    else :
        for line in ftp.list(path) :
            fn=line.strip().split()[-1]
            if fn.startswith(path_finder.species.lower()) :
                sel_path="%s/%s/ncrna"%(path,fn)
                break
    if sel_path==None :
        sys.stderr.write("ERROR: %s: fasta_ncrna\n"%\
                (path_finder.species))
        return False
    #
    sys.stdout.write("LIST %s\n"%sel_path)
    selected=None
    for line in ftp.list(sel_path) :
        fn=line.strip().split()[-1]
        if fn.endswith(".ncrna.fa.gz") :
            selected=fn
            break
    if selected==None :
        sys.stderr.write("WARN: No NCRNA_FN.\n")
        return False
    #
    out_fn="%s/%s"%(path_finder.home(),selected)
    if out_fn.endswith(".gz") :
        out_fn=out_fn[:-3]
    if not os.path.exists(out_fn) :
        sys.stdout.write("Downloading NCRNA %s\n"%selected)
        ftp.download("%s/%s"%(sel_path,selected),"%s.gz"%out_fn)
    else :
        sys.stdout.write("PASS: Found previous files (%s)\n"%out_fn)
    #
    if out_fn.endswith(".gz") :
        os.system("gzip -d %s"%out_fn)
    path_finder.ens_ncrna_fn=out_fn[:-3]
 
#    out_fn=path_finder.ncrna_fn()
#    if not os.path.exists(out_fn) :
#        sys.stdout.write("Downloading NCRNA %s\n"%selected)
#        ftp.download("%s/%s"%(sel_path,selected),"%s.gz"%out_fn)
#    else :
#        sys.stdout.write("PASS: Found previous files (%s)\n"%out_fn)
#    #
#    os.system("gzip -d -f %s.gz"%out_fn)
#    path_finder.ens_ncrna_fn=out_fn
    #
    return True

def update_sub_path(path_finder,ftp) :
    division=None
    sub=None
    for line in file(path_finder.meta_list_fn()) :
        if line.startswith("#") :
            continue
        unit_s=line.strip().split("\t")
        name=unit_s[0].replace(" ","_").lower()
        if name==path_finder.ens_formatted_species() :
            division=unit_s[1]
            sub=unit_s[12]
            break
    if sub!=None :
        ftp.core_db=sub
        path_finder.division=division
        path_finder.core_db=sub
        sys.stdout.write("INFO: CORE_DB: %s: %s\n"%\
                (path_finder.ens_formatted_species(),sub))
        return True
    else :
        sys.stdout.write("INFO: No CORE_DB: %s\n"%(path_finder.ens_formatted_species()))
        return False

def download_files(path_finder,run) :
    sys.stdout.write("INFO: Starting COPY (%s).\n"%(run.gtf_fn))
    short_name=path_finder.short_name()

    out=file(path_finder.info_down_fn(),'wt')
    out.write("RELEASE\t%s\n"%run.release)
    out.write("SPECIES\t%s\n"%run.species)
    out.write("SHORT_NAME\t%s\n"%short_name)
    out.write("GTF\t%s\n"%run.gtf_fn)
    out.write("DNA\t%s\n"%run.genome_fa)
    out.close()

    return True
# ------------------------------------------------------------------------------
#  Geneset
#
def read_attr_s(line) :
    attr_s=[]
    buf=[]
    is_quot=False
    for ch in line :
        if ch=='"' :
            is_quot=not is_quot
        if ch==';' and (not is_quot) :
            attr_s.append("".join(buf))
            buf=[]
            continue
        buf.append(ch)
    return attr_s

def get_mol_info(gtf_fn,feature=None) :
    for line in file(gtf_fn) :
        if line.startswith("#") :
            continue
        f=read_gtf_line(line,feature=feature)
        if f==None :
            continue
        yield f

def rebuild_chr_gene_set(path_finder,run) :
    if os.path.exists(path_finder.gtf_fn()) \
            and os.path.getsize(path_finder.gtf_fn())>30000 :
        sys.stdout.write("Found pre-build CHR_GTF_FN.\n")
        return True
    os.system("cp %s %s"%(run.gtf_fn,path_finder.gtf_fn()))
    return True

def rebuild_chr_fn(path_finder,run) :
    if os.path.exists(path_finder.chr_fn()) \
            and os.path.getsize(path_finder.chr_fn())>1000 :
        return True
    os.system("cp %s %s"%(run.genome_fa,path_finder.chr_fn()))
    return True

###
def rebuild_pep_fn(path_finder,run) :
    if os.path.exists(path_finder.pep_fn()) \
	    and os.path.getsize(path_finder.pep_fn())>1000 :
        return True
    os.system("cp %s %s"%(run.pep_fa,path_finder.pep_fn()))
    return True
###

def prepare_gene_id(path_finder):
    if os.path.exists(path_finder.gene_id_info()) :
        return True
    #
    coding_id_s=set()
    coding_mol_s={}
    #
    for mol in get_mol_info(path_finder.gtf_fn()) :
        if mol.gene_id not in  coding_id_s :
            f=Foo()
            f.tran_id_s=set()
            f.tran_id_s.add(mol.tran_id)
            coding_mol_s[mol.gene_id]=f
            coding_id_s.add(mol.gene_id)
        else :
            coding_mol_s[mol.gene_id].tran_id_s.add(mol.tran_id)
    #
    gene_id_s=coding_mol_s.keys()
    gene_id_s.sort()
    #
    out=file(path_finder.gene_id_info(),'wt')
    for gene_id in gene_id_s :
        gene=coding_mol_s[gene_id]
        tran_id_s=list(gene.tran_id_s)
        tran_id_s.sort()
        if len(tran_id_s)==0 :
            print "ERROR: No TRANSCRIPT_ID_s (%s)"%(gene_id)
            continue
        out.write(gene_id)
        out.write("\t")
        out.write("%d"%len(tran_id_s))
        out.write("\t")
        out.write(",".join(tran_id_s))
        out.write("\n")
    out.close()
    return True

def prepare_coding_gene_id(path_finder):
    if os.path.exists(path_finder.gene_id_coding_info()) :
        return True
    #
    sys.stdout.write("Starting to prepare CODING_GENE_ID.\n")
    #
    is_exon_only=False
    coding_id_s=set()
    for mol in get_mol_info(path_finder.gtf_fn(),"transcript") :
        if mol.transcript_biotype not in CODING_S :
            continue
        if mol.gene_id==None :
            sys.stderr.write("ERROR: Wrong LINE (%s)\n"%line)
            continue
        coding_id_s.add(mol.gene_id)
    if len(coding_id_s)< 10 :
        is_exon_only=True
        for mol in get_mol_info(path_finder.gtf_fn(),"exon") :
            if mol.transcript_biotype not in CODING_S :
                continue
            if mol.gene_id==None :
                sys.stderr.write("ERROR: Wrong LINE (%s)\n"%line)
                continue
            coding_id_s.add(mol.gene_id)

    coding_mol_s={}
    for gene_id in coding_id_s  :
        f=Foo()
        f.tran_id_s=set()
        f.tran_biotype_s=set()
        coding_mol_s[gene_id]=f
    #
    if not is_exon_only :
        for mol in get_mol_info(path_finder.gtf_fn(),"transcript") :
            if mol.gene_id not in coding_id_s :
                continue
            coding_mol_s[mol.gene_id].tran_id_s.add(mol.tran_id)
            coding_mol_s[mol.gene_id].tran_biotype_s.add(mol.transcript_biotype)
    else :
        for mol in get_mol_info(path_finder.gtf_fn(),"exon") :
            if mol.gene_id not in coding_id_s :
                continue
            coding_mol_s[mol.gene_id].tran_id_s.add(mol.tran_id)
            coding_mol_s[mol.gene_id].tran_biotype_s.add(mol.transcript_biotype)
    #
    gene_id_s=coding_mol_s.keys()
    gene_id_s.sort()
    #
    out=file(path_finder.gene_id_coding_info(),'wt')
    for gene_id in gene_id_s :
        gene=coding_mol_s[gene_id]
        tran_id_s=list(gene.tran_id_s)
        tran_id_s.sort()
        out.write(gene_id)
        out.write("\t")
        out.write("%d"%len(tran_id_s))
        out.write("\t")
        out.write(",".join(tran_id_s))
        out.write("\t")
        out.write(",".join(gene.tran_biotype_s))
        out.write("\n")
    out.close()
    return True

def build_gene_set(path_finder) :
    if os.path.exists(path_finder.info_geneset()) \
            and os.path.getsize(path_finder.info_geneset())>100 :
        sys.stdout.write("Found pre-built GENESET_S.\n")
        return True

    meta=Meta(path_finder)
    coding_id_s=meta.get_coding_id_s()
    sys.stdout.write("found %d protein_coding_genes.\n"%len(coding_id_s))
    #
    mask_s=[]
    noncoding_s=[]
    mt_s=[]
    pt_s=[]
    #
    out_coding=file(path_finder.coding_gtf(),'wt')
    out_etc=file(path_finder.t_r_mt_pt_rna_gtf(),'wt')
    out_mask=file(path_finder.mask_gtf(),'wt')
    for mol in get_mol_info(path_finder.gtf_fn()) :
        if mol.chr in MT_name_s :
            mt_s.append(mol.line)
            out_mask.write(mol.line)
        #
        if mol.chr in PT_name_s :
            pt_s.append(mol.line)
            out_mask.write(mol.line)
        #
        if mol.gene_id in coding_id_s :
            out_coding.write(mol.line)
        elif (mol.chr in MT_name_s) or (mol.chr in PT_name_s) :
            noncoding_s.append(mol.line)
        else :
            noncoding_s.append(mol.line)
            if mol.chr in MT_name_s :
                pass
            elif mol.chr in PT_name_s :
                pass
            else :
                out_mask.write(mol.line)
        #
        if mol.chr in MT_name_s \
                or mol.chr in PT_name_s \
                or mol.gene_biotype in T_R_RNA_S :
            out_etc.write(mol.line)

    out_mask.close()
    out_coding.close()
    out_etc.close()
    # Non-coding : All
    out=file(path_finder.non_coding_all_gtf(),'wt')
    for line in noncoding_s :
        out.write(line)
    out.close()
    # tRNA & rRNA
    out=file(path_finder.t_r_rna_gtf(),'wt')
    for line in noncoding_s :
        unit_s=line.strip().split("\t")
        if unit_s[1] not in T_R_RNA_S :
            continue
        out.write(line)
    out.close()
    # mtRNA
    out=file(path_finder.mt_rna_gtf(),'wt')
    for line in mt_s :
        out.write(line)
    out.close()
    # ptRNA
    out=file(path_finder.pt_rna_gtf(),'wt')
    for line in pt_s :
        out.write(line)
    out.close()
    # lincRNA
    out=file(path_finder.linc_gtf(),'wt')
    for line in noncoding_s :
        unit_s=line.strip().split("\t")
        if unit_s[1] not in LINC_RNA_S :
            continue
        out.write(line)
    out.close()
    # miRNA
    out=file(path_finder.mirna_gtf(),'wt')
    for line in noncoding_s :
        unit_s=line.strip().split("\t")
        if unit_s[1] not in MI_RNA_S :
            continue
        out.write(line)
    out.close()
    #
    to_rel_path=path_finder.to_rel_path
    #
    out=file(path_finder.info_geneset(),'wt')
    out.write("CODING\t%s\n"%(",".join(CODING_S)))
    out.write("NCRNA\t%s\n"%(','.join(NONCODING_S)))
    out.write("T_R_RNA\t%s\n"%(",".join(T_R_RNA_S)))
    out.write("LINC_RNA\t%s\n"%(",".join(LINC_RNA_S)))
    out.write("MI_RNA\t%s\n"%(",".join(MI_RNA_S)))
    out.write("#\n")
    out.write("GENESET\tSOURCE:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn()))
    out.write("GENESET\tALL:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn()))
    out.write("GENESET\tNONCODING:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path((path_finder.non_coding_all_gtf())))
    out.write("GENESET\tCODING:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.coding_gtf()))
    out.write("GENESET\tMASK:MT:PT:FULL:NONCODING:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.mask_gtf()))
    out.write("#\n")
    out.close()
    return True

def build_exonic_gene_set(path_finder) :
    info_fn=path_finder.info_geneset_exonic()
    if os.path.exists(info_fn) \
            and os.path.getsize(info_fn)>100 :
        sys.stdout.write("Found pre-built EXONIC_GENESET_S.\n")
        return True

    gtf_fn_s=(\
            (path_finder.gtf_fn(),path_finder.gtf_fn(tag='exonic')) ,\
            (path_finder.coding_gtf(),path_finder.coding_gtf(tag='exonic')) ,\
            (path_finder.t_r_mt_pt_rna_gtf(),path_finder.t_r_mt_pt_rna_gtf(tag='exonic')), \
            (path_finder.mask_gtf(),path_finder.mask_gtf(tag='exonic')), \
            (path_finder.non_coding_all_gtf(),path_finder.non_coding_all_gtf(tag='exonic')), \
            (path_finder.t_r_rna_gtf(),path_finder.t_r_rna_gtf(tag='exonic')),\
            (path_finder.mt_rna_gtf(),path_finder.mt_rna_gtf(tag='exonic')), \
            (path_finder.pt_rna_gtf(),path_finder.pt_rna_gtf(tag='exonic')), \
            (path_finder.linc_gtf(),path_finder.linc_gtf(tag='exonic')), \
            (path_finder.mirna_gtf(),path_finder.mirna_gtf(tag='exonic')), \
            )

    for gtf_fn,out_gff in gtf_fn_s :
        buf=[]
        buf.append("dexseq_prepare_annotation.py")
        buf.append(gtf_fn)
        buf.append(out_gff)
        os.system(" ".join(buf))

    to_rel_path=path_finder.to_rel_path

    out=file(path_finder.info_geneset_exonic(),'wt')
    out.write("GENESET\tALL:EXONIC:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn(tag='exonic')))
    out.write("GENESET\tNONCODING:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path((path_finder.non_coding_all_gtf(tag='exonic'))))
    out.write("GENESET\tCODING:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.coding_gtf(tag='exonic')))
    out.write("GENESET\tMASK:MT:PT:NONCODING:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.mask_gtf(tag='exonic')))
    out.close()

    return True

def prepare_repr_gene_set(path_finder,info) :
    repr_id_s=read_repr_id_gene(info)
    sys.stdout.write("Found %d repr transcripts.\n"%len(repr_id_s))

    out=file(path_finder.repr_gtf(),'wt')
    for line in file(path_finder.gtf_fn()) :
        unit_s=line.split("\t")
        tran_id=None
        for attr in unit_s[8].split(";") :
            attr=attr.strip()
            if attr=='' :
                continue
            key,val=attr.split(None,1)
            if key=='transcript_id' :
                tran_id=val[1:-1]
                break
        if tran_id==None :
            continue
        if tran_id not in repr_id_s :
            continue
        out.write(line)
    out.close()

    out=file(path_finder.coding_repr_gtf(),'wt')
    for line in file(path_finder.coding_gtf()) :
        unit_s=line.split("\t")
        tran_id=None
        for attr in unit_s[8].split(";") :
            attr=attr.strip()
            if attr=='' :
                continue
            key,val=attr.split(None,1)
            if key=='transcript_id' :
                tran_id=val
                break
        if tran_id==None :
            continue
        if tran_id.startswith('"') and tran_id.endswith('"') :
            tran_id=tran_id[1:-1]
        if tran_id not in  repr_id_s :
            continue
        out.write(line)
    out.close()
    
    out=file(path_finder.non_coding_repr_gtf(),'wt')
    for line in file(path_finder.non_coding_all_gtf()) :
        unit_s=line.split("\t")
        tran_id=None
        for attr in unit_s[8].split(";") :
            attr=attr.strip()
            if attr=='' :
                continue
            key,val=attr.split(None,1)
            if key=='transcript_id' :
                tran_id=val
                break
        if tran_id==None :
            continue
        if tran_id.startswith('"') and tran_id.endswith('"') :
            tran_id=tran_id[1:-1]
        if tran_id not in  repr_id_s :
            continue
        out.write(line)
    out.close()

# ------------------------------------------------------------------------------
# Build index
#
def build_gene(path_finder,db) :
    source_fn='%s/gene.txt.gz'%(path_finder.db_path())
    done_fn='%s.done'%source_fn
    if os.path.exists(done_fn) :
        return True
    gene_s=[]
    for line in os.popen("zcat %s"%source_fn) :
        gene=line.strip().split("\t")[:-2]
        val_s=[]
        for val in gene :
            if val=='\N' :
                val_s.append(None)
            else :
                val_s.append(val)
        gene_s.append(val_s)
    db.insert_gene_s(gene_s)
    db.commit()
    os.system("touch %s"%done_fn)

def build_transcript(path_finder,db) :
    source_fn='%s/transcript.txt.gz'%(path_finder.db_path())
    done_fn='%s.done'%source_fn
    if os.path.exists(done_fn) :
        return True
    transc_s=[]
    for line in os.popen("zcat %s"%(source_fn)) :
        transc=line.strip().split("\t")[:-2]
        val_s=[]
        for val in transc :
            if val=='\N' :
                val_s.append(None)
            else :
                val_s.append(val)
        transc_s.append(val_s)
    db.insert_transc_s(transc_s)
    db.commit()
    os.system("touch %s"%done_fn)

def build_translation(path_finder,db) :
    source_fn='%s/translation.txt.gz'%(path_finder.db_path())
    done_fn='%s.done'%source_fn
    if os.path.exists(done_fn) :
        return True
    transl_s=[]
    for line in os.popen("zcat %s"%(source_fn)) :
        transl=line.strip().split("\t")[:-2]
        val_s=[]
        for val in transl :
            if val=='\N' :
                val_s.append(None)
            else :
                val_s.append(val)
        transl_s.append(val_s)
    db.insert_transl_s(transl_s)
    db.commit()
    os.system("touch %s"%done_fn)

def build_exon(path_finder,db) :
    source_fn='%s/exon.txt.gz'%(path_finder.db_path())
    done_fn='%s.done'%source_fn
    if os.path.exists(done_fn) :
        return True
    in_s=[]
    for line in os.popen("zcat %s"%source_fn) :
        _in=line.strip().split("\t")[:-2]
        val_s=[]
        for val in _in :
            if val=='\N' :
                val_s.append(None)
            else :
                val_s.append(val)
        in_s.append(val_s)
    db.insert_exon_s(in_s)
    db.commit()
    os.system("touch %s"%done_fn)

def build_protein_feature(path_finder,db) :
    key="protein_feature"
    sys.stdout.write("INFO: Reading %s.\n"%key)
    n_unit_s=None
    in_s=[]
    source_fn='%s/%s.txt.gz'%(path_finder.db_path(),key)
    done_fn="%s.done"%source_fn
    if os.path.exists(done_fn) :
        return True
    if not os.path.exists(source_fn) :
        return False
    for line in os.popen('zcat %s'%source_fn) :
        _in=line.replace("\n","").split("\t")
        val_s=[]
        for val in _in :
            if val=='\N' :
                val_s.append(None)
            else :
                val_s.append(unicode(val, errors='ignore'))
        if n_unit_s==None :
            n_unit_s=len(val_s)
        if n_unit_s != len(val_s) :
            sys.stderr.write("ERROR: Failed to parse the next lines from %s\n"%key)
            sys.stderr.write("ERROR: %d <-> %d: %s\n"%(n_unit_s,len(val_s),line.strip()))
        else :
            in_s.append(val_s[:12])
    sys.stdout.write("INFO: Inserting %s.\n"%key)
    db.insert_in_s(key,in_s)
    db.commit()
    os.system("touch %s"%done_fn)
    return True

def build_extra(path_finder,db) :
    for key in DB_EXTRA_s :
        sys.stdout.write("INFO: Reading %s.\n"%key)
        n_unit_s=None
        in_s=[]
        source_fn='%s/%s.txt.gz'%(path_finder.db_path(),key)
        done_fn="%s.done"%source_fn
        if os.path.exists(done_fn) :
            continue
        for line in os.popen('zcat %s'%source_fn) :
            _in=line.replace("\n","").split("\t")
            val_s=[]
            for val in _in :
                if val=='\N' :
                    val_s.append(None)
                else :
                    val_s.append(unicode(val, errors='ignore'))
            if n_unit_s==None :
                n_unit_s=len(val_s)
            if n_unit_s != len(val_s) :
                sys.stderr.write("ERROR: Failed to parse the next lines from %s\n"%key)
                sys.stderr.write("ERROR: %d <-> %d: %s\n"%(n_unit_s,len(val_s),line.strip()))
            else :
                in_s.append(val_s)
#                print line
#                db.insert_in(key,val_s)
        sys.stdout.write("INFO: Inserting %s.\n"%key)
        db.insert_in_s(key,in_s)
        db.commit()
        os.system("touch %s"%done_fn)

def gen_gene2tran(path_finder,gene_s) :
    fn=path_finder.gene_to_tran()
    sys.stdout.write("INFO\tGenerating GENE2TRAN (%s)\n"%fn)
    out_map=file(fn,'wt')
    out_map.write("#GeneID\tTranID\tTranLength\tGeneType\tTranType\n")
    for gene in gene_s :
        for tran in gene.tran_s :
            out_map.write("%s"%(gene.gene_id))
            if tran.trans_name==None :
                out_map.write("\t-")
            else :
                out_map.write("\t%s"%(tran.trans_name))
            out_map.write("\t%d"%len(tran))
            if tran.transl_name==None :
                out_map.write("\t-")
            else :
                out_map.write("\t%s"%(tran.transl_name))
            out_map.write("\t%s"%(gene.gene_biotype))
            out_map.write("\t%s"%(tran.transcript_biotype))
            out_map.write("\n")
    out_map.close()

def gen_mol_goa(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating MOL_GOA\n")

    out_gene=file("%s/%s"%(path_finder.home(),ENS_GOA_GENE_FN),'wt')
    for gene in gene_s :
        out_gene.write(gene.gene_id)
        if len(gene.go_id_s)==0 :
            out_gene.write("\tNA\n")
        else :
            out_gene.write("\t%s\n"%(",".join(gene.go_id_s)))
    out_gene.close()

    out_tran=file("%s/%s"%(path_finder.home(),ENS_GOA_TRANSCRIPT_FN),'wt')
    for gene in gene_s :
        for tran in gene.tran_s :
            if tran.transl_name==None :
                continue
            #
            out_tran.write("%s"%(tran.trans_name))
            if len(tran.go_id_s)==0 :
                out_tran.write("\tNA\n")
            else :
                out_tran.write("\t%s\n"%(",".join(tran.go_id_s)))
    out_tran.close()

def gen_gene_desc(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating GENE_DESC\n")
    out=file("%s/%s"%(path_finder.home(),ENS_GENE_DESC_FN),'wt')
    for gene in gene_s :
        out.write("%s"%(gene.gene_id))
        if gene.desc==None :
            out.write("\t-")
        else :
            out.write("\t%s"%(gene.desc))
        out.write("\n")
    out.close()

def gen_gene_info(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating GENE_INFO\n")
    out=file(path_finder.gene_info(),'wt')
    out.write("#StableID\tGeneType\tGeneName\tGeneStatus")
    out.write("\tReprTranID\tReprTranType\tTranName\tReprTranLen")
    out.write("\tDESC\tGO\tInterPro\tEntrezGene\tRefSeq_mRNA\tGO_with_Anch")
    out.write("\tGeneChr\tGeneStrand\tGeneStart\tGeneEnd")
    out.write("\tTranChr\tTranStrand\tTranStart\tTranEnd")
    out.write("\n")
    for gene in gene_s :
        repr=gene.repr_tran()
        if repr==None :
            print "No REPR_TRAN (%s)"%(gene.id())
        #
        out.write("%s"%(gene.gene_id))
        out.write("\t%s"%(gene.gene_biotype))
        #
        if gene.symbol==None :
            out.write("\t-")
        else :
            out.write("\t%s"%(gene.symbol))
        #
        out.write("\t%s"%(gene.status))
        #
        out.write("\t%s\t%s\t%s\t%d"%\
                (repr.id(),repr.transcript_biotype,repr.symbol,len(repr)))
        #
        if gene.desc==None :
            out.write("\t-")
        else :
            out.write("\t%s"%(gene.desc))
        #
        if len(gene.go_id_s)==0 :
            out.write("\tNA")
        else :
            out.write("\t%s"%(",".join(gene.go_id_s)))
        #
        if len(gene.interpro_s)==0 :
            out.write("\tNA")
        else :
            interpro_s=list(gene.interpro_s)
            interpro_s.sort()
            out.write("\t%s"%(",".join(interpro_s)))
        #
        if gene.ncbi_gene_s==None :
            out.write("\tNA")
        else :
            ncbi_gene=gene.ncbi_gene_s[0]
            out.write("\t%s"%ncbi_gene.dbprimary_acc)
        #
        if repr.refseq_s==None :
            out.write("\tNA")
        else :
            refseq=repr.refseq_s[0]
            out.write("\t%s"%refseq.dbprimary_acc)
        #
        if gene.go_with_anch_id_s==None or  len(gene.go_with_anch_id_s)==0 :
            out.write("\tNA")
        else :
            out.write("\t%s"%(",".join(gene.go_with_anch_id_s)))
        #
        out.write("\t%s"%(gene.chr()))
        out.write("\t%s"%(gene.strand()))
        out.write("\t%d"%(gene.start()))
        out.write("\t%d"%(gene.end()))
        #
        out.write("\t%s"%(repr.chr()))
        out.write("\t%s"%(repr.strand()))
        out.write("\t%d"%(repr.start()))
        out.write("\t%d"%(repr.end()))
        #
        out.write("\n")
    out.close()

def gen_gene_features(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating GENE_FEATURES\n")

    out=file(path_finder.gene_feature(),'wt')
    #
    out.write("#StableID\tGeneType\tGeneName\tGeneStatus")
    for FEATURE in PROTEIN_FEATURE_S :
        out.write("\t%s"%FEATURE)
    out.write("\n")
    #
    for gene in gene_s :
        out.write("%s"%(gene.gene_id))
        out.write("\t%s"%(gene.biotype))
        #
        if gene.symbol==None :
            out.write("\t-")
        else :
            out.write("\t%s"%(gene.symbol))
        #
        out.write("\t%s"%(gene.status))
        #
        if len(gene.interpro_s)==0 :
            out.write("\tNA")
        else :
            interpro_s=list(gene.interpro_s)
            interpro_s.sort()
            out.write("\t%s"%(",".join(interpro_s)))
        #
        for FEATURE in PROTEIN_FEATURE_S :
            f_id_s=set()
            for tran in gene.tran_s :
                if tran.transl_id==None :
                    continue
                for feature in tran.feature_s[FEATURE] :
                    f_id_s.add(feature.hit_name)
            #
            if len(f_id_s)==0 :
                out.write("\tNA")
            else :
                f_id_s=list(f_id_s)
                f_id_s.sort()
                out.write("\t%s"%(",".join(f_id_s)))
        out.write("\n")
    out.close()

def gen_tran_info(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating TRAN_INFO\n")
    out=file(path_finder.tran_info(),'wt')
    out.write("#StableID\tTranType\tTranName\tTranLength")
    out.write("\tGeneID\tGeneType\tGeneName\tGeneStatus")
    out.write("\tDESC\tGO\tInterPro\tRefSeq_mRNA")
    out.write("\tEntrezGene\tProteinID")
    out.write("\tGeneChr\tGeneStrand\tGeneStart\tGeneEnd")
    out.write("\tTranChr\tTranStrand\tTranStart\tTranEnd")
    out.write("\n")
    for gene in gene_s :
        for tran in gene.tran_s :
            out.write("%s"%(tran.id()))
            out.write("\t%s"%(tran.transcript_biotype))
            #
            if gene.symbol==None :
                out.write("\t-")
            else :
                out.write("\t%s"%(gene.symbol))
            #
            out.write("\t%d"%len(tran))
            #
            out.write("\t%s\t%s"%(gene.gene_id,gene.gene_biotype))
            if gene.symbol==None :
                out.write("\t-")
            else :
                out.write("\t%s"%(gene.symbol))
            out.write("\t%s"%(gene.status))
            #
            if gene.desc==None :
                out.write("\t-")
            else :
                out.write("\t%s"%(gene.desc))
            #
            if tran.go_id_s==None or len(tran.go_id_s)==0 :
                out.write("\tNA")
            else :
                out.write("\t%s"%(",".join(tran.go_id_s)))
            #
            if tran.interpro_s==None or len(tran.interpro_s)==0 :
                out.write("\tNA")
            else :
                interpro_s=list(tran.interpro_s)
                interpro_s.sort()
                out.write("\t%s"%(",".join(interpro_s)))
            #
            if tran.refseq_s==None :
                out.write("\tNA")
            else :
                refseq=tran.refseq_s[0]
                out.write("\t%s"%refseq.dbprimary_acc)
            #
            if gene.ncbi_gene_s==None :
                out.write("\tNA")
            else :
                ncbi_gene=gene.ncbi_gene_s[0]
                out.write("\t%s"%ncbi_gene.dbprimary_acc)
            #
            if tran.protein==None :
                out.write("\t-")
            else :
                out.write("\t%s"%(tran.protein.stable_id))
            #
            out.write("\t%s"%(gene.chr()))
            out.write("\t%s"%(gene.strand()))
            out.write("\t%d"%(gene.start()))
            out.write("\t%d"%(gene.end()))
            #
            out.write("\t%s"%(tran.chr()))
            out.write("\t%s"%(tran.strand()))
            out.write("\t%d"%(tran.start()))
            out.write("\t%d"%(tran.end()))
            #
            out.write("\n")
    out.close()

def gen_tran_features(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating TRAN_FEATURES\n")

    out=file(path_finder.tran_feature(),'wt')
    #
    out.write("#StableID\tTranType\tTranName\tTranLength")
    for FEATURE in PROTEIN_FEATURE_S :
        out.write("\t%s"%FEATURE)
    out.write("\n")
    #
    for gene in gene_s :
        for tran in gene.tran_s :
            out.write("%s"%(tran.trans_name))
            out.write("\t%s"%(tran.biotype))
            #
            if gene.symbol==None :
                out.write("\t-")
            else :
                out.write("\t%s"%(gene.symbol))
            #
            out.write("\t%d"%(tran._len))
            #
            if len(tran.interpro_s)==0 :
                out.write("\tNA")
            else :
                interpro_s=list(tran.interpro_s)
                interpro_s.sort()
                out.write("\t%s"%(",".join(interpro_s)))
            #
            for FEATURE in PROTEIN_FEATURE_S :
                if tran.transl_id==None :
                    out.write("\tNA")
                    continue
                f_id_s=set()
                for feature in tran.feature_s[FEATURE] :
                    f_id_s.add(feature.hit_name)
                f_id_s=list(f_id_s)
                f_id_s.sort()
                out.write("\t%s"%(",".join(f_id_s)))
            #
            if tran.refseq_s==None :
                out.write("\tNA")
            else :
                refseq=tran.refseq_s[0]
                out.write("\t%s"%refseq.dbprimary_acc)

            out.write("\n")
    out.close()


def gen_gene_tran_map(path_finder,gene_s) :
    sys.stdout.write("INFO\tGenerating GEN_TRAN_MAP\n")
    out=file(path_finder.gene_tran_map(),'wt')
    for gene in gene_s :
        for tran in gene.tran_s :
            out.write("%s"%(gene.gene_id))
            out.write("\t")
            out.write("%s"%(tran.trans_name))
            out.write("\n")
    out.close()

def build_db(path_finder) :
    db=RineDB(path_finder.db_fn())
    db.init()
    build_gene(path_finder,db)
    build_transcript(path_finder,db)
    build_translation(path_finder,db)
    build_exon(path_finder,db)
    if not build_protein_feature(path_finder,db) :
        sys.stderr.write("ERROR: Failed to build protein_features.\n")
        return
    build_extra(path_finder,db)
    os.system("touch %s"%(path_finder.db_done()))
    return True

def sort_by_name(one,other) :
    if one.gene_id > other.gene_id :
        return 1
    if one.gene_id < other.gene_id :
        return -1
    return 0

def build_desc(path_finder,species) :
    gene_s=[]
    for gene in read_gtf_gene(path_finder.gtf_fn()) :
        gene_s.append(gene)
    gene_s.sort()
    #
    for gene in gene_s :
        for tran in gene.tran_s :
            tran.refseq_s=None
        #
        #repr=gene.repr_tran()
    #
    gen_gene2tran(path_finder,gene_s)
    gen_gene_desc(path_finder,gene_s)
    #
    gen_mol_goa(path_finder,gene_s)
    gen_gene_info(path_finder,gene_s)
    gen_tran_info(path_finder,gene_s)
    #
    #gen_gene_features(path_finder,gene_s)
    #gen_tran_features(path_finder,gene_s)
    #gen_gene_tran_map(path_finder,gene_s)
    #
    to_rel_path=path_finder.to_rel_path
    #
    out=file(path_finder.info_desc_fn(),'wt')
    out.write("#ENS_DB\t%s\n"%(path_finder.db_fn()))
    out.write("#ENS_GENE_TRANSCRIPT\t%s\n"%to_rel_path(path_finder.gene_to_tran()))
    out.write("#ENS_GENE_DESC\t%s/%s\n"%(path_finder.home(),ENS_GENE_DESC_FN))
    out.write("#ENS_GO_GENE\t%s/%s\n"%(path_finder.home(),ENS_GOA_GENE_FN))
    out.write("#ENS_GO_TRANSCRIPT\t%s/%s\n"%(path_finder.home(),ENS_GOA_TRANSCRIPT_FN))
    #out.write("ENS_GENE_INFO\t%s\n"%(path_finder.gene_info()))
    #out.write("ENS_GENE_FEATURE\t%s\n"%(path_finder.gene_feature()))
    #out.write("ENS_TRAN_INFO\t%s\n"%(path_finder.tran_info()))
    #out.write("ENS_TRAN_FEATURE\t%s\n"%(path_finder.tran_feature()))
    out.write("INFO\tGENES_TRANS\t%s\n"%to_rel_path(path_finder.gene_to_tran()))
    out.write("GENE\tINFO\t%s\n"%to_rel_path(path_finder.gene_info()))
    out.write("GENE\tFEATURE\t%s\n"%to_rel_path(path_finder.gene_feature()))
    out.write("TRANSCRIPT\tINFO\t%s\n"%to_rel_path(path_finder.tran_info()))
    out.write("TRANSCRIPT\tFEATURE\t%s\n"%to_rel_path(path_finder.tran_feature()))
    out.write("#\n")
    out.write("FEATURE\tGENE_TRAN:MAP\t%s\n"%to_rel_path(path_finder.gene_tran_map()))

    out.close()


def build_mapping(path_finder,species) :
    if os.path.exists(path_finder.gene_tran_map()) \
            and os.path.getsize(path_finder.gene_tran_map())>1000 :
        sys.stdout.write(\
                "INFO: Found GENE_TRAN_MAP (%s).\n"%\
                (path_finder.gene_tran_map()))
        return True

    tran_id_s=set()
    out=file(path_finder.gene_tran_map(),'wt')
    for mol in get_mol_info(path_finder.gtf_fn()) :
        if mol.tran_id in tran_id_s :
            continue
        tran_id_s.add(mol.tran_id)
        out.write("%s\t%s\n"%(mol.gene_id,mol.tran_id))
    out.close()

    return True


def build_repeat(path_finder,species) :
    sys.stdout.write("Starting to generate repeat bed file.\n")
    #
    info_dn=read_info_down(path_finder.info_down_fn())
    list_fn="%s.list"%(info_dn.dna_fn)
    dna_s=read_dna_s(list_fn)
    chr_name_s=set([dna.name for dna in dna_s.ref_dna_s()])
    #
    db=RineDB(path_finder.db_fn())
    #
    repeat_s=db.search_repeat()
    sys.stdout.write("REPEAT (%s)\n"%len(repeat_s))
    sys.stdout.write(",".join(chr_name_s))
    #
    i_repeat=0
    out=file(path_finder.repeat_bed_fn(),'wt')
    #out.write("Chr\tStart\tEnd\tName\tClass\tType\tConsensus\n")
    for repeat in repeat_s :
        if repeat.chr not in chr_name_s :
            continue
        out.write("%s\t%s\t%s"%(repeat.chr,repeat.start-1,repeat.end))
        out.write("\t%s"%(repeat.name))
        out.write("\t%s"%(repeat._class))
        out.write("\t%s"%(repeat._type))
        out.write("\t%s"%(repeat.consensus))
        out.write("\n")
        i_repeat+=1
    out.close()
    #
    out=file(path_finder.info_repeat_fn(),'wt')
    out.write("%d\n"%(i_repeat))
    out.close()


def build_bowtie_idx(fa_fn) :
    sys.stdout.write("INFO\tGenerating BOWTIE_IDX for %s\n"%fa_fn)
    short_name=fa_fn
    if short_name.endswith(".fa") :
        short_name=short_name[:-3]
    elif short_name.endswith(".fasta") :
        short_name=short_name[:-6]

    if not os.path.exists("%s.1.bt2"%short_name) \
            and not os.path.exists("%s.1.bt2"%short_name):
        sys.stdout.write("INFO: Starting to build BOWTIE2_IDX")
        os.system("%s %s %s"%(BOWTIE2_BUILD_EXEC,fa_fn,short_name))

    if not os.path.exists("%s.1.ebwt"%short_name) \
            and not os.path.exists("%s.1.ebwt"%short_name):
        sys.stdout.write("INFO: Starting to build BOWTIE1_IDX")
        os.system("%s %s %s"%(BOWTIE1_BUILD_EXEC,fa_fn,short_name))

    if not os.path.exists("%s.1.ebwt"%short_name) \
            or not os.path.exists("%s.1.bt2"%short_name) :
        msg="ERRROR\tFailed to build bowtie index (%s.1.ebwt).\n"%short_name
        sys.stderr.write(msg)
        return 

    return short_name


def build_bwa_idx(fa_fn) :
    bwt_fn="%s.bwt"%fa_fn
    if not os.path.exists(bwt_fn) :
        sys.stdout.write("INFO\tStarting to build index for bwa (%s)\n"%bwt_fn)
        os.system("%s index -a bwtsw %s"%(BWA_EXEC,fa_fn))
    if not os.path.exists("%s.bwt"%fa_fn) :
        msg="ERRROR\tFailed to build BWA INDEX (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return fa_fn

def build_star_idx(path_finder,fa_fn) :
    if not os.path.exists(path_finder.staridx_path()) :
        os.makedirs(path_finder.staridx_path())
    idx_fn="%s/SAindex"%path_finder.staridx_path()

    cmd=[]
    cmd.append(STAR_EXEC)
    cmd.append("--runMode")
    cmd.append("genomeGenerate")
    cmd.append("--genomeDir")
    cmd.append(path_finder.staridx_path())
    cmd.append("--genomeFastaFiles")
    cmd.append(fa_fn)
    cmd.append("--limitGenomeGenerateRAM")
    cmd.append(STAR_LIMITRAM)
    cmd.append("--runThreadN")
    cmd.append("%d"%STAR_N_CPU)
    if not os.path.exists(idx_fn) :
        sys.stdout.write("INFO\tStarting to build index for STAR\n")
        print " ".join(cmd)
        os.system(" ".join(cmd))
    if not os.path.exists(idx_fn) :
        msg="ERRROR\tFailed to build STAR INDEX (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return path_finder.staridx_path()

def build_fai_idx(fa_fn) :
    if not os.path.exists("%s.fai"%fa_fn) :
        sys.stdout.write("INFO\tStarting to build fasta index\n")
        os.system("%s faidx %s"%(SAMTOOLS_EXEC,fa_fn))
    if not os.path.exists("%s.fai"%fa_fn) :
        msg="ERRROR\tFailed to build fasta index (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return "%s.fai"%fa_fn

def build_blast_n_idx(fa_fn) :
    if not os.path.exists("%s.nin"%fa_fn) :
        sys.stdout.write("INFO\tStarting to build BLAST index\n")
        os.system("%s -in %s -dbtype nucl"%(BLAST_EXEC,fa_fn))
    if not os.path.exists("%s.nin"%fa_fn) :
        msg="ERRROR\tFailed to build BLAST index (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return fa_fn

def build_blast_p_idx(fa_fn) :
    if not os.path.exists("%s.pin"%fa_fn) :
        sys.stdout.write("INFO\tStarting to build BLAST index\n")
        os.system("%s -in %s -dbtype prot"%(BLAST_EXEC,fa_fn))
    if not os.path.exists("%s.pin"%fa_fn) :
        msg="ERRROR\tFailed to build BLAST index (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return fa_fn

def build_dict_idx(fa_fn) :
    dict_fn=None
    if fa_fn.endswith(".fa") :
        dict_fn="%s.dict"%(fa_fn[:-3])
    elif fa_fn.endswith(".fasta") :
        dict_fn="%s.dict"%(fa_fn[:-6])
    else :
        dict_fn="%s.dict"%(fa_fn)
    #
    if not os.path.exists(dict_fn) :
        sys.stdout.write("INFO\tStarting to build dict index. (%s)\n"%dict_fn)
        os.system("java -jar %s R=%s O=%s"%\
                (PICARD_CREATESEQUENCEDICT_EXEC,fa_fn,dict_fn))
    if not os.path.exists(dict_fn) or os.path.getsize(dict_fn)<100 :
        msg="ERRROR\tFailed to build fasta index (%s).\n"%fa_fn
        sys.stderr.write(msg)
        return 
    return "%s.fai"%fa_fn


def read_info_down(info_fn) :
    dn=Downloads()
    short_name=None
    for line in file(info_fn) :
        if line.startswith("#") :
            continue
        unit_s=line.strip().split("\t")
        tag=unit_s[0]
        if tag=="SHORT_NAME" :
            dn.short_name=line.strip().split()[1]
        elif tag=='SPECIES' :
            dn.species=line.strip().split()[1]
        elif tag=='RELEASE' :
            dn.release=line.strip().split()[1]
        elif tag=='GTF' :
            dn.ens_gtf_fn=unit_s[1]
        elif tag=='PEP' :
            dn.pep_fn=unit_s[1]
        elif line.startswith("DNA\t") :
            unit_s=line.strip().split()
            if unit_s[1].endswith("toplevel.fa.gz") :
                dn.dna_fn=unit_s[2]
        elif line.startswith("CDNA\t") :
            unit_s=line.strip().split()
            dn.cdna_fn=unit_s[1]
        elif line.startswith("NCRNA\t") :
            unit_s=line.strip().split()
            dn.ncrna_fn=unit_s[1]
    dn.uncompress()
    return dn

def chr_name_split(chr_name_s) :
    major_s=[]
    minor_s=[]

    if "1" in chr_name_s and "2" in chr_name_s :
        for name in chr_name_s :
            is_int=False
            try :
                i_name=int(name)
                is_int=True
            except ValueError:
                pass
            if is_int :
                major_s.append(name)
            elif name in SELECT_S :
                major_s.append(name)
            else :
                minor_s.append(name)
    else :
        for name in chr_name_s :
            major_s.append(name)
            minor_s.append(name)
    return major_s,minor_s

def read_dna_s(list_fn) :
    dna_s=DnaGroup()
    for line in file(list_fn) :
        if line.startswith("#") :
            continue
        unit_s=line.strip().split()
        dna=DNA(unit_s[0],level=unit_s[3].split(":")[1])
        dna.fn=unit_s[1]
        sub_s=unit_s[4].split(":")
        dna.source=sub_s[1]
        dna.length=int(sub_s[4])
        if len(unit_s)>5 :
            dna.type=unit_s[5]
        dna_s.append(dna)
    dna_s.sort()
    return dna_s

def build_index(path_finder) :
    if os.path.exists(path_finder.info_idx_fn()) :
        sys.stdout.write("INFO: Found INFO.IDX(%s).\n"%\
                (path_finder.info_idx_fn()))
        return True

    sys.stdout.write("INFO BUILD_INDEX\n")
    #
    chr_bowtie_fn=build_bowtie_idx(path_finder.chr_fn())
    chr_bwa_fn=build_bwa_idx(path_finder.chr_fn())
    chr_fai_fn=build_fai_idx(path_finder.chr_fn())
    chr_dict_fn=build_dict_idx(path_finder.chr_fn())
    chr_star_fn=build_star_idx(path_finder,path_finder.chr_fn())
    chr_blast_fn=build_blast_n_idx(path_finder.chr_fn())
    #
    rna_bowtie_fn=build_bowtie_idx(path_finder.rna_fn())
    rna_bwa_fn=build_bwa_idx(path_finder.rna_fn())
    rna_fai_fn=build_fai_idx(path_finder.rna_fn())
    rna_dict_fn=build_dict_idx(path_finder.rna_fn())
    rna_blast_fn=build_blast_n_idx(path_finder.rna_fn())
    #
    rna_coding_bowtie_fn=build_bowtie_idx(path_finder.rna_coding_fn())
    rna_coding_bwa_fn=build_bwa_idx(path_finder.rna_coding_fn())
    rna_coding_fai_fn=build_fai_idx(path_finder.rna_coding_fn())
    rna_coding_dict_fn=build_dict_idx(path_finder.rna_coding_fn())
    rna_coding_blast_fn=build_blast_n_idx(path_finder.rna_coding_fn())
    #
    pep_blast_fn=build_blast_p_idx(path_finder.pep_fn())
    if os.path.exists(path_finder.ncrna_fn()) \
            and os.path.getsize(path_finder.ncrna_fn())>0 :
        ncrna_bowtie_fn=build_bowtie_idx(path_finder.ncrna_fn())
        ncrna_bwa_fn=build_bwa_idx(path_finder.ncrna_fn())
        ncrna_fai_fn=build_fai_idx(path_finder.ncrna_fn())
        ncrna_dict_fn=build_dict_idx(path_finder.rna_coding_fn())
    #
    if rna_bowtie_fn==None  :
        sys.stderr.write("Failed to build RNA_BOWTIE.\n")
        return False
    if rna_bwa_fn==None :
        sys.stderr.write("Failed to build RNA_BWA.\n")
        return False
    #
    if chr_bowtie_fn==None :
        sys.stderr.write("Failed to build CHR_BOWTIE.\n")
        return False
    if chr_bwa_fn==None :
        sys.stderr.write("Failed to build CHR_BWA.\n")
        return False
    if chr_star_fn==None :
        sys.stderr.write("Failed to build CHR_STAR.\n")
        return False
    #
    if chr_fai_fn==None :
        sys.stderr.write("Failed to build CHR_FAI_FN (%s).\n"%\
                (path_finder.chr_fn()))
        return False
    if rna_fai_fn==None :
        sys.stderr.write("Failed to build RNA_FAI_FN (%s).\n"%\
                (path_finder.rna_fn()))
        return False
    #
    to_rel_path=path_finder.to_rel_path
    #
    out=file(path_finder.info_idx_fn(),'wt')
    #out.write("SHORT_NAME\t%s\n"%dn.short_name)
    #out.write("SOURCE_FAS\t%s\n"%dn.dna_fn)
    #
    #out.write("#\n")
    #out.write("CDNA_FAS\t%s\n"%path_finder.rna_fn())
    #out.write("CDNA_FAS_FAI\t%s\n"%rna_fai_fn)
    #out.write("CDNA_BOWTIEIDX\t%s\n"%(path_finder.rna_fn()[:-3]))
    #out.write("CDNA_BWAIDX\t%s\n"%(path_finder.rna_fn()))
    ##
    #out.write("#\n")
    #out.write("CDNA_CODING_FAS\t%s\n"%path_finder.rna_coding_fn())
    #out.write("CDNA_CODING_FAS_FAI\t%s\n"%rna_coding_fai_fn)
    #out.write("CDNA_CODING_BOWTIEIDX\t%s\n"%(path_finder.rna_coding_fn()[:-3]))
    #out.write("CDNA_CODING_BWAIDX\t%s\n"%(path_finder.rna_coding_fn()))
    ##
    #out.write("#\n")
    #out.write("NCRNA_FAS\t%s\n"%path_finder.ncrna_fn())
    #out.write("NCRNA_FAS_FAI\t%s\n"%ncrna_fai_fn)
    #out.write("NCRNA_BOWTIEIDX\t%s\n"%(path_finder.ncrna_fn()[:-3]))
    #out.write("NCRNA_BWAIDX\t%s\n"%(path_finder.ncrna_fn()))
    ##
    #out.write("#\n")
    #out.write("CHRO_FAS\t%s\n"%path_finder.chr_fn())
    #out.write("CHRO_FAS_FAI\t%s\n"%to_rel_path(rna_fai_fn))
    #out.write("CHRO_FAS_PATH\t%s\n"%path_finder.chr_path())
    #out.write("CHRO_FAS_SOURCE_PATH\t%s\n"%path_finder.chr_source_path())
    #out.write("CHRO_BOWTIEIDX\t%s\n"%chr_bowtie_fn)
    #out.write("CHRO_BWAIDX\t%s\n"%chr_bwa_fn)
    #out.write("CHRO_STARIDX\t%s\n"%chr_star_fn)
    #out.write("BOWTIEIDX\t%s\n"%chr_bowtie_fn)
    #out.write("BWAIDX\t%s\n"%chr_bwa_fn)
    #out.write("STARIDX\t%s\n"%chr_star_fn)
    #
    out.write("#\n")
    out.write("SEQ\tCHRO:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.chr_fn()))
    out.write("INDEX\tCHRO:FASTA\t%s\n"%to_rel_path(chr_fai_fn))
    out.write("INDEX\tCHRO:BWA\t%s\n"%to_rel_path(chr_bwa_fn))
    out.write("INDEX\tCHRO:BOWTIE\t%s\n"%to_rel_path(chr_bowtie_fn))
    out.write("INDEX\tCHRO:STAR\t%s\n"%to_rel_path(chr_star_fn))
    out.write("INDEX\tCHRO:DICT\t%s\n"%to_rel_path(chr_dict_fn))
    out.write("INDEX\tCHRO:BLAST\t%s\n"%to_rel_path(path_finder.chr_fn()))
    #
    out.write("#\n")
    #out.write("SEQ\tRNA:CDNA:ALL:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.rna_fn()))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:FASTA\t%s\n"%to_rel_path(rna_fai_fn))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:BWA\t%s\n"%to_rel_path(rna_bwa_fn))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:BOWTIE\t%s\n"%to_rel_path(rna_bowtie_fn))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:DICT\t%s\n"%to_rel_path(rna_dict_fn))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:FASTA\t%s\n"%to_rel_path(rna_fai_fn))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:BWA\t%s\n"%to_rel_path(rna_bwa_fn))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:BOWTIE\t%s\n"%to_rel_path(rna_bowtie_fn))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:DICT\t%s\n"%to_rel_path(rna_dict_fn))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:DICT\t%s\n"%to_rel_path(rna_dict_fn))
    out.write("INDEX\tPROTEIN:FULL:BLAST:FORMAT=BLAST\t%s\n"%to_rel_path(pep_blast_fn))
    out.write("INDEX\tPROTEIN:REPR:BLAST:FORMAT=BLAST\t%s\n"%to_rel_path(pep_blast_fn))
    #
    #out.write("#\n")
    #out.write("SEQ\tRNA:CDNA:CODING:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.rna_coding_fn()))
    #out.write("INDEX\tRNA:CDNA:CODING:FASTA\t%s\n"%to_rel_path(rna_coding_fai_fn))
    #out.write("INDEX\tRNA:CDNA:CODING:BWA\t%s\n"%to_rel_path(rna_coding_bwa_fn))
    #out.write("INDEX\tRNA:CDNA:CODING:BOWTIE\t%s\n"%to_rel_path(rna_coding_bowtie_fn))
    #out.write("INDEX\tRNA:CDNA:CODING:DICT\t%s\n"%to_rel_path(rna_coding_dict_fn))
    #
    #if os.path.exists(path_finder.ncrna_fn()) \
    #        and os.path.getsize(path_finder.ncrna_fn())>1000 :
    #    out.write("#\n")
    #    out.write("SEQ\tRNA:NCRNA:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.ncrna_fn()))
    #    out.write("INDEX\tRNA:NCRNA:CDNA:FASTA\t%s\n"%to_rel_path(ncrna_fai_fn))
    #    out.write("INDEX\tRNA:NCRNA:CDNA:BOWTIE\t%s\n"%to_rel_path(ncrna_bowtie_fn))
    #    out.write("INDEX\tRNA:NCRNA:CDNA:BWA\t%s\n"%to_rel_path(ncrna_bwa_fn))
    #    out.write("INDEX\tRNA:NCRNA:CDNA:DICT\t%s\n"%to_rel_path(ncrna_dict_fn))
    #
    out.close()

    return True

def export(path_finder,run) :

    #db=RineDB(path_finder.db_fn())
    #meta_ass_name=db.search_meta("assembly.name")
    #meta_ass_acc=db.search_meta("assembly.accession")
    #meta_tax_id=db.search_meta("species.taxonomy_id")
    #meta_name=db.search_meta("species.display_name")
    #meta_sci_name=db.search_meta("species.scientific_name")
    #meta_gencode=db.search_meta("gencode.version")

    #info_dn=read_info_down(path_finder.info_down_fn())
    #list_fn="%s.list"%(run.genome_fa)
    #dna_s=read_dna_s(list_fn)

    dna_s=[]
    for line in file(run.genome_fa) :
        if line.startswith(">") :
            dna_s.append(line[1:].strip().split()[0])

    to_rel_path=path_finder.to_rel_path

    out=file(path_finder.info(),'wt')
    out.write("NAME\t%s\n"%(path_finder.short_name()))
    out.write("SPECIES\t%s\n"%(run.species))
    out.write("SOURCE\t%s\t%s\n"%(run.source,run.release))
    out.write("PATH\t%s\n"%(path_finder.home()))

    #out.write("#\n")
    #out.write("META\tNAME\t%s\n"%(meta_name))
    #out.write("META\tSCI_NAME\t%s\n"%(meta_sci_name))
    #out.write("META\tNCBI_TAX_ID\t%s\n"%(meta_tax_id))
    #out.write("META\tASSEMBLY_NAME\t%s\n"%(meta_ass_name))
    #out.write("META\tASSEMBLY_ACC\t%s\n"%(meta_ass_acc))
    #if meta_gencode!=None :
    #    out.write("META\tGENCODE\t%s\n"%(meta_gencode))
    #else :
    #    out.write("META\tGENCODE\tNA\n")

    #out.write("#\n")
    #out.write("GENES\tSOURCE\tGTF\t%s\n"%to_rel_path(run.gtf_fn))
    #out.write("GENES\tALL\tGTF\t%s\n"%to_rel_path(path_finder.gtf_fn()))
    #out.write("GENES\tREPR\tGTF\t%s\n"%to_rel_path(path_finder.repr_gtf()))
    #out.write("GENES\tCODING\tGTF\t%s\n"%to_rel_path(path_finder.coding_gtf()))
    #out.write("GENES\tCODING_REPR\tGTF\t%s\n"%to_rel_path(path_finder.coding_repr_gtf()))
    #out.write("GENES\tNC_RNA\tGTF\t%s\n"%\
    #        to_rel_path((path_finder.non_coding_all_gtf())))
    #out.write("GENES\tMT_RNA\tGTF\t%s\n"%to_rel_path(path_finder.mt_rna_gtf()))
    #out.write("GENES\tPT_RNA\tGTF\t%s\n"%to_rel_path(path_finder.pt_rna_gtf()))
    #out.write("GENES\tTnR_RNA\tGTF\t%s\n"%to_rel_path(path_finder.t_r_rna_gtf()))
    #out.write("GENES\tTnRnMTnPT_RNA\tGTF\t%s\n"%\
    #        to_rel_path(path_finder.t_r_mt_pt_rna_gtf()))
    #out.write("GENES\tLINC_RNA\tGTF\t%s\n"%to_rel_path(path_finder.linc_gtf()))
    #out.write("GENES\tMI_RNA\tGTF\t%s\n"%to_rel_path(path_finder.mirna_gtf()))
    #out.write("GENES\tMASK\tGTF\t%s\n"%to_rel_path(path_finder.mask_gtf()))
    #
    out.write("#\n")
    out.write("GENESET\tSOURCE:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn()))
    out.write("GENESET\tALL:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn()))
    out.write("GENESET\tSOURCE:FULL:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn()))
    out.write("GENESET\tCODING:FULL:LEVEL=GENE:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.coding_gtf()))
    out.write("GENESET\tNONCODING:FULL:LEVEL=GENE:FORMAT=GTF\t%s\n"%\
            to_rel_path((path_finder.non_coding_all_gtf())))
    out.write("GENESET\tCODING:FULL:LEVEL=TRAN:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.coding_gtf()))
    out.write("GENESET\tNONCODING:FULL:LEVEL=TRAN:FORMAT=GTF\t%s\n"%\
            to_rel_path((path_finder.non_coding_all_gtf())))
    out.write("GENESET\tRNA:MASK:MT:PT:NONCODING:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.mask_gtf()))
    #
    out.write("#\n")
    out.write("GENESET\tALL:EXONIC:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.gtf_fn(tag='exonic')))
    out.write("GENESET\tNONCODING:EXONIC:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path((path_finder.non_coding_all_gtf(tag='exonic'))))
    out.write("GENESET\tCODING:EXONIC:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.coding_gtf(tag='exonic')))
    out.write("GENESET\tMASK:MT:PT:NONCODING:TOOL=DEXSEQ:FORMAT=GTF\t%s\n"%\
            to_rel_path(path_finder.mask_gtf(tag='exonic')))
    #
    out.write("#\n")
    for line in file(path_finder.info_desc_fn()) :
        if line.startswith("#") :
            continue
        out.write(line)
    out.write("#\n")
    #out.write("INFO\tGENES_TRANS\t%s\n"%to_rel_path(info.ens_gene_tran_info))
    #out.write("GENE\tINFO\t%s\n"%to_rel_path(info.ens_gene_info))
    #out.write("GENE\tFEATURE\t%s\n"%to_rel_path(info.ens_gene_feature))
    #out.write("TRANSCRIPT\tINFO\t%s\n"%to_rel_path(info.ens_tran_info))
    #out.write("TRANSCRIPT\tFEATURE\t%s\n"%to_rel_path(info.ens_tran_feature))
    #out.write("#\n")
    #out.write("PEP_FAS\t%s\n"%path_finder.to_rel_path(info_dn.pep_fn))
    #out.write("PEP_FAS_PATH\t[PATH]/pep\n")
    # CDNA
    #out.write("# =============================================\n")
    #out.write("# INDEX\n")
    #out.write("# ---------------------------------------------\n")
    #out.write("CDNA_FAS\t%s\n"%to_rel_path(path_finder.rna_fn()))
    #out.write("CDNA_FAS_FAI\t%s.fai\n"%to_rel_path(path_finder.rna_fn()))
    #out.write("CDNA_BOWTIEIDX\t%s\n"%to_rel_path(path_finder.rna_fn()[:-3]))
    #out.write("CDNA_BWAIDX\t%s\n"%to_rel_path(path_finder.rna_fn()))
    #
    #out.write("#\n")
    #out.write("CDNA_CODING_FAS\t%s\n"%to_rel_path(path_finder.rna_coding_fn()))
    #out.write("CDNA_CODING_FAS_FAI\t%s.fai\n"%to_rel_path(path_finder.rna_coding_fn()))
    #out.write("CDNA_CODING_BOWTIEIDX\t%s\n"%to_rel_path(path_finder.rna_coding_fn()[:-3]))
    #out.write("CDNA_CODING_BWAIDX\t%s\n"%to_rel_path(path_finder.rna_coding_fn()))
    # NCRNA
    #out.write("#\n")
    #out.write("NCRNA_FAS\t%s\n"%to_rel_path(path_finder.ncrna_fn()))
    #out.write("NCRNA_FAS_FAI\t%s.fai\n"%to_rel_path(path_finder.ncrna_fn()))
    #out.write("NCRNA_BOWTIEIDX\t%s\n"%to_rel_path(path_finder.ncrna_fn()[:-3]))
    #out.write("NCRNA_BWAIDX\t%s\n"%to_rel_path(path_finder.ncrna_fn()))
    # Chromosome
    #out.write("#\n")
    #out.write("CHRO_SOURCE_FAS\t%s\n"%to_rel_path(info.source_fas))
    #out.write("CHRO_FAS\t%s\n"%to_rel_path(path_finder.chr_fn()))
    #out.write("CHRO_FAS_FAI\t%s.fai\n"%to_rel_path(path_finder.chr_fn()))
    #out.write("CHRO_FAS_PATH\t%s\n"%to_rel_path(path_finder.chr_path()))
    #out.write("CHRO_BOWTIEIDX\t%s\n"%to_rel_path(path_finder.chr_fn()[:-3]))
    #out.write("CHRO_BWAIDX\t%s\n"%to_rel_path(path_finder.chr_fn()))
    #out.write("CHRO_STARIDX\t%s\n"%to_rel_path(path_finder.staridx_path()))
    #
    if os.path.exists(path_finder.info_idx_fn()) :
        for line in file(path_finder.info_idx_fn()) :
            unit_s=line.strip().split()
            tag=unit_s[0]
            if tag not in ('BOWTIEIDX','BWAIDX','STARIDX') :
                continue
            out.write("%s\t%s\n"%(tag,to_rel_path(unit_s[1])))
    #===========================================================================
    # Add New format
    # --------------------------------------------------------------------------
    out.write("#\n")
    out.write("SEQ\tCHRO:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.chr_fn()))
    out.write("SEQ\tRNA:FULL:CDNA:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.rna_fn()))
    #out.write("SEQ\tRNA:REPR:CDNA:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.rna_fn()))
    out.write("SEQ\tPROTEIN:FULL:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.pep_fn()))
    out.write("SEQ\tPROTEIN:REPR:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(path_finder.pep_fn()))
    #out.write("SEQ\tPROTEIN:FULL:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(info_dn.pep_fn))
    #out.write("SEQ\tPROTEIN:REPR:FASTA:FORMAT=FASTA\t%s\n"%to_rel_path(info_dn.pep_fn))
    # --------------------------------------------------------------------------
    #out.write("#\n")
    #if os.path.exists(path_finder.info_idx_fn()) :
    #    for line in file(path_finder.info_idx_fn()) :
    #        unit_s=line.strip().split()
    #        if unit_s[0]=='INDEX' :
    #            out.write(line)
    out.write("#\n")
    out.write("INDEX\tCHRO:FASTA\t%s\n"%to_rel_path(build_fai_idx(path_finder.chr_fn())))
    out.write("INDEX\tCHRO:BWA\t%s\n"%to_rel_path(build_bwa_idx(path_finder.chr_fn())))
    out.write("INDEX\tCHRO:BOWTIE\t%s\n"%to_rel_path(build_bowtie_idx(path_finder.chr_fn())))
    out.write("INDEX\tCHRO:STAR\t%s\n"%to_rel_path(build_star_idx(path_finder,path_finder.chr_fn())))
    out.write("INDEX\tCHRO:DICT\t%s\n"%to_rel_path(build_dict_idx(path_finder.chr_fn())))
    out.write("INDEX\tCHRO:BLAST\t%s\n"%to_rel_path(build_blast_n_idx(path_finder.chr_fn())))
    out.write("#\n")
    out.write("INDEX\tRNA:FASTA\t%s\n"%to_rel_path(build_fai_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:FASTA\t%s\n"%to_rel_path(build_fai_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:BWA\t%s\n"%to_rel_path(build_bwa_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:BOWTIE\t%s\n"%to_rel_path(build_bowtie_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:FULL:ALL:DICT\t%s\n"%to_rel_path(build_dict_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:FASTA\t%s\n"%to_rel_path(build_fai_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:BWA\t%s\n"%to_rel_path(build_bwa_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:BOWTIE\t%s\n"%to_rel_path(build_bowtie_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:DICT\t%s\n"%to_rel_path(build_dict_idx(path_finder.rna_fn())))
    out.write("INDEX\tRNA:CDNA:REPR:CODING:BLAST\t%s\n"%to_rel_path(build_blast_n_idx(path_finder.rna_fn())))
    out.write("INDEX\tPROTEIN:FULL:BLAST:FORMAT=BLAST\t%s\n"%to_rel_path(build_blast_p_idx(path_finder.pep_fn())))
    out.write("INDEX\tPROTEIN:REPR:BLAST:FORMAT=BLAST\t%s\n"%to_rel_path(build_blast_p_idx(path_finder.pep_fn())))
    # --------------------------------------------------------------------------
    #out.write("#\n")
    #out.write("#INDEX\tCHRO:CHIMERASCAN\t[PATH]/chimeraidx\n")
    #out.write("#INDEX\tCHRO:SNPiR:JUN95\t[PATH]/SNPiR/%s_jun95.chr.fa\n"%\
    #        (path_finder.short_name()))
    # --------------------------------------------------------------------------
    #out.write("#\n")
    #out.write("#VCF\tCOSMIC:VERSION=66\t[PATH]/vcf_s/b37_cosmic_V66_250713.vcf\n")
    #out.write("#VCF\tDBSNP:VERSION=132\t[PATH]/vcf_s/dbsnp_132_b37.leftAligned.vcf\n")
    #out.write("#VCF\tINDEL:GOLD:1000G\t[PATH]/vcf_s/1000G_phase1.indels.b37.vcf\n")
    # --------------------------------------------------------------------------
    out.write("#\n")
    out.write("FEATURE\tGENE:INFO\t%s\n"%to_rel_path(path_finder.gene_info()))
    out.write("FEATURE\tGENE:FEATURE\t%s\n"%to_rel_path(path_finder.gene_feature()))
    out.write("FEATURE\tTRANSCRIPT:INFO\t%s\n"%to_rel_path(path_finder.tran_info()))
    out.write("FEATURE\tTRANSCRIPT:FEATURE\t%s\n"%to_rel_path(path_finder.tran_feature()))
    out.write("FEATURE\tGENE_TRAN:MAP\t%s\n"%to_rel_path(path_finder.gene_tran_map()))
    #
    out.write("#FEATURE\tREPEAT:BED\t[PATH]/%s.repeat.bed\n"%(path_finder.short_name()))
    #out.write("#FEATURE\tGENE:MERGED\t[PATH]/SNPiR/hg19.annotation.new.txt\n")
    #out.write("#FEATURE\tEDIT:BED\t[PATH]/SNPiR/Human_AG_all_hg19.bed\n")
    #===========================================================================
    out.write("#\n")
    out.write("CHRO_MAJOR\t%s\n"%(",".join(dna_s)))
    out.write("#\n")
    out.write("CHRO_LIST\n")
    out.write("\n".join(dna_s))
    out.write("\n")
    out.close()

def interpret(argv) :
    optlist, args = getopt.getopt(argv, 'g:t:s:r:h:v:p:')
    run=Run()
    for param,val in optlist :
        if param=='-g' :
            run.genome_fa=val
        elif param=='-t' :
            run.gtf_fn=val
        #
        elif param=='-r' :
            run.source=val
        elif param=='-s' :
            run.species=val.replace(" ","_")
            run.name_s=val.replace("_"," ").split()
        elif param=='-v' :
            run.release=val
        #
        elif param=='-h' :
            run.ref_home=val
        elif param=='-p' :
            run.pep_fa=val
    return run
 
def usage(out) :
    out.write("""Rine Local Auto-building System
    Version 1.0

    python %s [Options]

    -g [GenomeFasta]
    -t [Ensembl GTF]

    -r [Source]
    -s [Species]
    -v [Release]

    -h [Ref Home]
    -p [ProteinFasta]

    """%__file__)
    out.write("\n")

def main() :
    if len(sys.argv)<4 :
        usage(sys.stdout)
        return
    #
    run=interpret(sys.argv[1:])
    #
    if run.species==None :
        print "ERROR: No SPECIES_INFO.\n"
        return False
    #
    path_finder=PathFinder(run.species,run.release,run.name_s)
    path_finder._ref_home=run.ref_home
    path_finder._source=run.source
    sys.stdout.write("HOME\t%s\n"%path_finder.home())
    #
    if not prepare_path(path_finder) :
        print "Failed to PREPARE PATH."
        return False
    # Rebuild CHR_FN
    if not rebuild_chr_fn(path_finder,run) :
        sys.stderr.write("ERROR: failed to rebuild CHR_FN")
        return
    if not rebuild_pep_fn(path_finder,run) :
        sys.stderr.write("ERROR: failed to rebuild PEP_FN")
        return
    if not rebuild_chr_gene_set(path_finder,run) :
        return 
    # Select GENE_ID
    if not prepare_gene_id(path_finder):
        sys.stderr.write("ERROR: failed to select GENE_ID.")
        return
    if not prepare_coding_gene_id(path_finder):
        sys.stderr.write("ERROR: failed to select CODING_GENE_ID.")
        return
    # Prepare GENESET
    if not build_gene_set(path_finder) :
        return False
    if not build_exonic_gene_set(path_finder) :
        return False
    # Rebuild CDNA_FN
    if not rebuild_cdna_fn(path_finder) :
        sys.stdout.write("ERROR: failed to rebuild CDNA_FN\n")
        return
    # Build index
    if not build_index(path_finder) :
        sys.stdout.write("ERROR: Failed to build index.\n")
        return
    if not os.path.exists(path_finder.info_desc_fn()) :
        build_desc(path_finder,run.species)
    #
    if not build_mapping(path_finder,run.species) :
        return 
    #
    #if not os.path.exists(path_finder.info_repeat_fn()) :
    #    build_repeat(path_finder,run.species)
    #
    export(path_finder,run)
    # prepare_repr_gene_set(path_finder,info)
    #

if __name__=='__main__' :
    main()

