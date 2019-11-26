
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fasta_fn_dic = {
    'floweringtimeA':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/floweringTime_A.faa',
    'arabidopsis':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Arabidopsis_thaliana.fasta',
    'carica':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Carica_papaya.fasta',
    'korea':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Korea_Pear.fasta',
    'malus':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Malus_domestica.fasta',
    'medicago':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Medicago_truncatula.fasta',
    'populus':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Populus_trichocarpa.fasta',
    'pyrus':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Pyrus_bretschneideri.fasta',
    'vitis':'/BiO/sgpark/Projects/RDA_Pear/AdvancedAnalysis/OrthoMCL/peptides/Vitis_vinifera.fasta'}

flowering_time_a_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/A_Gene/A_vs_Vitis.xls'
    }
flowering_time_b_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/B_Gene/B_vs_Vitis.xls'
    }
flowering_time_c_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/C_Gene/C_vs_Vitis.xls'
    }
flowering_time_d_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/D_Gene/D_vs_Vitis.xls'
    }
flowering_time_e_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/flowering_time/E_Gene/E_vs_Vitis.xls'
    }
lignin_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/Lignin/Lignin_vs_Vitis.xls'
    }
sugar_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/sugar_patyway/sugar_vs_Vitis.xls'
    }
tf_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/Transcription_factor/TF_vs_Vitis.xls'
    }
volatile_a_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/A_Gene/A_vs_Vitis.xls'
    }
volatile_b_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/B_Gene/B_vs_Vitis.xls'
    }
volatile_c_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/C_Gene/C_vs_Vitis.xls'
    }
volatile_d_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/D_Gene/D_vs_Vitis.xls'
    }
volatile_e_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/E_Gene/E_vs_Vitis.xls'
    }
volatile_f_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/F_Gene/F_vs_Vitis.xls'
    }
volatile_g_blast_fn_dic = {
    'arabidopsis':'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Arabidopsis.xls',
    'carica'     :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Carica.xls',
    'korea'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Korea.xls',
    'malus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Malus.xls',
    'medicago'   :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Medicago.xls',
    'populus'    :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Populus.xls',
    'pyrus'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Pyrus.xls',
    'vitis'      :'/BiO/BioPeople/boram/Projects/pear_peach/volatile/G_Gene/G_vs_Vitis.xls'
    }












def load_blast(_blast_fn_dic):
    blast_fn_dic = dict()
    for queryname, xls_fn in _blast_fn_dic.items():
        if os.path.isfile(xls_fn):
            blast_fn_dic[queryname] = xls_fn
        else:
            print('Could not found : {0}'.format(xls_fn))
            sys.exit()
    return blast_fn_dic


class OrthologousFinder():
    def __init__(self, targetname, blast_fn_dic, fasta_fn_dic):
        self._dic = dict()
        self.targetname = targetname
        self.blast_screener(blast_fn_dic)
        #
        self.fasta_fn_dic = fasta_fn_dic

    def blast_screener(self, blast_fn_dic):
        self.querynames = list()
        for queryname, xls_fn in blast_fn_dic.items():
            self.querynames.append(queryname)
            self.blast_parser(queryname, xls_fn)
        self.fill_dic(self.querynames)

    def blast_parser(self, queryname, xls_fn):
        for line in open(xls_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['queryNum']:
                #['queryNum', 'hit_num', 'hsp_num',
                # 'queryID', 'e-value', 'bit-score', 'queryLen',
                # 'targetLen', 'targetAcc', 'targetDesc',
                # 'queryCov', 'targetCov', 'alignLen', 'alignSim',
                # 'queryStr', 'queryEnd', 'targetStr', 'targetEnd',
                # 'gapInfo', 'frameInfo', 
                # 'querySeq', 'comparison', 'targetSeq', 'targetDB']
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            q_id = None
            if queryname in ['arabidopsis']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['carica']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['korea']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['malus']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['medicago']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['populus']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['pyrus']:
                q_id = items[idx_dic['queryID']].split()[0]
            elif queryname in ['vitis']:
                q_id = items[idx_dic['queryID']].split()[0]
            #
            t_id = None
            if self.targetname in ['floweringtimeA']:
                t_id = items[idx_dic['targetDesc']].split()[0]
            else:
                t_id = items[idx_dic['targetDesc']].split()[0]
            #
            evalue = float(items[idx_dic['e-value']])
            bit = float(items[idx_dic['bit-score']])
            queryCov = float(items[idx_dic['queryCov']])
            targetCov = float(items[idx_dic['targetCov']])
            alignLen = int(items[idx_dic['alignLen']])
            alignSim = float(items[idx_dic['alignSim']])
            #
            q_s = int(items[idx_dic['queryStr']])
            q_e = int(items[idx_dic['queryEnd']])
            t_s = int(items[idx_dic['targetStr']])
            t_e = int(items[idx_dic['targetEnd']])
            #
            #print(q_id, t_id, evalue, bit, queryCov, targetCov, alignLen, alignSim)
            if self.filtering(evalue, queryCov, targetCov, alignLen, alignSim):
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('q_s', q_s)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('q_e', q_e)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('t_s', t_s)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('t_e', t_e)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('alignLen', alignLen)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('alignSim', alignSim)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('e-value', evalue)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('bit-score', bit)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('queryCov', queryCov)
                self._dic.setdefault(t_id, {}).setdefault(queryname, {}).setdefault(q_id, {}).setdefault('targetCov', targetCov)
            #
        #

    def filtering(self, evalue, queryCov, targetCov, alignLen, alignSim):
        if float(evalue) > 1e-5:
            return 0
        if float(queryCov) < 0.0:
            return 0
        if float(targetCov) < 50.0:
            return 0
        if float(alignLen) < 0:
            return 0
        if float(alignSim) < 0.0:
            return 0
        return 1

    def fill_dic(self, querynames):
        for t_id in self._dic:
            for queryname in self.querynames:
                self._dic[t_id].setdefault(queryname, {})

    def write_table(self):
        outfn = '{0}.orthologous.xls'.format(self.targetname)
        outfh = open(outfn, 'w')
        headers = ['target_id']
        headers.append('count.not0')
        headers.append('sum.num')
        for queryname in self.querynames:
            headers.append('num.{0}'.format(queryname))
        for queryname in self.querynames:
            headers.append('members.{0}'.format(queryname))
        outfh.write('{0}\n'.format('\t'.join(headers)))
        for t_id, sub_dic in self._dic.items():
            items = [t_id]
            _count_not0 = 0
            _sum_num = 0
            for queryname in self.querynames:
                if not len(sub_dic[queryname]) in [0]:
                    _count_not0 += 1
                _sum_num += len(sub_dic[queryname])
            items.append(str(_count_not0))
            items.append(str(_sum_num))
            for queryname in self.querynames:
                items.append(str(len(sub_dic[queryname])))
            for queryname in self.querynames:
                items.append(','.join(sub_dic[queryname]))
            outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()
        print(outfn)
        return outfn

    def write_seq(self, target_id):
        if not target_id:
            return 0
        #
        best_dic = dict()
        for queryname, query_id_dic in self._dic[target_id].items():
            if args.best_selection_method in ['alignLen','alignSim','bit-score','queryCov','targetCov']:
                for query_id, align_dic in sorted(query_id_dic.items(), key=lambda (k,v): (v[args.best_selection_method],k), reverse=True):
                    q_s = align_dic['q_s']
                    q_e = align_dic['q_e']
                    t_s = align_dic['t_s']
                    t_e = align_dic['t_e']
                    value = align_dic[args.best_selection_method]
                    print(queryname, query_id, q_s, q_e, t_s, t_e, value)
                    best_dic.setdefault(queryname, [query_id, q_s, q_e, t_s, t_e])
            elif args.best_selection_method in ['e-value']:
                for query_id, align_dic in sorted(query_id_dic.items(), key=lambda (k,v): (v[args.best_selection_method],k)):
                    q_s = align_dic['q_s']
                    q_e = align_dic['q_e']
                    t_s = align_dic['t_s']
                    t_e = align_dic['t_e']
                    value = align_dic[args.best_selection_method]
                    print(queryname, query_id, q_s, q_e, t_s, t_e, value)
                    best_dic.setdefault(queryname, [query_id, q_s, q_e, t_s, t_e])
        #
        outfn = '{0}.{1}.{2}.orthologous.fasta'.format(self.targetname, target_id, args.best_selection_method)
        outfh = open(outfn, 'w')
        for queryname, query_info_s in best_dic.items():
            print(queryname, query_info_s)
            query_id = query_info_s[0]
            q_s = query_info_s[1]
            q_e = query_info_s[2]
            t_s = query_info_s[3]
            t_e = query_info_s[4]
            for record in SeqIO.parse(open(self.fasta_fn_dic[queryname]), 'fasta'):
                if query_id in record.id:
                    query_region_id = '{0}({1}:{2}-{3})'.format(queryname, query_id, q_s, q_e)
                    query_region_seq = record.seq[q_s-1:q_e-1+1]
                    print(queryname, query_region_id, len(query_region_seq), query_region_seq)
                    region_record = SeqRecord(query_region_seq, id=query_region_id, description='')
                    SeqIO.write(region_record, outfh, 'fasta')
                else:
                    pass
            #
        outfh.close()
        print(outfn)
        return outfn


def main(args):
    if args.targetname in ['floweringtimeA']:
        blast_fn_dic = load_blast(flowering_time_a_blast_fn_dic)
    elif args.targetname in ['floweringtimeB']:
        blast_fn_dic = load_blast(flowering_time_b_blast_fn_dic)
    elif args.targetname in ['floweringtimeC']:
        blast_fn_dic = load_blast(flowering_time_c_blast_fn_dic)
    elif args.targetname in ['floweringtimeD']:
        blast_fn_dic = load_blast(flowering_time_d_blast_fn_dic)
    elif args.targetname in ['floweringtimeE']:
        blast_fn_dic = load_blast(flowering_time_e_blast_fn_dic)
    elif args.targetname in ['lignin']:
        blast_fn_dic = load_blast(lignin_blast_fn_dic)
    elif args.targetname in ['sugar']:
        blast_fn_dic = load_blast(sugar_blast_fn_dic)
    elif args.targetname in ['transcriptfactor']:
        blast_fn_dic = load_blast(tf_blast_fn_dic)
    elif args.targetname in ['volatileA']:
        blast_fn_dic = load_blast(volatile_a_blast_fn_dic)
    elif args.targetname in ['volatileB']:
        blast_fn_dic = load_blast(volatile_b_blast_fn_dic)
    elif args.targetname in ['volatileC']:
        blast_fn_dic = load_blast(volatile_c_blast_fn_dic)
    elif args.targetname in ['volatileD']:
        blast_fn_dic = load_blast(volatile_d_blast_fn_dic)
    elif args.targetname in ['volatileE']:
        blast_fn_dic = load_blast(volatile_e_blast_fn_dic)
    elif args.targetname in ['volatileF']:
        blast_fn_dic = load_blast(volatile_f_blast_fn_dic)
    elif args.targetname in ['volatileG']:
        blast_fn_dic = load_blast(volatile_g_blast_fn_dic)

    orthologous = OrthologousFinder(args.targetname, blast_fn_dic, fasta_fn_dic)
    orthologous_table = orthologous.write_table()
    orthologous_seq = orthologous.write_seq(args.targetid)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--best_selection_method', choices=('alignLen','e-value','alignSim','bit-score','queryCov','targetCov'),
            default='alignLen')
    parser.add_argument('--targetname', choices=('floweringtimeA','floweringtimeB','floweringtimeC', 'floweringtimeD','floweringtimeE',
                                                 'lignin', 'sugar','transcriptfactor',
                                                 'volatileA', 'volatileB', 'volatileC', 'volatileD', 'volatileE', 'volatileF', 'volatileG'),
            default='floweringtimeB')
    parser.add_argument('--targetid', default='AT4G24210.1')
    args = parser.parse_args()
    main(args)
